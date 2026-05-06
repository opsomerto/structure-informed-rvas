"""
Phase 2 of the Harmonic Mean P-value (HMP) pipeline for structure-informed RVAS.

Reads the master p_values.h5 produced by q_scan_test_hmp.build_hmp_master_h5
and computes FDR and FWER corrections jointly over all (M × N) entries
(M neighborhoods × N clusters).
"""

import os
import numpy as np
import pandas as pd
import h5py
from logger_config import get_logger

logger = get_logger(__name__)

_NEUTRAL = {'tstat': 0.0, 'hill': 0.0, 'pval': 1.0, 'hillp': 1.0}


def _prepare_fdr_filters(df_fdr_filter):
    """Prepare optional filtering criteria."""
    if df_fdr_filter is None:
        return None, None

    uniprot_filter_list = np.unique(df_fdr_filter['uniprot_id'])
    aa_pos_filters = None

    if 'aa_pos' in df_fdr_filter.columns:
        aa_pos_filters = {}
        for uid in uniprot_filter_list:
            aa_pos_filters[uid] = set(
                df_fdr_filter.loc[df_fdr_filter['uniprot_id'] == uid, 'aa_pos'].values
            )

    return uniprot_filter_list, aa_pos_filters


def _load_all_stats_hmp(results_dir, uniprot_filter_list, aa_pos_filters, pval_file):
    """
    Load observed and null statistics from the master HMP HDF5 file.

    Entries are indexed by the 'entry_ids' dataset written by build_hmp_master_h5;
    each key has the form '{uniprot_id}.{cluster}'.  The aa_pos filter uses
    np.isin to correctly handle non-contiguous positions (e.g. from --select-nbhds).

    Returns:
        df_stats:        DataFrame [uniprot_id, cluster, aa_pos, n_a_tstat,
                                    mean_betahat, n_betahat] sorted by n_a_tstat
        null_stats_dict: {entry_key: ndarray (n_rows, n_sims)}
        entry_keys:      list of keys in null_stats_dict (same order)
        n_sims:          int
    """
    to_concat = []
    null_stats_dict = {}
    n_sims = None

    h5_path = os.path.join(results_dir, pval_file)
    with h5py.File(h5_path, 'r') as fid:
        if 'entry_ids' not in fid:
            raise ValueError(
                "Master HDF5 file missing 'entry_ids' index. "
                "Was it created by build_hmp_master_h5?"
            )
        entry_ids = [e.decode('ascii') for e in fid['entry_ids'][:]]
        logger.info(f'Reading {len(entry_ids)} entries from {h5_path}')

        for eid in entry_ids:
            dot_idx = eid.index('.')
            uniprot_id = eid[:dot_idx]
            cluster = eid[dot_idx + 1:]

            if uniprot_filter_list is not None and uniprot_id not in uniprot_filter_list:
                continue

            stat_data = fid[eid][:, 0]                    # (n_res,)
            beta_info = fid[f'{eid}_mean_beta'][:]         # (n_res, 2)
            null_stat = fid[f'{eid}_null_n_a_tstat'][:]   # (n_res, n_sims)

            aa_pos_key = f'{eid}_aa_pos'
            aa_pos = (
                fid[aa_pos_key][:].flatten()
                if aa_pos_key in fid
                else np.arange(1, len(stat_data) + 1)
            )

            # Apply aa_pos filter using np.isin (robust to non-contiguous positions)
            if aa_pos_filters is not None and uniprot_id in aa_pos_filters:
                mask = np.isin(aa_pos, list(aa_pos_filters[uniprot_id]))
                stat_data = stat_data[mask]
                beta_info = beta_info[mask]
                null_stat = null_stat[mask, :]
                aa_pos = aa_pos[mask]

            if len(stat_data) == 0:
                continue

            if n_sims is None:
                n_sims = null_stat.shape[1]

            to_concat.append(pd.DataFrame({
                'uniprot_id': uniprot_id,
                'cluster': cluster,
                'aa_pos': aa_pos,
                'n_a_tstat': stat_data,
                'mean_betahat': beta_info[:, 0],
                'n_betahat': beta_info[:, 1],
            }))
            null_stats_dict[eid] = null_stat

    if not to_concat:
        raise ValueError("No entries found for FDR computation. Check filters and input data.")

    logger.info('Concatenating and sorting by n_a_tstat')
    df_stats = pd.concat(to_concat, ignore_index=True)
    df_stats = df_stats.sort_values('n_a_tstat').reset_index(drop=True)

    return df_stats, null_stats_dict, list(null_stats_dict.keys()), n_sims


def _hmp_compute_false_discoveries(df_stats, null_stats_dict, entry_keys, n_sims,
                                    large_threshold):
    """
    Count expected false discoveries from the pooled null distribution.

    Uses strict < comparison (large_threshold = neutral value) so all
    non-neutral observed entries enter the FDR pool.
    """
    logger.info('Computing false discoveries for HMP FDR')
    mask = df_stats['n_a_tstat'] < large_threshold

    null_pool = []
    for i, key in enumerate(entry_keys):
        if len(entry_keys) > 100 and i % 100 == 0:
            logger.debug(f'Pooling null stats: entry {i}/{len(entry_keys)}')
        ns = null_stats_dict[key].flatten()
        null_pool.extend(ns[ns < large_threshold])

    null_pool = np.sort(np.array(null_pool))

    false_discoveries = np.empty(len(df_stats))
    if np.any(mask):
        false_discoveries[mask] = (
            np.searchsorted(null_pool, df_stats['n_a_tstat'][mask], side='right')
            / n_sims
        )
    if np.any(~mask):
        false_discoveries[~mask] = float(df_stats.shape[0])

    return false_discoveries


def _hmp_compute_fwer(df_stats, null_stats_dict, entry_keys, n_sims, neutral,
                       chunk_size=50):
    """
    Compute empirical FWER via the minimum-statistic method.

    For each simulation the minimum statistic across all entries is tracked;
    FWER for an observed value x is the fraction of simulations where that
    minimum is <= x.

    Args:
        neutral: neutral value (used to initialize the running minimum so that
                 it is only updated by entries with signal)
    """
    logger.info('Computing empirical FWER for HMP')
    min_stats_per_sim = np.full(n_sims, neutral)

    for chunk_start in range(0, len(entry_keys), chunk_size):
        chunk_end = min(chunk_start + chunk_size, len(entry_keys))
        chunk_keys = entry_keys[chunk_start:chunk_end]
        null_chunk = np.vstack([null_stats_dict[k] for k in chunk_keys])
        min_stats_per_sim = np.minimum(min_stats_per_sim, np.min(null_chunk, axis=0))

    p_obs = df_stats['n_a_tstat'].values.reshape(-1, 1)
    fwer = np.mean(min_stats_per_sim.reshape(1, -1) <= p_obs, axis=1)
    return fwer


def _hmp_apply_corrections(df_stats, false_discoveries, fwer):
    """Apply FDR and FWER corrections."""
    logger.info('Applying FDR and FWER corrections')

    df_stats = df_stats.copy()
    df_stats['false_discoveries_avg'] = false_discoveries
    df_stats['fdr'] = [x / (i + 1) for i, x in enumerate(false_discoveries)]
    df_stats['fdr'] = df_stats['fdr'][::-1].cummin()[::-1]
    df_stats['fwer'] = fwer

    return df_stats[['uniprot_id', 'cluster', 'aa_pos', 'n_a_tstat',
                      'fdr', 'fwer', 'mean_betahat', 'n_betahat']]


def summarize_results_hmp(df_results, fdr_cutoff, fwer_cutoff=0.05):
    """Log a brief summary of significant hits."""
    top_hits = df_results.loc[
        df_results.groupby(['uniprot_id', 'cluster'])['fdr'].idxmin()
    ]
    sig_fdr = top_hits[top_hits['fdr'] < fdr_cutoff].sort_values('n_a_tstat')
    sig_fwer = top_hits[top_hits['fwer'] < fwer_cutoff].sort_values('n_a_tstat')

    logger.info('')
    logger.info(
        f'{len(sig_fdr)} / {len(top_hits)} protein-cluster pairs '
        f'significant at FDR < {fdr_cutoff}'
    )
    logger.info(
        f'{len(sig_fwer)} / {len(top_hits)} protein-cluster pairs '
        f'significant at FWER < {fwer_cutoff}'
    )
    if len(sig_fdr) > 0:
        logger.info(f'Top 20 FDR hits:\n{sig_fdr.head(20).to_string()}')
    if len(sig_fwer) > 0:
        logger.info(f'Top FWER hits:\n{sig_fwer.to_string()}')


def q_compute_fdr_hmp(results_dir, fdr_cutoff, df_fdr_filter, reference_dir, pval_file,
                       stat_method='hill'):
    """
    Compute FDR and FWER corrections over the HMP master h5.

    FDR is computed jointly across all (M × N) neighborhood-cluster entries
    (M neighborhoods, N clusters).  The large_threshold equals the neutral
    value so all non-neutral entries enter the FDR pool.

    Args:
        results_dir:    directory containing pval_file
        fdr_cutoff:     FDR threshold for significance reporting
        df_fdr_filter:  optional DataFrame [uniprot_id (, aa_pos)] to restrict entries
        reference_dir:  directory with reference files (for gene name annotation)
        pval_file:      HDF5 filename (typically 'p_values.h5')
        stat_method:    'tstat', 'hill', or 'pval' — determines large_threshold

    Returns:
        DataFrame [uniprot_id, cluster, aa_pos, n_a_tstat, fdr, fwer,
                   mean_betahat, n_betahat (, gene_name)]
    """
    logger.info('HMP Phase 2: computing FDR and FWER')

    neutral = _NEUTRAL.get(stat_method, 0.0)
    large_threshold = neutral  # strict < comparison in false_discoveries

    uniprot_filter_list, aa_pos_filters = _prepare_fdr_filters(df_fdr_filter)

    df_stats, null_stats_dict, entry_keys, n_sims = _load_all_stats_hmp(
        results_dir, uniprot_filter_list, aa_pos_filters, pval_file
    )

    false_discoveries = _hmp_compute_false_discoveries(
        df_stats, null_stats_dict, entry_keys, n_sims, large_threshold
    )

    fwer = _hmp_compute_fwer(
        df_stats, null_stats_dict, entry_keys, n_sims, neutral
    )

    df_results = _hmp_apply_corrections(df_stats, false_discoveries, fwer)

    # Optional gene name annotation
    guide_path = os.path.join(reference_dir, 'protein_sequence_guide.tsv')
    if os.path.exists(guide_path):
        df_gene = pd.read_csv(guide_path, sep='\t')
        df_results = df_results.merge(
            df_gene[['gene_name', 'uniprot_id']], how='left', on='uniprot_id'
        )
    else:
        logger.warning(
            f'protein_sequence_guide.tsv not found at {guide_path}; '
            'gene_name column will not be added'
        )

    summarize_results_hmp(df_results, fdr_cutoff)
    return df_results
