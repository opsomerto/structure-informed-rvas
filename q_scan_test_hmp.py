"""
Phase 1 of the Harmonic Mean P-value (HMP) pipeline for structure-informed RVAS.

Streams K per-trait p_values.h5 files (produced by q_scan_test), computes the
harmonic mean statistic per cluster across traits, and writes a master
p_values.h5 file consumed by q_empirical_fdr_hmp.py for FDR correction.

Entry keys in the master h5 have the form '{uniprot_id}.{cluster}' (dot
separator; UniProt accessions never contain dots).  An 'entry_ids' dataset
indexes all keys so the FDR reader can enumerate them without inspecting key
names.
"""

import os
import numpy as np
import pandas as pd
import h5py
from utils import write_dataset
from logger_config import get_logger

logger = get_logger(__name__)

_NEUTRAL = {'tstat': 0.0, 'hill': 0.0, 'pval': 1.0, 'hillp': 1.0}


def _h5_path(trait):
    """Return the per-trait p_values.h5 path (relative to current working directory)."""
    return f'q3dnt_results/ukbb_{trait}_pval_all-nbhd_gp250506/p_values.h5'


def _build_cluster(traits, cluster, master_path, neutral, n_sims_expected=None,
                   min_variants=10):
    """
    Stream all trait files for one cluster and accumulate harmonic mean statistics.

    A neighborhood m in trait k is considered *valid* if it has at least
    min_variants variants in the neighborhood (read from the _mean_beta dataset).
    This includes neighborhoods whose observed statistic is neutral (no signal but
    enough data to test), which is the correct behaviour: a valid null result
    should count in the denominator of the harmonic mean.

    For each (protein p, neighborhood m) the HM observed statistic is:

        HM_obs[m]  = K_valid[m] / sum_{k: n_betahat[k,m] >= min_variants} (1 / obs[k,m])

    Note on neutral statistics and method choice:
    - pval/hillp (neutral=1.0): a valid-but-neutral entry contributes 1/1.0 = 1.0 to
      the sum, pulling the HM toward 1.0 (non-significant). Well-behaved.
    - hill/tstat (neutral=0.0): a valid-but-neutral entry contributes 1/0 = inf,
      which after clamping makes HM_obs collapse to 0 (neutral) for that neighborhood.
      This is overly conservative; pval/hillp methods are preferred for the HMP pipeline.

    Args:
        traits:           list of trait strings for this cluster
        cluster:          cluster label string (must not contain a dot)
        master_path:      absolute path to the master HDF5 file (opened with 'a')
        neutral:          neutral value (0.0 for hill/tstat, 1.0 for hillp/pval)
        n_sims_expected:  expected null simulation count (None = infer from first file)
        min_variants:     minimum number of variants in the neighborhood for a valid test

    Returns:
        (entry_ids_written, n_sims)
    """
    S_obs = {}          # p -> (n_res,)          sum of 1/obs for valid traits
    S_null = {}         # p -> (n_res, n_sims)   sum of 1/null for valid traits
    count_valid = {}    # p -> (n_res,) int32     number of valid traits per nbhd
    aa_pos_store = {}   # p -> (n_res,) int       amino acid positions
    n_sims = n_sims_expected

    for k_idx, trait in enumerate(traits):
        h5path = _h5_path(trait)
        if not os.path.exists(h5path):
            logger.warning(f'Trait file not found, skipping: {h5path}')
            continue

        logger.debug(f'  Trait {k_idx + 1}/{len(traits)}: {trait}')
        with h5py.File(h5path, 'r') as fid:
            protein_ids = [k for k in fid.keys() if '_' not in k]
            for p in protein_ids:
                null_key = f'{p}_null_n_a_tstat'
                if null_key not in fid:
                    continue

                obs_stat = fid[p][:, 0]       # (n_res,)
                null_stat = fid[null_key][:]  # (n_res, n_sims)

                if n_sims is None:
                    n_sims = null_stat.shape[1]

                aa_pos_key = f'{p}_aa_pos'
                aa_pos = (
                    fid[aa_pos_key][:].flatten()
                    if aa_pos_key in fid
                    else np.arange(1, len(obs_stat) + 1)
                )

                # Validity: neighborhoods with sufficient variant count.
                # Uses n_betahat from the _mean_beta dataset (col 1).
                # This includes valid-but-neutral neighborhoods (enough variants
                # but no signal), unlike the old obs_stat != neutral criterion.
                mean_beta_key = f'{p}_mean_beta'
                if mean_beta_key in fid:
                    n_betahat = fid[mean_beta_key][:, 1]  # (n_res,)
                    vmask = n_betahat >= min_variants      # (n_res,)
                else:
                    # Fallback if _mean_beta is absent: use old criterion
                    vmask = obs_stat != neutral

                # Initialize accumulators on first encounter for this protein
                if p not in S_obs:
                    S_obs[p] = np.zeros(len(obs_stat))
                    S_null[p] = np.zeros((len(obs_stat), n_sims))
                    count_valid[p] = np.zeros(len(obs_stat), dtype=np.int32)
                    aa_pos_store[p] = aa_pos
                elif len(obs_stat) != len(S_obs[p]):
                    logger.warning(
                        f'{p}: shape mismatch for trait {trait} '
                        f'(initialized {len(S_obs[p])}, got {len(obs_stat)}); '
                        f'skipping this protein for this trait'
                    )
                    continue

                count_valid[p] += vmask.astype(np.int32)
                with np.errstate(divide='ignore', invalid='ignore'):
                    S_obs[p] += np.where(vmask, 1.0 / obs_stat, 0.0)
                    S_null[p] += np.where(
                        vmask[:, np.newaxis],
                        1.0 / null_stat,
                        0.0,
                    )

    if n_sims is None:
        logger.warning(f'No valid data found for cluster {cluster}')
        return [], n_sims_expected

    # Compute harmonic means and write to master h5
    entry_ids_written = []
    with h5py.File(master_path, 'a') as fid:
        for p in sorted(S_obs.keys()):
            keep = count_valid[p] > 0
            if not keep.any():
                continue

            cv_k = count_valid[p][keep].astype(float)  # (n_keep,)

            with np.errstate(divide='ignore', invalid='ignore'):
                hm_obs = cv_k / S_obs[p][keep]                       # (n_keep,)
                hm_null = cv_k[:, np.newaxis] / S_null[p][keep, :]   # (n_keep, n_sims)

            # Clamp non-finite values (0/0, k/inf) to neutral
            hm_obs = np.where(np.isfinite(hm_obs), hm_obs, neutral)
            hm_null = np.where(np.isfinite(hm_null), hm_null, neutral)

            key = f'{p}.{cluster}'
            write_dataset(fid, key, hm_obs.reshape(-1, 1))
            write_dataset(fid, f'{key}_null_n_a_tstat', hm_null)
            # mean_beta slot: col 0 = 0.0 (no beta for HMP), col 1 = count_valid
            write_dataset(fid, f'{key}_mean_beta',
                          np.column_stack([np.zeros(int(keep.sum())), cv_k]))
            write_dataset(fid, f'{key}_aa_pos',
                          aa_pos_store[p][keep].reshape(-1, 1))
            entry_ids_written.append(key)
            logger.debug(f'  Wrote {key}: {int(keep.sum())} neighborhoods')

        # Update the entry_ids index dataset
        existing = (
            [e.decode('ascii') for e in fid['entry_ids'][:]]
            if 'entry_ids' in fid else []
        )
        if 'entry_ids' in fid:
            del fid['entry_ids']
        fid.create_dataset(
            'entry_ids',
            data=[e.encode('ascii') for e in existing + entry_ids_written],
        )

    return entry_ids_written, n_sims


def build_hmp_master_h5(trait_cluster_file, results_dir, stat_method='hill',
                        min_variants=10):
    """
    Build the master p_values.h5 from per-trait h5 files.

    Reads a TSV with columns 'trait' and 'cluster', streams each trait's
    q3dnt_results/ukbb_{trait}_pval_all-nbhd_gp250506/p_values.h5, and writes cluster-level harmonic mean
    statistics to {results_dir}/p_values.h5.

    Args:
        trait_cluster_file: path to TSV with columns 'trait' and 'cluster'
        results_dir:        output directory; master h5 written here
        stat_method:        'tstat', 'hill', 'pval', or 'hillp' — determines neutral value
        min_variants:       minimum variants in a neighborhood to count a trait as valid
                            (default 10, matching the q_scan_test threshold)

    Returns:
        path to the master h5 file
    """
    os.makedirs(results_dir, exist_ok=True)
    master_path = os.path.join(results_dir, 'p_values.h5')
    if os.path.exists(master_path):
        os.remove(master_path)

    df_tc = pd.read_csv(trait_cluster_file, sep='\t')
    if not {'trait', 'cluster'}.issubset(df_tc.columns):
        raise ValueError("trait_cluster_file must contain columns 'trait' and 'cluster'")

    neutral = _NEUTRAL.get(stat_method, 0.0)
    clusters = df_tc['cluster'].unique()
    logger.info(
        f'HMP Phase 1: {len(clusters)} clusters, {len(df_tc)} trait rows, '
        f'stat_method={stat_method}, neutral={neutral}, min_variants={min_variants}'
    )

    n_sims = None
    for cluster in clusters:
        traits = df_tc.loc[df_tc['cluster'] == cluster, 'trait'].tolist()
        logger.info(f'Cluster {cluster}: combining {len(traits)} traits')
        ids, n_sims = _build_cluster(
            traits, cluster, master_path, neutral, n_sims, min_variants
        )
        logger.info(f'  -> {len(ids)} entries written')

    logger.info(f'Master h5 complete: path={master_path}, n_sims={n_sims}')
    return master_path
