"""Portable study-directed and study-withheld dsRNASeeker supervised transfer.

Uses the SAME fitting functions as modules.supervised; never substitutes new ADPS.
Train/test datasets must have precomputed candidate evidence and genuine RIP labels.
Output is exploratory for previously inspected study targets.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import platform
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
from modules import supervised as sup

FEATURES = (
    'orientation_adps', 'annotation_adps', 'case_expression_adps',
    'energy_adps', 'interface_adps', 'case_editing_adps', 'RI_adps',
)
PANELS = {
    'raw7': list(FEATURES),
    'core3': ['orientation_adps', 'annotation_adps', 'case_expression_adps'],
    'no_orientation6': [x for x in FEATURES if x != 'orientation_adps'],
    'structure5': ['orientation_adps', 'annotation_adps', 'case_expression_adps',
                   'energy_adps', 'interface_adps'],
}
WEIGHTS = ('class_balanced', 'family_class_balanced')


def _hash(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def _read(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f'Input file not found: {path}')
    return pd.read_csv(path, sep='\t' if path.suffix.lower() in ('.tsv', '.txt') else ',', low_memory=False)


def _resolve(base: Path, value: object) -> Path:
    p = Path(str(value)).expanduser()
    return (p if p.is_absolute() else base / p).resolve()


def _numeric_score(df: pd.DataFrame, *, study: str, path: Path,
                   summary: Path | None = None) -> tuple[pd.Series, str, str]:
    """Use a real recorded ADPS score; join supplied original summary only if absent.

    Never derive integrated ADPS from component means or weighted guesses.
    """
    candidates = [c for c in ('adaptive_priority_score', 'original_ADPS') if c in df]
    if candidates:
        score = pd.to_numeric(df[candidates[0]], errors='coerce')
        if len(candidates) == 2:
            secondary = pd.to_numeric(df[candidates[1]], errors='coerce')
            if not np.isfinite(secondary.to_numpy(float)).all() or not np.allclose(score, secondary, atol=1e-12, rtol=0):
                raise ValueError(f'{study}: conflicting ADPS columns in {path}; no silent selection')
        origin = f'audit_table:{candidates[0]}'
    else:
        if summary is None:
            raise ValueError(
                f'{study}: missing canonical ADPS in {path}. '
                'Provide original_summary in manifest with a real dsRNASeeker '
                'adaptive_priority_score table; available columns: '+str(list(df.columns))+'; '
                'a fabricated substitute is prohibited.')
        s = _read(summary)
        if 'pair_id' not in s or 'adaptive_priority_score' not in s:
            raise ValueError(f'{study}: original summary lacks pair_id/adaptive_priority_score: {summary}')
        if s.pair_id.isna().any() or s.pair_id.astype(str).duplicated().any():
            raise ValueError(f'{study}: ambiguous duplicate pair_id in original summary: {summary}')
        joined = df[['pair_id']].merge(s[['pair_id', 'adaptive_priority_score']], on='pair_id',
                    how='left', validate='one_to_one', sort=False, indicator=True)
        if not joined['_merge'].eq('both').all():
            raise ValueError(f'{study}: source summary missing labeled pair_id in {summary}')
        score = pd.to_numeric(joined['adaptive_priority_score'], errors='coerce')
        origin = f'original_summary:{summary};sha256={_hash(summary)}'
    if not np.isfinite(score.to_numpy(float)).all():
        raise ValueError(f'{study}: non-finite original ADPS for {path}; origin={origin}')
    # Do not auto-replace an audited frozen ADPS merely because its source differs.
    status = 'NO_SOURCE_COMPARISON'
    if summary is not None:
        if not summary.is_file():
            raise FileNotFoundError(summary)
        s = _read(summary)
        if 'pair_id' in s and 'adaptive_priority_score' in s:
            if s.pair_id.astype(str).duplicated().any():
                raise ValueError(f'{study}: duplicate source summary pair IDs: {summary}')
            joined = df[['pair_id']].merge(s[['pair_id', 'adaptive_priority_score']], on='pair_id',
                 how='left',validate='one_to_one',sort=False,indicator=True)
            source_score = pd.to_numeric(joined['adaptive_priority_score'],errors='coerce')
            mismatches = int((~np.isclose(score, source_score, rtol=0, atol=1e-8, equal_nan=False)).sum())
            status = f'SOURCE_MATCH' if mismatches == 0 else f'SOURCE_DRIFT_{mismatches}_OF_{len(df)}'
    return score.reset_index(drop=True), origin, status


def _load_studies(manifest: Path):
    registry = _read(manifest)
    required = {'study_id', 'audit_table', 'label_definition', 'analysis_variant'}
    if not required.issubset(registry.columns):
        raise ValueError(f'Manifest {manifest} missing: {sorted(required - set(registry.columns))}')
    if len(registry) < 3:
        raise ValueError('Source-only selection requires at least three distinct independent study families')
    if registry['study_id'].isna().any() or registry['study_id'].astype(str).duplicated().any():
        raise ValueError('Duplicate or empty study_id in manifest')
    if 'study_family' in registry.columns:
        families = registry['study_family'].astype(str)
        if families.eq('nan').any() or families.duplicated().any():
            raise ValueError('Manifest contains duplicated study_family; exclude same-family overlap before transfer')
    frames: dict[str, pd.DataFrame] = {}
    records: list[dict[str, Any]] = []
    for row in registry.itertuples(index=False):
        study = str(row.study_id)
        path = _resolve(manifest.parent, row.audit_table)
        df = _read(path)
        fields = {'pair_id', 'group_id', 'label', *FEATURES}
        missing = fields - set(df.columns)
        if missing:
            raise ValueError(f'{study}: missing {sorted(missing)} in {path}')
        if df['pair_id'].isna().any() or df['pair_id'].astype(str).str.strip().eq('').any() or df['pair_id'].astype(str).duplicated().any():
            raise ValueError(f'{study}: missing or duplicate pair IDs in {path}')
        raw = df['label']
        y = pd.to_numeric(raw, errors='coerce')
        invalid = raw.notna() & raw.astype(str).str.strip().ne('') & y.isna()
        if invalid.any() or not y.dropna().isin((0, 1)).all():
            raise ValueError(f'{study}: invalid or fractional binary label in {path}')
        d = df.loc[y.notna()].copy().reset_index(drop=True)
        d['label'] = y[y.notna()].astype(int).to_numpy()
        npos = int(d['label'].sum())
        nneg = int(len(d) - npos)
        if min(npos, nneg) < 2:
            raise ValueError(f'{study}: cannot evaluate with {npos} positives and {nneg} negatives')
        groups = d['group_id']
        if groups.isna().any() or groups.astype(str).str.strip().eq('').any() or groups.astype(str).str.lower().isin(('nan', 'none', 'null')).any():
            raise ValueError(f'{study}: invalid connected-locus group_id')
        d['group_id'] = groups.astype(str)
        for c in FEATURES:
            x = pd.to_numeric(d[c], errors='coerce')
            if x.notna().sum() == 0 or not np.isfinite(x.dropna().to_numpy(float)).all():
                raise ValueError(f'{study}: feature {c} invalid/all missing in {path}')
            d[c] = x
        source = getattr(row, 'original_summary', None)
        source_path = (_resolve(manifest.parent, source) if isinstance(source, str)
                       and source.strip() and source.lower().strip() != 'nan' else None)
        score, score_origin, source_status = _numeric_score(
            d, study=study, path=path, summary=source_path)
        d['original_ADPS'] = score.to_numpy(float)
        d['study_family'] = study  # Must be independent unique study family, enforced above.
        for name, count in (('expected_positive', npos), ('expected_negative', nneg),
                            ('expected_groups', d['group_id'].nunique())):
            want = getattr(row, name, np.nan)
            if pd.notna(want) and int(want) != int(count):
                raise ValueError(f'{study}: {name} expected {want}, observed {count}')
        frames[study] = d
        records.append({'study_id': study, 'label_definition': row.label_definition,
                        'analysis_variant': row.analysis_variant, 'path': str(path),
                        'sha256': _hash(path), 'candidate_universe': len(df),
                        'n_labeled': len(d), 'positive': npos, 'negative': nneg,
                        'groups': d['group_id'].nunique(),
                        'adps_col_in_input': 'adaptive_priority_score' if 'adaptive_priority_score' in df else ('original_ADPS' if 'original_ADPS' in df else 'MISSING_JOINED_FROM_SUMMARY'),
                        'adps_source_origin':score_origin, 'adps_source_concordance':source_status})
    return frames, records


def _metrics(y, scores):
    a = np.asarray(scores, dtype=float)
    target = np.asarray(y, dtype=int)
    if not np.isfinite(a).all() or set(np.unique(target)) != {0, 1}:
        raise ValueError('Invalid binary performance inputs')
    return {'roc_auc': float(roc_auc_score(target, a)),
            'average_precision': float(average_precision_score(target, a))}


def _fit_native(train: pd.DataFrame, test: pd.DataFrame, panel: str,
                weighting: str, seed: int) -> np.ndarray:
    cols = PANELS[panel]
    y = train['label'].astype(int).reset_index(drop=True)
    x = train[cols].apply(pd.to_numeric, errors='coerce')
    families = train['study_family'].astype(str).reset_index(drop=True)
    if weighting == 'family_class_balanced':
        weights = sup._family_class_balanced_weights(y, families)
    elif weighting == 'class_balanced':
        weights = sup._class_balanced_weights(y)
    else:
        raise ValueError(f'Unknown weighting policy: {weighting}')
    model = sup._model(int(seed), model_type='legacy_l2', C=1.0, l1_ratio=None)
    model, converged = sup._fit_checked(model, x, y, weights)
    if not converged:
        raise RuntimeError(f'Native estimator failed to converge: {panel}, {weighting}')
    score = model.predict_proba(test[cols].apply(pd.to_numeric, errors='coerce'))[:, 1]
    if not np.isfinite(score).all():
        raise RuntimeError('Native model produced non-finite scores')
    return score


def _inner_selection(source: dict[str, pd.DataFrame], seed: int):
    if len(source) < 2:
        raise ValueError('Not enough source families for study-wise inner selection')
    rows = []
    for panel, weight in itertools.product(PANELS, WEIGHTS):
        for k, (held, test) in enumerate(source.items()):
            train = pd.concat([d for sid, d in source.items() if sid != held], ignore_index=True)
            score = _fit_native(train, test, panel, weight, seed + k)
            m = _metrics(test['label'], score)
            rows.append({'panel': panel, 'weight': weight, 'inner_held_out': held,
                         'inner_auc': m['roc_auc'], 'inner_ap': m['average_precision'],
                         'inner_ap_lift': m['average_precision'] - float(test.label.mean()),
                         'inner_train_n': len(train), 'inner_test_n': len(test)})
    detail = pd.DataFrame(rows)
    ranked = detail.groupby(['panel', 'weight'], as_index=False).agg(
        mean_inner_auc=('inner_auc', 'mean'), mean_inner_ap_lift=('inner_ap_lift', 'mean'),
        min_inner_auc=('inner_auc', 'min'), n_sources=('inner_held_out', 'nunique'))
    ranked = ranked.sort_values(['mean_inner_auc', 'mean_inner_ap_lift', 'panel', 'weight'],
                                ascending=[False, False, True, True], kind='stable').reset_index(drop=True)
    return detail, ranked


def _group_bootstrap(test, score, adps, *, seed: int, nboot: int):
    if nboot == 0:
        return {'bootstrap_valid': 0, 'roc_auc_ci_low': np.nan, 'roc_auc_ci_high': np.nan,
                'delta_auc_ci_low': np.nan, 'delta_auc_ci_high': np.nan}
    rng = np.random.default_rng(seed)
    grp = test['group_id'].astype(str).to_numpy()
    indexes = [np.flatnonzero(grp == g) for g in np.unique(grp)]
    y = test['label'].to_numpy(dtype=int)
    auc, delta = [], []
    for _ in range(nboot):
        ix = np.concatenate([indexes[j] for j in rng.integers(len(indexes), size=len(indexes))])
        if np.unique(y[ix]).size != 2:
            continue
        a = float(roc_auc_score(y[ix], score[ix]))
        b = float(roc_auc_score(y[ix], adps[ix]))
        auc.append(a)
        delta.append(a-b)
    def q(x):
        return (float(np.quantile(x, .025)), float(np.quantile(x, .975))) if x else (np.nan, np.nan)
    a, b = q(auc), q(delta)
    return {'bootstrap_valid': len(auc), 'roc_auc_ci_low': a[0], 'roc_auc_ci_high': a[1],
            'delta_auc_ci_low': b[0], 'delta_auc_ci_high': b[1]}


def run_transfer(args):
    import sklearn
    manifest = Path(args.manifest).expanduser().resolve()
    out = Path(args.output_dir).expanduser().resolve()
    if args.bootstrap < 0:
        raise ValueError('Bootstrap count cannot be negative')
    if out.exists() and (out / 'COMPLETE.txt').exists() and not args.overwrite:
        raise FileExistsError(f'Completed output already exists; choose new output-dir or --overwrite: {out}')
    studies, recs = _load_studies(manifest)  # Fail before writing results.
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(recs).to_csv(out / 'frozen_inputs.csv', index=False)
    adps_rows = []
    for sid, d in studies.items():
        adps_rows.append({'study_id': sid, 'positives': int(d.label.sum()), 'negatives': int(len(d)-d.label.sum()),
                          **_metrics(d['label'], d['original_ADPS'])})
    pd.DataFrame(adps_rows).to_csv(out / 'adps_metrics.csv', index=False)
    prediction_rows = []
    if args.mode in ('directed', 'both'):
        rows = []
        for i, (source_id, target_id) in enumerate(itertools.permutations(studies.keys(), 2)):
            train, test = studies[source_id], studies[target_id]
            for panel in ('raw7', 'core3'):
                score = _fit_native(train, test, panel, 'class_balanced', args.seed + i)
                m = _metrics(test['label'], score)
                rows.append({'train': source_id, 'test': target_id, 'panel': panel,
                             'weight': 'class_balanced', 'n_train': len(train), 'n_test': len(test),
                             'positive': int(test.label.sum()), 'negative': int(len(test)-test.label.sum()),
                             'ap_baseline': float(test.label.mean()), **m})
                prediction_rows.append(pd.DataFrame({'evaluation': 'directed', 'train': source_id,
                    'test': target_id, 'panel': panel, 'pair_id': test.pair_id,
                    'group_id': test.group_id, 'label': test.label, 'prediction': score,
                    'original_ADPS': test.original_ADPS}))
        pd.DataFrame(rows).to_csv(out/'directed_transfer_metrics.csv', index=False)
    if args.mode in ('source-only', 'both'):
        out_rows, candidates, inner_folds = [], [], []
        for i, (target_id, test) in enumerate(studies.items()):
            source = {k: v for k, v in studies.items() if k != target_id}
            detail, rank = _inner_selection(source, args.seed)
            detail.insert(0, 'outer_target', target_id)
            rank.insert(0, 'outer_target', target_id)
            inner_folds.append(detail)
            candidates.append(rank)
            choice = rank.iloc[0]
            tr = pd.concat(list(source.values()), ignore_index=True)
            score = _fit_native(tr, test, str(choice.panel), str(choice.weight), args.seed + i)
            m = _metrics(test['label'], score)
            adps_m = _metrics(test['label'], test['original_ADPS'])
            bootstrap = _group_bootstrap(test, score, test.original_ADPS.to_numpy(float),
                                         seed=args.seed + i, nboot=args.bootstrap)
            out_rows.append({'held_out': target_id, 'train_families': ';'.join(source.keys()),
                             'selected_panel': choice.panel, 'selected_weight': choice.weight,
                             'source_only_inner_mean_auc': float(choice.mean_inner_auc),
                             'n_train': len(tr), 'n_test': len(test),
                             'positive': int(test.label.sum()), 'negative': int(len(test)-test.label.sum()),
                             'ap_baseline': float(test.label.mean()),
                             'roc_auc': m['roc_auc'], 'average_precision': m['average_precision'],
                             'adps_auc': adps_m['roc_auc'], 'adps_ap': adps_m['average_precision'],
                             'delta_auc_vs_adps': m['roc_auc'] - adps_m['roc_auc'],
                             'delta_ap_vs_adps': m['average_precision'] - adps_m['average_precision'],
                             **bootstrap, 'status': 'EXPLORATORY_KNOWN_TARGETS'})
            prediction_rows.append(pd.DataFrame({'evaluation': 'source_only',
                'train': ';'.join(source.keys()), 'test': target_id, 'panel': choice.panel,
                'pair_id': test.pair_id, 'group_id': test.group_id, 'label': test.label,
                'prediction': score, 'original_ADPS': test.original_ADPS}))
            print(f'{target_id}: selected {choice.panel}/{choice.weight} from sources; '
                  f'target AUC={m["roc_auc"]:.4f}, AP={m["average_precision"]:.4f}', flush=True)
        pd.DataFrame(out_rows).to_csv(out/'outer_source_only_transfer_metrics.csv', index=False)
        pd.concat(candidates, ignore_index=True).to_csv(out/'inner_source_only_candidates.csv', index=False)
        pd.concat(inner_folds, ignore_index=True).to_csv(out/'inner_source_only_folds.csv', index=False)
    pd.concat(prediction_rows, ignore_index=True).to_csv(out/'all_predictions.csv', index=False)
    meta = {'entrypoint': 'python main.py supervised-transfer',
            'model_function': 'modules.supervised._model + _fit_checked',
            'model_source': str(Path(sup.__file__).resolve()),
            'model_sha256': _hash(Path(sup.__file__).resolve()),
            'transfer_module_source': str(Path(__file__).resolve()),
            'transfer_module_sha256': _hash(Path(__file__).resolve()),
            'manifest': str(manifest), 'manifest_sha256': _hash(manifest),
            'mode': args.mode, 'n_studies': len(studies), 'seed': args.seed,
            'bootstrap': args.bootstrap, 'python_library': platform.python_version(),
            'sklearn': sklearn.__version__, 'pandas': pd.__version__,
            'study_selection': 'source families only; no outer target labels until evaluation',
            'biological_data_edited': False, 'original_adps_fabricated': False,
            'interpretation': 'Previously inspected cohorts: retrospective development, not fresh validation',
            'note': 'Direct-study transfer, source-only selected transfer and pooled LOSO are different estimands'}
    (out/'run_metadata.json').write_text(json.dumps(meta, indent=2)+'\n')
    (out/'COMPLETE.txt').write_text('Native dsRNASeeker supervised-transfer completed. Prior cohorts were examined previously; exploratory results only.\n')
    print('[COMPLETE]', out, flush=True)
    return out


def add_supervised_transfer_args(parser: argparse.ArgumentParser):
    parser.add_argument('--manifest', required=True, help='Portable TSV/CSV: study_id, audit_table, label_definition, analysis_variant; optional expected counts, study_family')
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--mode', choices=('directed', 'source-only', 'both'), default='both')
    parser.add_argument('--bootstrap', type=int, default=500)
    parser.add_argument('--seed', type=int, default=20261011)
    parser.add_argument('--overwrite', action='store_true', help='Permit recomputing in an existing completed output directory')


def cli():
    parser = argparse.ArgumentParser(description=__doc__)
    add_supervised_transfer_args(parser)
    args = parser.parse_args()
    path = run_transfer(args)
    print('[COMPLETE]', path)


if __name__ == '__main__':
    cli()
