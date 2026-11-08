import csv
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import shutil
import statistics

repo = Path(__file__).resolve().parents[2]
view = repo / 'runs/beeline/leap_corrected_view_20261105'
out = repo / 'source_data/supplementary_tables/beeline_pseudotime_20261105'
if out.exists():
    raise FileExistsError(out)
out.mkdir()
metrics = ['AP', 'AUPRC', 'AUROC', 'EarlyPrecision', 'EPR', 'ACC', 'Precision', 'Recall', 'F1']


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def rank(rows, suffix):
    cases = defaultdict(list)
    for row in rows:
        cases[row['Case']].append(row)
    for values in cases.values():
        for metric in metrics:
            ordered = sorted(values, key=lambda row: -float(row[metric]))
            positions = defaultdict(list)
            for i, row in enumerate(ordered, 1):
                positions[float(row[metric])].append(i)
            for row in values:
                row[metric + suffix] = statistics.mean(positions[float(row[metric])])


def summaries(rows, suffix):
    answer = []
    for scope in ['all', 'curated', 'synthetic']:
        groups = defaultdict(list)
        for row in rows:
            if scope == 'all' or row['DatasetType'] == scope:
                groups[row['Method']].append(row)
        for method, values in sorted(groups.items()):
            item = dict(Method=method, Scope=scope, Cases=len(values))
            for metric in metrics:
                item[metric] = statistics.mean(float(row[metric]) for row in values)
                item[metric + suffix] = statistics.mean(float(row[metric + suffix]) for row in values)
            answer.append(item)
    return answer


baseline = read(repo / 'results/beeline/ap_sensitivity/case_metrics.tsv')
leap = [row for p in sorted((repo / 'runs/beeline/leap_ordering_20261105').glob('shard_*/case_metrics.tsv')) for row in read(p)]
assert len(leap) == 500 and len({row['Case'] for row in leap}) == 500
replacement = {row['Case']: row for row in leap}
updated = []
leap_changes = []
for old in baseline:
    row = {key: value for key, value in old.items() if key not in ('AP_rank', 'AUPRC_rank')}
    if row['Method'] == 'LEAP':
        corrected = replacement[row['Case']]
        for name in list(row):
            if name in corrected and name not in ['Case', 'DatasetType', 'DatasetGroup', 'Method']:
                row[name] = corrected[name]
        leap_changes.append(dict(Case=row['Case'], DatasetType=row['DatasetType'],
                                 OldAP=old['AP'], CorrectedAP=row['AP'], OldAUPRC=old['AUPRC'], CorrectedAUPRC=row['AUPRC']))
    updated.append(row)
rank(updated, '_rank11')
write(out / 'corrected11_case_metrics.tsv', updated)
write(out / 'corrected11_summary.tsv', summaries(updated, '_rank11'))
write(out / 'LEAP_corrected_case_metrics.tsv', leap)
write(out / 'LEAP_before_after.tsv', leap_changes)
write(view / 'results/beeline/evaluation/case_metrics.tsv', updated)
write(view / 'results/beeline/ap_sensitivity/case_metrics.tsv', updated)
write(view / 'results/beeline/ap_sensitivity/method_summary.tsv', summaries(updated, '_rank11'))
for filename in ['runtime_all_methods.tsv', 'runtime_summary.tsv']:
    path = view / 'results/beeline/runtime' / filename
    write(path, [row for row in read(path) if row['Method'] != 'LEAP'])
(view / 'results/beeline/metrics').rename(view / 'results/beeline/metrics_historical')
for case in ['HSC-2000-3-30', 'VSC-2000-9']:
    source = next((repo / 'runs/beeline/leap_ordering_20261105').glob('shard_*/networks/LEAP__' + case + '.csv'))
    shutil.copyfile(source, view / 'results/beeline/networks' / source.name)
mode = read(repo / 'source_data/supplementary_tables/beeline_static_20261105/mode_case_metrics.tsv')
joint = [dict(row) for row in updated if row['Method'] != 'multiCSN'] + [dict(row) for row in mode]
rank(joint, '_rank12')
assert len(joint) == 6000
write(out / 'joint12_case_metrics.tsv', joint)
write(out / 'joint12_summary.tsv', summaries(joint, '_rank12'))
scenarios = []
scenario_ranks = []
for mode_name in ['static', 'pseudotime']:
    selected = [dict(row) for row in joint if not row['Method'].startswith('multiCSN') or
                row['Method'] == f'multiCSN ({mode_name})']
    rank(selected, '_rank11')
    scenarios.extend(dict(row, Scenario=mode_name) for row in summaries(selected, '_rank11'))
    scenario_ranks.extend({key: row[key] for key in ['Case', 'Method', 'DatasetType', 'DatasetGroup',
                                                    'AP', 'AUPRC', 'AP_rank11', 'AUPRC_rank11']} |
                          {'Scenario': mode_name} for row in selected)
write(out / 'scenario_summary.tsv', scenarios)
write(out / 'scenario_case_ranks.tsv', scenario_ranks)
shuffle_root = repo / 'runs/beeline/pseudotime_shuffle_20261105'
shuffle = read(shuffle_root / 'case_metrics.tsv')
assert len(shuffle) == 1500 and len({(row['Case'], row['ShuffleSeed']) for row in shuffle}) == 1500
for row in shuffle:
    path = shuffle_root / 'networks' / (row['Case'] + '__' + row['ShuffleSeed'] + '.csv')
    with path.open() as handle:
        row['SelectedEdges'] = sum(1 for _ in handle) - 1
write(out / 'state_order_controls.tsv', shuffle)
mode_map = {(row['Case'], row['Method']): row for row in mode}
paired = []
for case in sorted(replacement):
    static = mode_map[(case, 'multiCSN (static)')]
    pseudo = mode_map[(case, 'multiCSN (pseudotime)')]
    controls = [row for row in shuffle if row['Case'] == case]
    item = dict(Case=case, DatasetType=static['DatasetType'], DatasetGroup=static['DatasetGroup'], ShuffleFits=len(controls))
    for metric in ['AP', 'AUPRC']:
        item['Static_' + metric] = static[metric]
        item['Pseudotime_' + metric] = pseudo[metric]
        item['Shuffled_' + metric] = statistics.mean(float(row[metric]) for row in controls)
        item['PseudotimeMinusStatic_' + metric] = float(pseudo[metric]) - float(static[metric])
        item['PseudotimeMinusShuffled_' + metric] = float(pseudo[metric]) - item['Shuffled_' + metric]
    item['EmptyShuffleFits'] = sum(row['SelectedEdges'] == 0 for row in controls)
    paired.append(item)
write(out / 'paired_mode_controls.tsv', paired)
control_summary = []
for scope in ['all', 'curated', 'synthetic']:
    values = [row for row in paired if scope == 'all' or row['DatasetType'] == scope]
    controls = [row for row in shuffle if scope == 'all' or row['DatasetType'] == scope]
    item = dict(Scope=scope, Cases=len(values), ShuffleFits=len(controls),
                EmptyShuffleFits=sum(row['SelectedEdges'] == 0 for row in controls))
    for name in ['Static_AP', 'Pseudotime_AP', 'Shuffled_AP', 'Static_AUPRC', 'Pseudotime_AUPRC', 'Shuffled_AUPRC']:
        item[name] = statistics.mean(float(row[name]) for row in values)
    for metric in ['AP', 'AUPRC']:
        for comparator in ['Static', 'Shuffled']:
            name = 'PseudotimeMinus' + comparator + '_' + metric
            item[name] = statistics.mean(float(row[name]) for row in values)
            item[name + '_Wins'] = sum(float(row[name]) > 1e-12 for row in values)
            item[name + '_Ties'] = sum(abs(float(row[name])) <= 1e-12 for row in values)
    control_summary.append(item)
write(out / 'control_summary.tsv', control_summary)
write(out / 'model_control_summary.tsv', [dict(Model=model, Cases=len(values), **{name:statistics.mean(float(row[name]) for row in values) for name in ['PseudotimeMinusStatic_AP','PseudotimeMinusShuffled_AP','PseudotimeMinusStatic_AUPRC','PseudotimeMinusShuffled_AUPRC']}) for model in sorted({row['DatasetGroup'] for row in paired}) for values in [[row for row in paired if row['DatasetGroup'] == model]]])
print('Corrected 11-method summary:')
for row in sorted([row for row in summaries(updated, '_rank11') if row['Scope']=='all'], key=lambda row:-row['AP']):
    print(row['Method'], round(row['AP'],6), round(row['AUPRC'],6), row['AP_rank11'], row['AUPRC_rank11'])
print('Control summary:', control_summary)
