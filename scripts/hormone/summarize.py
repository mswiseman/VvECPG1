#!/usr/bin/env python3
"""Reconstruct the prespecified ten-test focal family from existing Wald P-values.

This does not refit DESeq2. MLE and apeglm estimates are retained in separate columns;
neither estimate replaces the Wald P-value. Python standard library only.
"""
from pathlib import Path
import csv
import math

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'results/hormone'
OUT.mkdir(parents=True, exist_ok=True)

def read(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def bh(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    result = [0.0] * len(values)
    running = 1.0
    for rank in range(len(order), 0, -1):
        i = order[rank - 1]
        running = min(running, values[i] * len(order) / rank)
        result[i] = running
    return result

full = read(ROOT / 'inputs/hormone/focal_gene_de_40_contrasts.tsv')
plot = read(ROOT / 'inputs/hormone/focal_24h_vs_control.tsv')
assert len(full) == 40 and all(row['gene_id'] == 'Vitvi05_01chr01g05750' for row in full)
primary = [row for row in full if row['condition'].endswith('_24h') and row['batch'] == '1']
assert len(primary) == 10 and len({row['condition'] for row in primary}) == 10
q = bh([float(row['pvalue']) for row in primary])
plot_by_condition = {row['condition']: row for row in plot}
for row, adjusted in zip(primary, q):
    original = plot_by_condition[row['condition']]
    assert row['control_samples'] == 'PPPJ8P_94,PPPJ8P_95,PPPJ8P_96'
    assert row['n_treated'] == row['n_control'] == '3'
    assert math.isclose(adjusted, float(original['padj_10_tests_BH']), abs_tol=1e-12)
    for left, right in [('pvalue', 'pvalue'), ('log2FC', 'log2FC'),
                        ('log2FC_shrunk_apeglm', 'log2FC_shrunk'), ('padj_BH', 'padj_genomewide')]:
        assert math.isclose(float(row[left]), float(original[right]), rel_tol=1e-12, abs_tol=1e-14)
    row['padj_focal_10_tests_BH'] = adjusted
    row['significant_focal_24h'] = adjusted < 0.05
    row['mean_batch1_normalized_count'] = original['mean']
    row['SE_batch1_normalized_count'] = original['se']
significant = {row['condition'] for row in primary if row['significant_focal_24h']}
assert significant == {'ABA_7p5mgL_24h', 'GA3_1mgL_24h', 'MeJA_100uM_24h'}
with (OUT / 'focal_24h_primary_tests.tsv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(primary[0]), delimiter='\t')
    writer.writeheader()
    writer.writerows(primary)
with (OUT / 'verification.txt').open('w') as handle:
    handle.write('PASS: ten 24-h focal BH values, original Wald P-values, MLE/apeglm estimates,\n'
                 'control identities, replicate counts and the three significant responses verified.\n')
print('Verified ten time-matched focal tests; three pass BH < 0.05.')
