import os
#!/usr/bin/env python3
"""Run selection for the focal-gene tissue/stage atlas pilot, from the curated study table.

curated_studies.tsv (user-approved 2026-11-06, option 2: draft fixes + replacement leaf/root/stem studies):
  tissue_class_override  relabel every run of the study to this class (e.g. skin+pulp -> whole berry)
  keep_classes           comma list; keep only runs classified into these classes (blank = all)
  control_only           yes -> keep only runs whose own sample text says control/mock/untreated
  cultivar_override      cultivar for every run of the study (used where ENA leaves it blank; source in note)
  exclude_reason         non-empty -> study left out entirely (reason kept in the table)
  max_runs               cap per study x class (spread across samples by taking every k-th run)
Runs must be QC-pass (V. vinifera, Illumina, >= 10 M reads, not small RNA) and, unless the study is
control_only, not mention a treatment in their own sample text.
run_overrides.tsv: run-level metadata corrections (regex on title fields -> cultivar/species/tissue_class, with source).
class_basis (2026-11-08): 'submitter' / 'curated relabel' / 'inferred' / 'interpretation flagged'; submitter_class keeps the
keyword class from submitter metadata. Inferred groups get their OWN class name and are never merged into a submitter class.
Output: selection.tsv (one row per run), selection_studies.tsv, run_list.txt.
"""
import pandas as pd

D = os.path.join(os.environ.get("VVECPG1_UPSTREAM_ROOT", "."), 'results/sra_atlas/metadata')
r = pd.read_csv(f"{D}/ena_runs_classified.tsv", sep="\t", dtype=str).fillna("")
cur = pd.read_csv(f"{D}/curated_studies.tsv", sep="\t", dtype=str).fillna("")
r = r[(r.qc_pass == "True") & r.study_accession.isin(cur.study_accession)].copy()
r["read_count"] = r.read_count.astype(float)

TREAT = (r"infect|inocul|pathogen|mildew|botrytis|virus|stress|drought|salt|cold|heat|freez|\buv\b|treat|"
         r"aba\b|abscisic|ethylene|jasmon|meja|gibberell|ga3|auxin|naa\b|cytokinin|elicit|chitosan|wound|"
         r"mutant|transgen|overexpress|silenc|deficien|water deficit|flood|hypox|shade|girdl|herbiv|beetle|fung")
CTRL = r"control|mock|untreated|\bck\b|non-?treated|\b0 ?h\b|wild ?type|\bwt\b"
own = (r.sample_title + " " + r.sample_description + " " + r.experiment_title).str.lower()
r["own_ctrl"] = own.str.contains(CTRL)
r["own_treat"] = own.str.contains(TREAT) & ~r.own_ctrl

picked = []
for _, c in cur.iterrows():
    if c.get("exclude_reason", ""):
        print(f"excluded {c.study_accession}: {c.exclude_reason[:80]}")
        continue
    g = r[r.study_accession == c.study_accession].copy()
    g["submitter_class"] = g.tissue_class          # keyword class from submitter metadata, before any override
    g["class_basis"] = c.get("class_basis", "") or ("curated relabel" if c.tissue_class_override else "submitter")
    if c.tissue_class_override:
        g["tissue_class"] = c.tissue_class_override
    if c.get("cultivar_override", ""):
        g["cultivar"] = c.cultivar_override
    if c.keep_classes:
        g = g[g.tissue_class.isin([k.strip() for k in c.keep_classes.split(",")])]
    g = g[g.own_ctrl] if c.control_only == "yes" else g[~g.own_treat]
    g = g[g.tissue_class != "unclassified"]
    for cls, gc in g.groupby("tissue_class"):
        gc = gc.sort_values("run_accession")
        cap = int(c.max_runs)
        k = max(1, len(gc) // cap)
        picked.append(gc.iloc[::k].head(cap).assign(curation_note=c.note))
sel = pd.concat(picked).drop_duplicates("run_accession")
# run-level corrections (run_overrides.tsv): regex on the named title fields -> cultivar / species
sel["species"] = sel.scientific_name
sel["metadata_override"] = ""
ovp = f"{D}/run_overrides.tsv"
import os, re
if os.path.exists(ovp):
    for o in pd.read_csv(ovp, sep="\t", dtype=str).fillna("").itertuples():
        hit = pd.Series(False, index=sel.index)
        for fld in o.match_fields.split(","):
            hit |= sel[fld.strip()].str.contains(o.pattern, regex=True)
        if o.study_accession not in ("", "*"):
            hit &= sel.study_accession == o.study_accession
        if o.cultivar:
            sel.loc[hit, "cultivar"] = o.cultivar
        if o.species:
            sel.loc[hit, "species"] = o.species
        if getattr(o, "tissue_class", ""):
            sel.loc[hit, "tissue_class"] = o.tissue_class
            sel.loc[hit, "class_basis"] = getattr(o, "class_basis", "") or "curated relabel"
        sel.loc[hit, "metadata_override"] = o.source
        print(f"override '{o.pattern}': {int(hit.sum())} runs ({', '.join(sorted(sel.loc[hit, 'study_accession'].unique()))})")
cols = ["tissue_class", "submitter_class", "class_basis", "study_accession", "run_accession", "sample_accession", "species", "cultivar", "tissue_type",
        "dev_stage", "sample_title", "library_layout", "read_count", "instrument_model", "study_title",
        "experiment_title", "curation_note", "metadata_override"]
sel = sel[cols].sort_values(["tissue_class", "study_accession", "run_accession"])
sel.to_csv(f"{D}/selection.tsv", sep="\t", index=False)
sel.run_accession.to_csv(f"{D}/run_list.txt", index=False, header=False)
ss = sel.groupby(["tissue_class", "study_accession"]).agg(
    runs=("run_accession", "size"), cultivar=("cultivar", lambda v: "; ".join(sorted(set(v) - {""}))[:40]),
    note=("curation_note", "first")).reset_index()
ss.to_csv(f"{D}/selection_studies.tsv", sep="\t", index=False)
missing = set(cur.study_accession) - set(sel.study_accession)
print(f"selection: {len(sel)} runs, {sel.study_accession.nunique()} studies; curated studies with 0 runs: {sorted(missing)}")
print(sel.groupby("tissue_class").agg(runs=("run_accession", "size"), studies=("study_accession", "nunique")).to_string())
