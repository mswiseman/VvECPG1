import os
#!/usr/bin/env python3
"""Classify ENA Vitis RNA-seq runs into tissue/stage classes for the focal-gene atlas pilot.

Input : results/sra_atlas/metadata/ena_vitis_rnaseq_runs.tsv (ENA portal read_run, tax_tree(3603)).
Output: ena_runs_classified.tsv (+ tissue column, treatment flag, QC-pass flag) and a study x tissue summary.
Keyword rules on tissue_type, then sample/experiment/study titles and descriptions. First match wins,
specific classes (somatic embryo, callus) before generic ones (leaf). Rules are deliberately simple —
the output is a candidate list for manual curation, not a final annotation.
"""
import re
import pandas as pd

D = os.path.join(os.environ.get("VVECPG1_UPSTREAM_ROOT", "."), 'results/sra_atlas/metadata')
r = pd.read_csv(f"{D}/ena_vitis_rnaseq_runs.tsv", sep="\t", dtype=str).fillna("")

RULES = [  # (class, regex)
    ("non-embryogenic callus", r"not embryogenic|non-?embryogenic"),
    ("somatic embryo", r"somatic embryo|embryogen(?!ic callus)|proembryo|\bSE\b"),
    ("embryogenic callus", r"embryogenic callus|embryogenic cultur|\bPEM\b"),
    ("callus / cell culture", r"callus|calli|cell suspension|suspension cell|cell culture|cell line"),
    ("seed / zygotic embryo", r"\bseeds?\b|zygotic|endosperm|seed coat"),
    ("pollen", r"pollen"),
    ("flower / inflorescence", r"flower|inflorescen|anther|stamen|pistil|ovar|carpel|petal|bloom|anthesis"),
    ("berry skin", r"\bskin|exocarp|epicarp"),
    ("berry pulp / flesh", r"pulp|flesh|mesocarp"),
    ("berry / fruit", r"berr(y|ies)|fruit|grape cluster|bunch|veraison|pericarp"),
    ("tendril", r"tendril"),
    ("root", r"\broots?\b|rootstock root|root tip"),
    ("bud", r"\bbuds?\b|dormant bud|latent bud|meristem"),
    ("stem / cane / xylem / phloem", r"\bstems?\b|cane|xylem|phloem|wood|trunk|internode|cambium|bark|shoot apex"),
    ("leaf", r"\bleaf|leaves|lamina|foliar|petiole"),
    ("shoot / seedling / whole plant", r"shoots?\b|seedling|plantlet|whole plant|in vitro plant"),
]
TREAT = r"infect|inocul|pathogen|mildew|botrytis|erysiphe|plasmopara|virus|stress|drought|salt|cold|heat|" \
        r"freez|uv|treat|hormone|aba\b|abscisic|ethylene|methyl jasmon|meja|gibberell|ga3|auxin|naa\b|" \
        r"cytokinin|brassinost|elicit|chitosan|wound|mutant|transgen|overexpress|knock|silenc|deficien|" \
        r"nitrogen|phosphat|potassium|water deficit|flood|hypox|shade|girdl|grafting"

def classify(row):
    t = row.tissue_type.lower()
    txt = " ".join([row.sample_title, row.sample_description, row.experiment_title]).lower()
    st = row.study_title.lower()
    for label, src in (("tissue_type", t), ("sample text", txt), ("study title", st)):
        if not src.strip():
            continue
        src = re.sub(r"(no|without|removed|minus) seeds?|seeds? (removed|excluded)|seedless|seed trace", " ", src)
        for cls, rx in RULES:
            if re.search(rx, src, re.I):
                return cls, label
    return "unclassified", ""

r[["tissue_class", "tissue_source"]] = r.apply(lambda x: pd.Series(classify(x)), axis=1)
alltxt = (r.tissue_type + " " + r.sample_title + " " + r.sample_description + " " + r.experiment_title + " " +
          r.study_title).str.lower()
r["treatment_mentioned"] = alltxt.str.contains(TREAT, regex=True)
r["read_count"] = pd.to_numeric(r.read_count, errors="coerce")
r["qc_pass"] = ((r.instrument_platform == "ILLUMINA") & (r.read_count >= 10e6) &
                (~r.library_selection.str.contains("size fractionation|miRNA", case=False)) &
                (r.scientific_name == "Vitis vinifera"))
r.to_csv(f"{D}/ena_runs_classified.tsv", sep="\t", index=False)

q = r[r.qc_pass]
s = (q.groupby("tissue_class")
       .agg(runs=("run_accession", "size"), studies=("study_accession", "nunique"),
            runs_no_treatment_words=("treatment_mentioned", lambda v: int((~v).sum())),
            cultivars=("cultivar", lambda v: v[v != ""].nunique()))
       .sort_values("runs", ascending=False))
s.to_csv(f"{D}/tissue_class_summary.tsv", sep="\t")
print(f"all runs {len(r)}; V. vinifera {int((r.scientific_name=='Vitis vinifera').sum())}; QC-pass {len(q)}; "
      f"studies {q.study_accession.nunique()}")
print(s.to_string())
print("\nscientific names:", r.scientific_name.value_counts().head(8).to_dict())
print("platforms (QC-pass):", q.instrument_model.value_counts().head(6).to_dict())
print("tissue_type filled (QC-pass):", int((q.tissue_type != "").sum()))
