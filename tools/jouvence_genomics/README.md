# Jouvence genomics research workflow

The public function `analyze_without_bam()` in `tools/analyze_dante_vcf.py`
performs a reusable VCF-only screen. It accepts SNP, indel, CNV and SV VCFs,
but never interprets an absent VCF record as a reliable reference genotype.

## What the no-BAM path adds

1. Exact reviewed ClinVar pathogenic/likely pathogenic matches.
2. GENCODE v19 principal-transcript consequences for variants in GenCC disease genes.
3. Conservative frameshift, stop-gain and canonical-splice pLoF prioritisation.
4. Separate missense/in-frame VUS output; these are not counted as pathogenic carriers.
5. Global and non-Finnish European gnomAD exome frequencies.
6. Existing phase-block handling for possible compound heterozygosity.

Loss of function is considered only when the gene has an autosomal-recessive
GenCC association and at least two reviewed ClinVar P/LP LoF precedents. This is
a transparent research heuristic, not a clinical ACMG classification. Exact
ClinVar benign/likely-benign alleles are excluded, and an exact ClinVar VUS can
never be promoted to the strict novel-pLoF tier.

## Public references

Prepare the compact reusable reference files once:

```python
from pathlib import Path

from tools.jouvence_genomics.references import prepare_no_bam_references

prepare_no_bam_references(Path("data/reference"))
```

Run Python through `uv`, for example:

```shell
uv run python -c 'from pathlib import Path; from tools.jouvence_genomics.references import prepare_no_bam_references; prepare_no_bam_references(Path("data/reference"))'
```

The analysis additionally needs the existing GRCh37 ClinVar VCF, ClinVar
variant summary, GenCC submissions, GENCODE v19 GTF and carrier-frequency file.

## Interpretation tiers

- `candidate_plof_requires_confirmation`: rare, technically credible pLoF in a
  supported recessive LoF disease gene. It must still be reviewed and confirmed.
- `protein_altering_vus`: rare missense, start/stop-loss or in-frame change.
  It is explicitly not counted as a causal or carrier finding.
- exact ClinVar P/LP matches retain the existing, separate reporting path.

## When BAM/CRAM becomes available

The VCF-only result remains reusable. A BAM-aware layer should add local
realignment, read inspection, callable-region evidence, CNV/SV read support and
variant re-calling. It should confirm or reject candidates rather than silently
changing the consequence or inheritance rules implemented here.
