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

## Combined VCF and couple analysis

For one SNP/indel VCF, call the Python function directly:

```python
from pathlib import Path
from tools.analyze_dante_vcf import analyze_sample
from tools.jouvence_genomics.couple import compare_couple, render_couple_report

sample = analyze_sample(
    Path("sample.vcf.gz"), Path("data/reference"), Path("analysis/sample"),
    sample_name="Sample",
)
couple = compare_couple(first_parent_result, sample)
report = render_couple_report(couple)
```

The equivalent CLI uses `--vcf sample.vcf.gz --sample-name Sample`; do not
also pass `--snp-vcf` or `--indel-vcf`. CNV/SV inputs are optional. A single
sample and GRCh37 contig lengths are checked before analysis. Inputs have
SHA-256 hashes and reference filenames in the JSON provenance.

Couple screening compares recessive disease identifiers as well as genes.
Different alleles can form a relevant pair. VUS pairs are retained for review
but receive no disease-risk estimate. A conditional 25% biallelic transmission
requires one quality-passing P/LP carrier haplotype in each parent and a
shared disease. Unknown phase leaves transmission unresolved. Penetrance,
allele-specific severity and difficult genomic regions require clinical review.
`inspect_partner_loci()` records explicit counterpart calls or unassessed
absence; it does not infer reliable homozygous-reference genotypes.

The novel coding screen currently covers autosomal-recessive genes using one
selected transcript per gene. Novel dominant/X-linked, deep intronic and
regulatory variants are outside its scope. GenCC/ClinVar snapshots, gnomAD
v2.1.1 exome frequencies and v4.1 gene-carrier frequencies remain distinct
reference sources. NFE is an ancestry proxy, not a French population panel.

For raw (`FILTER=.`) files, `audit_candidate_quality()` attaches INFO metrics
and checks available SNP/indel annotations against the Broad GATK hard-filter
recommendations. Failed candidate calls are retained for technical review,
and do not enter inheritance or couple risk calculations. Missing metrics
remain explicit. Nearby indels (within 20 bp) and indels >=50 bp receive a
separate complex-region tier requiring joint haplotype/read validation;
these flags do not prove that a call is false or benign.
