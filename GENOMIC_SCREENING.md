# Jouvence — exploratory genomic screening

The analysis combines three distinct evidence layers:

- ClinVar GRCh37 for reviewed germline pathogenic/likely-pathogenic variants;
- GenCC for gene–disease validity and mode of inheritance;
- published gnomAD v4.1 carrier frequencies for autosomal-recessive genes.

No single database is sufficient: ClinVar is variant-level, GenCC is
gene–disease-level, and gnomAD is population-level. OMIM is broader for
Mendelian phenotypes but its bulk data are licensed.

The code is split into small, documented modules:

- `vcf.py`: VCF parsing, QC, genotype and phase fields;
- `catalogs.py`: ClinVar, GenCC, GENCODE and carrier-frequency loaders;
- `analysis.py`: matching, inheritance, cis/trans logic and CNV/SV annotation;
- `reporting.py`: reproductive-risk estimates and Markdown/JSON output.

`tools/analyze_dante_vcf.py` is only the thin command-line entry point. Run it
with `uv`:

```bash
UV_CACHE_DIR=/tmp/jouvence_uv_cache uv run python tools/analyze_dante_vcf.py \
  --snp-vcf /path/to/sample.snp.vcf.zip \
  --indel-vcf /path/to/sample.indel.vcf.gz \
  --cnv-vcf /path/to/sample.cnv.vcf \
  --sv-vcf /path/to/sample.sv.vcf \
  --depth-summary /path/to/depth_summary.txt \
  --clinvar-small data/reference/clinvar_grch37.vcf.gz \
  --clinvar-structural data/reference/clinvar_variant_summary.txt.gz \
  --gencc data/reference/gencc_submissions.tsv \
  --gencode data/reference/gencode.v19.annotation.gtf.gz \
  --carrier-frequencies data/reference/ar_gene_carrier_frequencies_gnomad_v4.zip \
  --output analysis/dante
```

The workflow is local: it does not send genes or genotypes to an external API.

Phase is used only when the VCF contains a phased genotype (`|`) and a shared
`PID` or `PS` block. Two heterozygous variants with a common recessive disease
identifier are then classified as cis or trans. Otherwise their phase remains
unknown. Two different variants are never labelled homozygous.

CNV/SV calls are compared with GRCh37 genes and reviewed P/LP ClinVar intervals.
A pathogenic interval match requires the same event type and at least 80%
reciprocal overlap. Gene overlap alone is retained for technical review but is
not called pathogenic. Dante's files lack read support, confidence intervals,
caller metadata and useful quality scores, so every structural signal requires
confirmation from BAM/CRAM or FASTQ.

Repeat expansions require a supported caller output or BAM/CRAM/FASTQ. They
cannot be inferred reliably from the four variant-only Dante VCFs.

The output is research-only. A candidate must be checked in the BAM/CRAM and
confirmed in an accredited clinical laboratory before clinical or reproductive
use. The script does not assess PGT-M eligibility by itself.
