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

## Rapports HTML français

Le [guide de qualité](../../docs/rapport_genomique/report.md) définit le contenu
attendu des rapports individuels et du couple. Sa version HTML autonome est
`docs/rapport_genomique/guide.html`.

`build_genomic_reports()` réutilise les JSON figés, vérifie les SHA-256 des VCF,
enrichit les observations avec les fiches cliniques versionnées, extrait le
panel ciblé PGx/nutrition/traits et génère trois HTML autonomes. Les appels
inconnus restent inconnus ; aucun diplotype, posologie, groupe ABO définitif
ou ascendance n'est fabriqué. Le filtre de navigation ne change aucun calcul.

```python
from pathlib import Path
from tools.jouvence_genomics.html_reporting import build_genomic_reports

reports = build_genomic_reports(
    Path("analysis/parent_a/results.json"),
    Path("analysis/parent_b/results.json"),
    Path("analysis/couple/results.json"),
    Path("analysis/rapports_html"),
    skill_scripts=Path("/chemin/vers/html-report/scripts"),
    report_date="2026-10-10",
)
```

Le renderer provient du guide HTML local et expose `render_report`. Utiliser
`uv run --with Markdown==3.8.2 --with beautifulsoup4==4.13.4 --with plotly==6.3.0`
pour son environnement ; aucun runtime Plotly n'est embarqué lorsque le rapport
ne contient pas de graphiques. Les détails et filtres utilisent du JavaScript
local. Les HTML nominaux, caches, tableaux et captures restent sous `analysis/`
ou `output/`, exclus de Git.

Les connaissances génériques `clinical_*_cards.json` et `trait_markers.json`
contiennent les sources et doivent être revues lorsqu'une recommandation ou
classification évolue. Les cartes décrivent une maladie et sa prise en charge
conditionnelle ; elles n'attribuent pas ce diagnostic à un porteur ou un VUS.

Le panel de traits vérifie DP/GQ, AD, compatibilité de référence et FILTER.
Les homozygotes doivent être cohérents avec les lectures. CYP3A5 *3 et un
allèle ABO O courant correspondent ici à la référence : un site absent ne
permet donc aucune conclusion à leur sujet. Les gènes PGx complexes, les allèles
HLA, RhD et l'ascendance nécessitent des méthodes supplémentaires.

### Première lecture et annexe

Le rendu v3 sépare `patient_view.py` (sélection éditoriale et cartes simples)
de l'analyse et de l'annexe exhaustive. Les VUS et les sites inconnus ne sont
pas des cartes de maladie. Les marqueurs de qualité insuffisante sont écartés
de la synthèse ; `FILTER=.` reste visible comme observation provisoire, même
quand les autres contrôles disponibles sont satisfaisants. Le panel comprend
57 sites et 27 profils sourcés, pas 57 diagnostics ou 27 diplotypes.

### Ascendance et scores polygéniques : calculs facultatifs et explicites

Ces fonctions ne requièrent pas de BAM et ne transmettent aucun génome. Leurs
téléchargements portent uniquement sur des références publiques. Les caches,
génotypes extraits et résultats doivent rester dans les répertoires ignorés.
Utiliser `uv run --with numpy==2.5.3 --with scikit-learn==1.9.1` pour le calcul
d'ascendance, en complément des dépendances du rendu.

`ancestry.py` fournit `prepare_ancestry_reference()` et
`analyze_ancestry()`. Le panel Peddy/1000 Genomes GRCh37 est épinglé par commit
et SHA-256. L'entraînement et les références indépendantes sont séparés ; les
références sont masquées pour reproduire l'absence de 0/0 dans un VCF de variants.
Les contrôles incluent plis indépendants, sous-panels et chromosomes pairs/impairs.
La sortie est une affinité à cinq groupes de référence, jamais un pourcentage
de mélange, une ethnicité ou une nationalité. Le nuage PCA et le résultat
doivent porter la même empreinte de modèle, pas seulement les mêmes fichiers.

`polygenic.py` fournit `prepare_polygenic_references()`,
`audit_polygenic_scores()` et `calibrate_polygenic_bounds()`. Les poids PGS
Catalog GRCh37, métadonnées, licences et empreintes sont conservés. Tout dosage
absent reste inconnu ; des bornes explicitent cette incertitude. Les loci liés
à l'X exigent une convention de ploïdie validée. Une référence 1000 Genomes
peut calibrer un sous-ensemble de recherche, sans le transformer en PGS clinique
validé ou en risque absolu. Aucun percentile n'est affiché si la qualification
échoue. La règle de synthèse est strictement `abs(z) > 1`, ou des bornes
entièrement au-delà de ±1 ; elle s'applique aux deux queues.

Transmettre les JSON à `build_genomic_reports()` via `ancestry_results` et
`polygenic_results`, indexés par le `sample_name` exact, et le nuage public via
`ancestry_reference_projection`. Le constructeur vérifie les empreintes des
VCF et de l'analyse pour empêcher une attribution à une autre personne. Les
JSON complets restent téléchargeables dans l'HTML autonome.
