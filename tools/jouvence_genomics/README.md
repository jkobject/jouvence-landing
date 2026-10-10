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
58 sites et 28 profils sourcés, pas 58 diagnostics ou 28 diplotypes.
ACTN3 est un marqueur de recherche sur les fibres musculaires, pas un test
d'endurance individuelle. Le REF GRCh37 de rs1815739 est T ; ne pas lui appliquer
le REF C de GRCh38.

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
JSON volumineux restent dans le dossier privé de l'HTML ; les petits audits
et les CSV des panels sont embarqués. Conserver le dossier complet pour partager
un rapport avec toutes ses pièces, sans publier ces données dans Git.

### Imports Dante et panels v4

L'import compare des **sources**, il ne certifie pas les conclusions commerciales.
`read_dante_html()` lit les tableaux des méthodes, pas les compteurs animés.
`extract_dante_pdf_document()` conserve toutes les pages, les liens, les textes
et leur SHA-256 ; `compare_dante_pdf_genotypes()` confronte les allèles explicites
aux VCF à des positions GRCh37 vérifiées. Les seuls appels réseau concernent
des rsID et références publics. Les génomes restent locaux.

```python
from tools.jouvence_genomics.dante_html import read_dante_html, compare_reported_genotypes
from tools.jouvence_genomics.document_import import extract_dante_pdf_document, compare_dante_pdf_genotypes

external = read_dante_html([Path("private/report.html")], expected_profile="Profil attendu")
comparison = compare_reported_genotypes(external, [Path("private/sample.vcf.gz")])
document = extract_dante_pdf_document(Path("private/sleep.pdf"), Path("analysis/imports"))
pdf_comparison = compare_dante_pdf_genotypes(
    [document], [Path("private/sample.vcf.gz")], Path("data/reference/pdf_coordinates"),
)
```

Pour les PDF : `uv run --with pypdf==6.20.0 --with pdfplumber==0.11.10`.
Un résultat absent, un profil étoilé non reconstruit ou une rubrique « Not
currently supported » ne deviennent jamais « normal » ou « négatif ». Le brin
des PDF n'est pas attesté : les concordances sont littérales, avec contrôle de
compatibilité au complément en cas de divergence, pas une validation clinique.
`provider_topics.json` et `provider_result_labels.json` sont des textes français
génériques, pas des génotypes ni résultats personnels.

Les comparaisons fournisseur/VCF sont conservées dans `dante_comparison.json`
et `pdf_comparison.json`, à part du rapport personnel. Elles ne sont ni rendues
dans l'HTML, ni embarquées dans ses pièces jointes, ni ajoutées aux CSV des panels.
Les génotypes du fournisseur restent attribués à ce fournisseur ; les critères
de qualité et conclusions des analyses locales ne sont pas modifiés.

Le rapport personnel n'affiche aucun résultat importé non reproduit
indépendamment, même si le document fournit un percentile hors ±1 sigma.
Les PDF/HTML, les scores incomplets et leurs blocages restent dans l'audit
privé séparé, pas dans les panels, CSV ou pièces jointes du rapport personnel.
Une somme recalculée sur quelques variants n'est pas un PGS complet ; la
reproduction arithmétique et la qualification d'un résultat sont deux contrôles
distincts.

`supplementary_results` de `build_genomic_reports()` associe à chaque nom
d'échantillon les chemins JSON `dante_html`, `dante_comparison`, `pdf_reports`
(liste), `pdf_comparison`, `associations` et `exact_polygenic`. Les SHA-256
des VCF et documents sont contrôlés avant assemblage. Ne jamais associer un
document externe à une personne par la simple proximité de son nom de fichier.

Le dashboard offre filtres de catégories, recherche, explications françaises,
animations respectant `prefers-reduced-motion`, tableaux progressifs et CSV UTF-8.
Les champs supplémentaires restent dans le CSV ; les cellules pouvant être des
formules sont neutralisées. Les scores non qualifiés restent entièrement hors
rapport ; leurs audits se consultent séparément. Pour les gros modèles, un aperçu limité est séparé
du CSV complet des observations explicites ; les sites inconnus restent dans
l'audit du modèle. Aucun dosage manquant n'est complété.

### Essais de reproduction indépendants, sans BAM

`reproduction.py` fournit trois fonctions publiques : `evaluate_site_filter()`,
`audit_pgs_reproduction()` et `verify_subset_with_plink()`. Elles distinguent
explicitement les appels observés, les dosages inconnus et la vérification
arithmétique. Le filtre original reste la règle par défaut ; l'expérience
`gatk_hard_filter_2025` applique des seuils publiés aux annotations disponibles,
sans modifier les VCF ni les critères du rapport clinique.

```python
from tools.jouvence_genomics.reproduction import audit_pgs_reproduction, verify_subset_with_plink

audit = audit_pgs_reproduction(
    [Path("private/sample.vcf.gz")],
    [Path("data/reference/PGS000001_hmPOS_GRCh37.txt.gz")],
    output_directory=Path("analysis/reproduction"),
)
validation = verify_subset_with_plink(
    Path("analysis/reproduction/PGS000001.json"),
    Path("private/tools/plink2"), Path("analysis/reproduction/plink"),
)
```

Une sous-somme vérifiée par PLINK n'est ni un PGS complet, ni un percentile,
ni une reproduction du chiffre commercial. Les sites inconnus restent bornés,
jamais imputés comme référence. Un score exploitable demande ensuite une
calibration et une validation prédictive adaptées, distinctes de ce test.

`pharmacogenomic_reproduction.py` fournit `prepare_pharmcat_runtime()` et
`reproduce_pharmacogenomics()`. Le runtime public PharmCAT/JDK est épinglé par
SHA-256 et isolé dans le cache ; cet adaptateur est actuellement macOS arm64.
Il projette seulement des SNP GRCh37 existants vers GRCh38, avec mapping unique,
réciproque et contrôle des références publiques. Il exécute le matcher officiel
en mode non phasé, sans reporter clinique ni recommandation de posologie.
Cet essai accepte seulement les GT diploïdes. Les échecs FORMAT/FT sont
exclus ; les filtres individuels absents/non évalués restent explicitement signalés.

```python
from tools.jouvence_genomics.pharmacogenomic_reproduction import reproduce_pharmacogenomics

pgx_audit = reproduce_pharmacogenomics(
    Path("analysis/sample/results.json"),
    Path("data/reference/pharmcat"), Path("analysis/sample/pharmcat_audit"),
)
```

Utiliser `uv run --with pyliftover==0.4.1` pour cet adaptateur. Les références
publiques sont préparées avant la lecture des génotypes ; aucun génome n'est
transmis. Les SNP manquants, indels, HLA et CNV restent inconnus. Un diplotype
compatible avec des appels incomplets n'est pas un diplotype confirmé.

Ces deux essais acceptent actuellement des entrées **GRCh37** vérifiées. Pour
un futur VCF GRCh38, le moteur et les poids devront correspondre à cet assemblage.
Pour un VCF GRCh37, les poids PGS Catalog GRCh37 évitent un liftOver inutile.
Un liftOver ne remappe pas les lectures, ne rappelle aucun génotype et ne
récupère pas les données perdues. Même un nouveau VCF GRCh38 de variants seuls
ne documente pas ses positions absentes : demander si possible un gVCF avec
blocs de référence et métriques de qualité, ou un VCF tous-sites réellement
génotypé aux loci ciblés. Un BAM/CRAM indexé permettrait un rappel ultérieur.

### Catalogue d'associations et références d'ascendance alternatives

`association_catalog.py` télécharge le GWAS Catalog versionné, sélectionne des
domaines annoncés, remappe les rsID GRCh38 vers GRCh37 et lit les VCF une seule
fois. Chaque ligne garde PMID, étude GCST, phénotype, p-value, orientation,
effet publié et qualité de l'appel. Une association n'est pas une mutation
mendélienne, un résultat causal ou un score : ce module n'additionne aucun effet.
Les associations, sites, études et phénotypes doivent être comptés séparément.
Les codes et métadonnées de référence sont conservés pour vérifier la fraîcheur.

`analyze_ancestry(..., reference_loader=...)` accepte une référence alternative
déjà préparée avec matrice 0/1/2/NaN, sites GRCh37, labels, manifeste et sources.
Les contrôles par données indépendantes, masques variant-only, sous-panels et
chromosomes gardent les mêmes critères. Aucun groupe n'est déduit d'une nationalité,
et une référence plus grande qui échoue aux contrôles ne remplace pas un résultat
précédent par une conclusion plus précise non étayée.

`hgdp_reference.py` fournit `prepare_hgdp_reference()` (à exécuter une fois)
et `load_hgdp_reference()` (cache compact réutilisable). Les requêtes ciblent
seulement les SNP publics, par plages BCF indexées ; BCF/CSI, chaîne d'assemblage,
métadonnées et parenté sont liés par générations/empreintes. Les apparentés des
auteurs sont exclus ; groupes larges et noms de populations fines restent séparés.
Utiliser `uv run --with numpy==2.5.3 --with scikit-learn==1.9.1
--with pysam==0.24.1 --with pyliftover==0.4.1 --with zstandard==0.25.0`.

```python
from tools.jouvence_genomics.hgdp_reference import prepare_hgdp_reference, load_hgdp_reference
from tools.jouvence_genomics.ancestry import analyze_ancestry

prepare_hgdp_reference(Path("data/reference/ancestry_hgdp"), Path("data/reference/ancestry"))
outputs = analyze_ancestry(
    [Path("analysis/sample/results.json")], Path("data/reference/ancestry_hgdp"),
    Path("analysis/sample/ancestry"), reference_loader=load_hgdp_reference,
)
```

La diversité géographique ne se résume pas au nombre de personnes. Les très grands
panels français ou de biobanques peuvent demander un accord d'accès ; ils ne sont
pas remplacés par une prétendue nationalité tirée d'un petit groupe French/Basque.

### Extension tous thèmes, indépendante des résultats personnels

Les fonctions ci-dessous sont **opt-in** : elles élargissent l'audit sans changer
les seuils cliniques, le petit panel de traits ou les rapports déjà produits.
Le périmètre public est défini avant lecture des génotypes. « Inventorié »,
« référence préparée », « génotype observé », « calcul effectué » et « résultat
qualifié » sont des états différents. Un libellé GWAS, une colonne d'annotateur
ou un modèle PGS n'est pas un test clinique indépendant.

#### GWAS Catalog complet

`full_gwas.py` conserve toutes les lignes du snapshot, sans restriction aux
sept thèmes historiques. Un seuil de publication `p <= 5e-8` sélectionne les
loci à rechercher ; les lignes exclues restent exportées avec leur raison.
Les coordonnées sont contrôlées en GRCh37 par rsID, chaîne GRCh38→37 et REF
public. La référence dbSNP151Common est ancienne et limitée aux SNP communs :
les rsID non résolus restent inconnus, pas négatifs.

```python
from pathlib import Path
from tools.jouvence_genomics.association_catalog import download_gwas_catalog
from tools.jouvence_genomics.full_gwas import (
    prepare_exhaustive_gwas, prepare_full_gwas_resources,
    map_exhaustive_gwas, audit_exhaustive_gwas,
)

download = download_gwas_catalog(Path("data/reference/gwas_download"))
catalogue = prepare_exhaustive_gwas(download, Path("data/reference/gwas_all"))
resources = prepare_full_gwas_resources(Path("data/reference/gwas_mapping"))
catalogue = map_exhaustive_gwas(catalogue, resources)
audit = audit_exhaustive_gwas(
    [Path("private/sample.vcf.gz")], catalogue,
    Path("analysis/sample/gwas_all"), sample_name="Sample",
)
```

Sorties privées : résumé JSON, observations SQLite et CSV gzip complet, avec
les colonnes source, la qualité et les états inconnus/exclus. Les formules de
tableur sont neutralisées dans le CSV. Les p-values sont celles des études,
pas des tests réalisés sur la personne. Aucun effet n'est additionné ou traduit
en risque personnel. Un dossier d'audit existant n'est pas écrasé.
Utiliser `uv run --with biopython --with pyliftover`.

#### Tous les modèles du PGS Catalog, calcul à budget explicite

`pgs_inventory.py` fige les sept tables du catalogue et confronte leurs IDs à
la liste officielle complète. Le téléchargement des poids est borné par trois
paramètres obligatoires, sur tous les thèmes, par taille puis ID et jamais
d'après un résultat personnel. Licences restrictives/inconnues, interactions,
formats non pris en charge et discordances sont différés avec leur raison.
Une licence restrictive dans les métadonnées **ou** le fichier de poids est
éliminatoire ; un ancien header sans licence peut utiliser les termes explicites
des métadonnées. Ce filtre prudent ne remplace pas une revue juridique.

```python
from pathlib import Path
from tools.jouvence_genomics.pgs_inventory import (
    prepare_pgs_inventory, prepare_bounded_pgs_weights, audit_pgs_inventory,
)

inventory = prepare_pgs_inventory(Path("data/reference/pgs_inventory"))
selection = prepare_bounded_pgs_weights(
    inventory, Path("data/reference/pgs_inventory/weights"),
    max_variants_per_model=5_000, max_total_variants=2_000_000,
    max_download_bytes=1_000_000_000,
)
audit = audit_pgs_inventory(
    [Path("analysis/sample/results.json")], selection,
    Path("analysis/sample/pgs_all_strict"),
)
```

Les poids et VCF sont vérifiés par SHA avant calcul. Le filtre original reste
la règle par défaut. `site_filter_policy="gatk_hard_filter_2025"` est un essai
séparé, non une qualification clinique. Le moteur actuel accepte GRCh37 et
les modèles autosomiques additifs simples. Les dosages absents restent inconnus ;
seuls des sous-sommes/bornes et éventuels scores bruts complets sont produits,
jamais un percentile, z-score, risque absolu ou chiffre importé du fournisseur.
Augmenter les budgets est possible, mais augmente surtout la mémoire et les
sorties ; la liste des modèles différés permet de décider explicitement du
prochain lot. Un nouveau dossier évite de mélanger personnes ou configurations.

#### OpenCRAVAT local : annotations variant et gène

`annotation_resources.py` inventorie la boutique et ses colonnes déclarées,
puis prépare un moteur OpenCRAVAT versionné dans un cache dédié avec `uv`.
Chaque base conserve sa version et ses conditions : la licence du logiciel
n'autorise pas automatiquement toutes ses données.

Les six annotateurs de cette première politique sont ClinVar, GWAS Catalog,
PharmGKB, HPO, ClinGen et CGD. Les autres annotateurs restent inventoriés, pas
prétendument installés ou exécutés. CADD/REVEL/SpliceAI et les termes contradictoires
des fichiers AlphaMissense demandent une décision de licence distincte. Une
erreur de module arrête le run, même si OpenCRAVAT retourne un succès global.

`coding_variant_positions()` définit l'union de tous les CDS/start/stop
GENCODE19 et leurs abords ±8 bases avant lecture du VCF. Une base GWAS publique
validée peut y ajouter ses positions. Les lignes VCF sélectionnées sont copiées
intactes dans un sous-VCF privé, puis annotées par lots de 5 000 allèles maximum
via `annotate_vcf_positions()`. Le moteur prend seulement les appels diploïdes
explicitement observés, avec REF vérifiée, projection réciproque et chaînes
épinglées ; aucune absence n'est complétée. Les appels mal filtrés ou de qualité
insuffisante n'acquièrent aucune validité clinique grâce à une annotation.
Y/MT, SV symboliques et séquences de plus de 1 000 bases restent hors moteur.

`annotate_coding_variants()` orchestre ces étapes dans un dossier neuf et
conserve un résumé courant/final, les comptes et les empreintes des exports.
Un résumé `running` ne signifie jamais que le panel entier a été analysé.

```python
from pathlib import Path
from tools.jouvence_genomics.annotation_resources import (
    prepare_annotation_runtime, annotate_coding_variants,
)

runtime = prepare_annotation_runtime(Path("data/reference/opencravat"))
summary = annotate_coding_variants(
    Path("analysis/sample/results.json"),
    Path("data/reference/gencode.v19.annotation.gtf.gz"), runtime,
    Path("analysis/sample/annotations_all"),
    association_database=Path(catalogue["sqlite_path"]),  # Facultatif.
)
```

Sur macOS, le réseau est interdit pendant l'annotation privée, la télémétrie
est désactivée, et les sous-processus ne peuvent pas écrire dans la configuration
OpenCRAVAT partagée préexistante. Cet adaptateur n'est pas encore portable à
tous les hôtes. Les sorties officielles TSV/SQLite et l'audit JSON restent
privés ; ni classification ACMG automatique, ni diplotype, ni diagnostic n'est
produit. Utiliser `uv run --with pyyaml --with pyliftover --with twobitreader`.

Le champ natif GWAS `risk_allele_match` décrit le variant envoyé au moteur,
pas le dosage de la personne : OpenCRAVAT reçoit ici REF/ALT sans le GT.
Pour vérifier le portage, utiliser les observations GT et
`alternate_allele_present` du JSON audité. Ne pas convertir ce champ natif
en résultat patient ou en score ; les TSV/SQLite sont des annotations brutes.

Le moteur stock peut laisser des erreurs de mapper malgré un retour global
normal, voire compenser une ligne perdue par un identifiant dupliqué. L'adaptateur
contrôle donc stderr et l'ensemble exact des identifiants, pas seulement le
nombre de lignes. Les connexions SQLite sont explicitement fermées.

`derive_pharmgkb_runtime(original_runtime, new_directory)` est une option
explicite pour le bug de catégorie vide du module PharmGKB 2026.07.01 : version
et SHA du code fournisseur doivent correspondre exactement. Elle produit un
runtime séparé avec un correctif d'une ligne, diff, licence et empreintes ;
les données et le runtime d'origine restent inchangés. Aucun site personnel
n'est inclus dans ce runtime réutilisable. Un échec de mapper démontré peut
faire l'objet de `engine_exclusions` privées à l'appel `annotate_vcf_positions()`
avec TSV/log/SHA/version. Le même job, l'UID et le tuple normalisé doivent
correspondre exactement à INPUT/ERROR : l'allèle reste explicitement non annotable,
jamais « normal ». Un nouvel échec inattendu interrompt toujours l'analyse.

#### Référence d'ascendance AADR moderne, plus dense

`aadr_reference.py` fournit `prepare_aadr_reference()` et
`load_aadr_reference()` pour AADR v66.p1, Dataverse version 14, CC0. Les données
anciennes/pseudo-haploïdes, représentations répétées et apparentés proches sont
exclus. Les plages téléchargées portent des empreintes distinctes : elles ne
sont jamais présentées comme une vérification MD5 du fichier complet.

```python
from pathlib import Path
from tools.jouvence_genomics.aadr_reference import prepare_aadr_reference, load_aadr_reference
from tools.jouvence_genomics.ancestry import analyze_ancestry

prepare_aadr_reference(
    Path("data/reference/ancestry_aadr"), Path("data/reference/ancestry_hgdp"),
    Path("private/tools/plink2"),
)
outputs = analyze_ancestry(
    [Path("analysis/sample/results.json")], Path("data/reference/ancestry_aadr"),
    Path("analysis/sample/ancestry_aadr"), reference_loader=load_aadr_reference,
)
```

La préparation exige les références HGDP déjà vérifiées et PLINK2 explicitement
fourni. Le loader n'admet dans le classificateur que les labels publiés HGDP/1KG,
jamais ceux dérivés du seul pays. Le nombre du cache public et celui de
l'entraînement sont donc séparés. L'élagage LD, KING, le contrôle d'orientation
par génotypes publics, les validations masquées et la sensibilité 4/8 axes sont
des contrôles distincts. Une affinité de référence stable reste une affinité,
pas une décomposition d'ascendance ou une origine française.
