# Reprendre une analyse génomique

Ce document est un point d'entrée générique pour les futurs agents. Il ne contient ni résultat individuel ni inventaire d'un couple. Les [API détaillées](../../tools/jouvence_genomics/README.md) et le [guide de qualité](report.md) restent les références pour le calcul et l'interprétation. Les exemples Python sont exécutés depuis la racine du dépôt ; les liens documentaires suivent son arborescence.

## 1. Vérifier les entrées avant de calculer

Utiliser `uv` et les versions enregistrées dans le run ou les fonctions de préparation. Vérifier les tests des modules modifiés ; ne pas actualiser toutes les références pour une simple reconstruction du rendu.

Pour chaque personne, conserver un dossier de run distinct sous `analysis/` avec identité d'échantillon, assemblage, fichiers source, SHA-256, références/licences, paramètres et versions du code exécuté. Un changement d'entrée, de référence ou de méthode crée une nouvelle révision ; il ne réétiquette pas les anciens résultats.

Les `results.json` du moteur actuel enregistrent les SHA des VCF de petites variations, mais seulement les chemins et tailles de certaines références ; ils ne figent pas toutes les empreintes des références ni du code exécuté. Un manifeste complémentaire est donc requis pour les futurs runs : entrées, toutes les références utilisées, code exécuté (SHA des fichiers et commit), dépendances, paramètres et sorties avec leurs SHA. Pour un ancien run, marquer les preuves manquantes ; une empreinte calculée aujourd'hui ne prouve pas quels octets avaient été utilisés hier.

Le chemin principal actuel accepte un VCF **GRCh37**, mono-échantillon, combinant SNP et indels. Des fichiers SNP/indel séparés passent par `analyze_without_bam()` et son argumentaire documenté ; ne pas les concaténer sans contrôle. `analyze_sample()` accepte aussi `cnv_vcf` et `sv_vcf`, facultatifs. Les répétitions et résumés de profondeur ont une entrée dédiée dans le CLI, pas un paramètre inventé pour cette fonction.

Un site absent d'un VCF de variants est **inconnu**, jamais `0/0`. Un gVCF ou un VCF tous-sites doit documenter ses appels de référence et sa couverture. Un BAM/CRAM n'est pas nécessaire pour ce chemin ; il permettrait des contrôles de lectures et un rappel ultérieur, mais aucun rappel BAM n'est exécuté par `analyze_sample()`. Un liftOver ne rappelle pas les variants et ne récupère pas les lectures perdues. Ne pas passer un futur VCF GRCh38 à un adaptateur GRCh37 en changeant son libellé.

## 2. Choisir le périmètre public avant les génotypes

Préparer ou vérifier les références publiques une fois, dans un cache séparé. Garder URL, assemblage, version, taille, SHA et checksum publié lorsqu'il existe. Inventorié, téléchargé, audité et qualifié sont quatre états différents.

| Besoin | Point d'entrée | Limite à conserver |
|---|---|---|
| Analyse individuelle sans BAM | `analyze_sample()` / `analyze_without_bam()` | ClinVar et candidats moléculaires sont des niveaux distincts ; les références requises ne sont pas toutes installées par une seule fonction |
| Comparaison de couple | `compare_couple()` | Même maladie et allèles/haplotypes pertinents, pas seulement même gène ; phase inconnue et VUS sans chiffre de maladie |
| PGx, nutrition et traits ciblés | `audit_traits()` | Génotype explicite et qualité ; pas de diplotype complet ou posologie déduits de sites manquants |
| GWAS tous thèmes | `full_gwas.py` | Associations/studies/traits/sites séparés ; aucun effet additionné en score ou diagnostic |
| PGS tous thèmes | `pgs_inventory.py` puis `audit_pgs_inventory()` | Budgets explicites, modèles choisis avant les résultats ; sous-somme vérifiée ≠ PGS qualifié |
| Annotation variant/gène | `annotation_resources.py` | Runtime versionné, contrôles des logs/UID ; un annotateur inventorié n'est pas nécessairement exécuté |
| Prédictions missense/indels | `personal_predictors.py`, `cadd_predictions.py` | Jointures exactes, qualité d'origine conservée ; pas de nouvelle classification clinique |
| Ascendance | `ancestry.py` et loaders de référence | Validation indépendante et même modèle PCA ; affinité de référence, pas ethnicité/nationalité |

Les signatures, dépendances et exemples complets sont dans le README. Ne pas sélectionner seulement les modèles ou loci positifs chez une personne pour prétendre à un panel exhaustif.

### Compter le catalogue, pas promettre une détection complète

Avant de publier un nombre de conditions, gènes ou variants, recalculer les dénominateurs depuis le snapshot public exact **et les filtres du code exécuté**, avec leurs empreintes. Dédupliquer les identifiants de maladie, gènes et allèles normalisés séparément ; distinguer les modes de transmission, les types de variants, les associations et les assertions multiples. Conserver la formule, les exclusions et la table qui permettent de retrouver chaque compteur, indépendamment des résultats positifs d'une personne.

Séparer catalogue inventorié, périmètre réellement parcouru, appels explicites utilisables après QC et régions non évaluables. L'absence d'un site dans un VCF de variants n'établit ni couverture ni résultat négatif. Un catalogue de maladies n'est pas un test cliniquement validé pour chacune d'elles ; les limites CNV/SV, répétitions, phase et régions complexes restent explicites.

### Deux profils d'usage, pas deux critères cliniques

Le profil réutilisable conserve sa politique restrictive. `license_profile="personal_noncommercial"` pour les PGS et `usage_profile="personal_noncommercial"` pour les prédicteurs sont des options **explicites**, réservées à la recherche personnelle privée réellement non commerciale. L'absence de vente ne suffit pas ; la recherche produit commerciale n'est pas couverte.

Consulter le [registre des conditions](../../tools/jouvence_genomics/resource_use_review.json). Distinguer code, modèles et tables : une licence du logiciel ne licencie pas toutes ses données. Les termes sur mesure/académiques seulement ou inconnus restent à examiner ; aucun accord d'accès n'est accepté automatiquement. Les notices contradictoires sont conservées. Une ressource en attente d'accès reste différée, pas analysée.

## 3. Figement et rendu : réutiliser les fonctions existantes

Exemple minimal pour un VCF combiné, avec références locales déjà préparées :

```python
from pathlib import Path
from tools.analyze_dante_vcf import analyze_sample
from tools.jouvence_genomics.couple import compare_couple

first = analyze_sample(
    Path("private/parent_a.vcf.gz"), Path("data/reference"),
    Path("analysis/run/parent_a"), sample_name="parent_a",
)
second = analyze_sample(
    Path("private/parent_b.vcf.gz"), Path("data/reference"),
    Path("analysis/run/parent_b"), sample_name="parent_b",
)
couple = compare_couple(first, second)  # Dictionnaire à figer dans le run privé.
```

L'analyse enregistre `results.json`. Conserver aussi le dictionnaire du couple en JSON privé, sans modifier ses observations. Pour reconstruire uniquement le rendu, fournir ces trois JSON à `build_genomic_reports()` ; ne pas refaire les génotypes ou les annotations de recherche pour cette étape.

```python
from tools.jouvence_genomics.html_reporting import build_genomic_reports

reports = build_genomic_reports(
    Path("analysis/run/parent_a/results.json"),
    Path("analysis/run/parent_b/results.json"),
    Path("analysis/run/couple/results.json"),
    Path("analysis/run/rapports"),
    skill_scripts=Path("/chemin/local/vers/html-report/scripts"),
)
```

`refresh_traits=False` ne réutilise un cache que si les empreintes des entrées, du panel et du décodeur concordent ; sans cache, le panel est calculé. Les résultats facultatifs `ancestry_results`, `polygenic_results` et `supplementary_results` sont indexés par le `sample_name` exact et leurs SHA doivent correspondre. Ne pas attribuer un fichier par son seul nom ou voisinage dans un dossier.

Le texte d'interprétation a sa propre empreinte (`interpretation_definition_sha256`)
et une révision `report_version`, distinctes du run génomique. Modifier ces
textes n'invalide pas les génotypes du panel. En revanche, corriger le décodeur
(par exemple FT ou ploïdie) exige une réextraction du panel ciblé avec
`refresh_traits=True` ; ne pas présenter l'ancien cache comme vérifié par la
nouvelle règle. Ce n'est pas une relance ClinVar/pLoF/PGS.

Pour chaque association utile, utiliser la structure commune « résultat /
signification / ampleur publiée et comparateur / action / limite ». Les
coefficients du GWAS Catalog ne sont pas automatiquement des OR. Les pQTL
de protéines ne prédisent pas un chronotype ; une liste de SNP associés
n'est pas une somme de risques. Réexporter les audits complets avec
`export_association_audit()` pour ajouter ces explications aux CSV, en
conservant leur identité, tous les génotypes inconnus et les valeurs source.

### Ce qui est affiché et ce qui reste dans l'audit

- Une observation peut être intéressante sans constituer un diagnostic. Les fiches de maladie et suites possibles restent conditionnelles ; aucun traitement ou PGT décidé à partir d'un VUS ou d'une prédiction seule.
- Les VUS et appels techniques peuvent être inspectés dans l'annexe, jamais comme des cartes de maladie ou des portages établis. Les inconnus ne deviennent pas « rien détecté ».
- Les valeurs importées du fournisseur non reproduites/qualifiées et les comparaisons fournisseur/VCF restent hors HTML patient, panels, annexes, CSV et pièces jointes. Conserver leurs audits privés séparément, sans modifier les conclusions locales pour les faire correspondre.
- Les PGS incomplets, sous-sommes, bornes non qualifiées et prédictions personnelles AlphaMissense/REVEL/CADD restent des audits de recherche séparés. Ils ne sont pas ajoutés automatiquement au rapport.
- La synthèse d'un PGS déjà qualifié applique `abs(z) > 1`, ou des bornes entièrement hors de ±1 sigma, dans les deux queues. Ce seuil ne qualifie pas le modèle ; un percentile n'est pas un risque absolu. Aucune extrapolation depuis quelques SNP ou une capture.
- Les extensions exhaustives GWAS/OpenCRAVAT/PGS ne sont pas intégrées au rapport simplement parce qu'un nouveau JSON existe. Vérifier explicitement les types et critères que le renderer accepte.

### Panels et CSV : annoncer le bon dénominateur

Séparer variants/allèles uniques, gènes, maladies, associations, études, transcrits et lignes exportées. Un aperçu progressif n'est pas le CSV complet et un CSV d'observations n'est pas tous les sites du modèle.

Les CSV de panels reflètent uniquement les contenus admis dans le rapport. Les CSV de recherche restent privés et distincts : AlphaMissense/REVEL `*.matched.csv` contient les SNV appariés, par score/transcrit/observation ; le JSONL conserve aussi les inconnus/exclus. Le CSV CADD conserve les indels observés du panel, appariés ou non, avec toutes les lignes publiques correspondantes. Des scores CADD contradictoires restent inconnus, sans moyenne ni choix arbitraire.

Avant remise : vérifier SHA/compteurs/CSV, critères de qualification, liens et sources, navigation/recherche/détails, mobile et absence de transmission réseau. La vérification d'un rendu ou d'une jointure ne vaut pas reproduction du séquençage. L'envoi à une personne nécessite un destinataire et les pièces autorisés explicitement ; ne jamais envoyer le dossier d'audits entier par défaut.

## 4. Conditionner et livrer les rapports

Un HTML peut être lisible seul hors ligne tout en contenant des liens vers des CSV annexes. Pour partager le rapport et ses exports statiques, utiliser la fonction existante avec les HTML choisis explicitement et un dossier de sortie absent ou vide. Dans la suite de l'exemple précédent, `reports` contient les chemins retournés par le renderer :

```python
from pathlib import Path
from tools.jouvence_genomics.report_delivery import prepare_report_delivery

delivery_manifest = prepare_report_delivery(
    {
        "Parent A": reports["personne_1"],
        "Couple": reports["couple"],
    },
    Path("output/livraison_run"),
)
```

La fonction crée un ZIP par rapport avec ses annexes statiques référencées, une notice et un manifeste d'intégrité ; elle ne collecte aucun audit par voisinage et n'envoie rien. Décompresser chaque ZIP dans un **dossier distinct**, puis ouvrir son HTML en conservant l'arborescence des annexes. Vérifier les listes, CRC et SHA du manifeste avant partage. Les personnes dont les rapports sont inclus et le destinataire doivent être autorisés explicitement ; ne pas envoyer les analyses privées complètes par défaut.

## 5. Sauvegarder sans effacer le Mac

Une sauvegarde privée autorisée copie des données ; elle n'autorise ni publication, ni inférence distante, ni suppression locale. La destination exacte, le projet de facturation et les sources sont des paramètres vérifiés, pas des valeurs reprises d'un ancien mémo. Ne pas utiliser un préfixe existant du KB pour un nouveau backup.

API de sauvegarde :

```python
from pathlib import Path
from tools.jouvence_genomics.cloud_backup import (
    prepare_genomics_backup, upload_genomics_backup, verify_genomics_backup,
)

def backup_verified_sources(
    sources: dict[str, Path], work_directory: Path,
    destination: str, billing_project: str,
) -> Path:
    """Copier les sources autorisées ; retourner le chemin de verified.json.

    Les labels sont les chemins relatifs de restauration. Le dossier de travail
    est vide et hors des sources ; la destination est un nouveau préfixe GCS.
    Aucun original n'est modifié ou supprimé.
    """
    manifest = prepare_genomics_backup(sources, work_directory)
    receipt = upload_genomics_backup(
        manifest, destination, billing_project=billing_project,
    )
    return verify_genomics_backup(
        manifest, receipt, billing_project=billing_project, read_back=True,
    )
```

`sources` est un dictionnaire de chemins relatifs de restauration et de sources locales explicites : références, analyses et VCF externes seulement s'ils sont autorisés. Inclure les cibles des liens symboliques ; une dépendance externe non incluse arrête la préparation. Choisir `work_directory` vide et hors des sources archivées. Fermer les écrivains SQLite et toute génération d'analyse avant l'archivage ; les fichiers WAL/SHM restent inclus. Le runtime utilise `zstandard`, `google-crc32c` via `uv` et le client `gcloud` authentifié ; figer leurs versions dans la provenance d'exécution. Vérifier en direct la confidentialité et la politique de chiffrement du bucket ; ne pas supposer qu'une clé CSEK est requise ou permise.

Le manifeste lie les sources aux archives tar.zst, à leurs empreintes et aux objets cloud. Le reçu et `verified.json` doivent confirmer tailles, CRC32C, générations et relecture complète avec SHA des archives ; un upload lancé n'est pas une sauvegarde vérifiée. Conserver tous les originaux locaux. Les liens symboliques relatifs sont préparés dans les copies, pas en modifiant les originaux.

Après restauration, les chemins absolus des anciens runtimes/manifests peuvent être invalides. Conserver la provenance historique et revalider ou reconstruire les chemins locaux avec les fonctions `prepare_*` adaptées ; vérifier les entrées/exports avant de recalculer. Un cache restauré n'est pas encore un runtime opérationnel, et le backup ne prouve ni l'interface ni les calculs cliniques reproduits.
