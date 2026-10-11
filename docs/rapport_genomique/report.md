<section id="objectif" markdown="1">

## Ce qu'un rapport génomique de qualité doit permettre

<p class="takeaway">Comprendre un résultat, juger sa solidité et savoir quelle étape vient ensuite.</p>

Le rapport doit comporter une synthèse courte, une lecture détaillée par personne,
une comparaison de couple et une annexe technique inspectable. Chaque observation
intéressante doit mener à une fiche explicative : maladie associée, manifestations,
transmission, niveau de preuve, qualité de l'appel, limites et suites possibles.
Le français est la langue de lecture ; les symboles de gènes et noms internationaux
des variants restent inchangés pour pouvoir retrouver les sources.

La première page répond à quatre questions : quels résultats sont suffisamment
documentés pour être discutés en génétique clinique ; quels signaux demandent une
validation ; quel risque de couple est établi ou reste indéterminé ; quelles données
manquent pour aller plus loin. Une absence de résultat est accompagnée du périmètre
réellement testé, jamais présentée comme une garantie de santé.

<p class="caveat">La prédiction d'un effet moléculaire, la classification d'un variant et le diagnostic d'une personne sont trois conclusions différentes.</p>

</section>

<section id="sources-fournisseur" markdown="1">

## Comparer un rapport commercial sans recopier ses conclusions

Les documents PDF et HTML sont des sources attribuées au fournisseur, pas une
validation indépendante. Lire les tableaux des méthodes et les valeurs source,
pas des compteurs animés sauvegardés ou la longueur des barres. Conserver
assemblage, modèle PGS exact, page, texte original et empreinte du document.
Un modèle public ressemblant à celui du fournisseur n'est pas nécessairement
le même : mêmes loci et noms de maladie ne suffisent pas.

Chaque comparaison distingue concordance littérale, discordance, génotype
absent du VCF, allèles non documentés et orientation non attestée. Une
discordance ne prouve pas lequel des deux appels est biologiquement exact ;
elle impose de comparer les versions, puis les lectures ou une confirmation
indépendante si l'information est importante. Ces comparaisons sont conservées
dans un audit technique séparé, pas dans le rapport personnel ni ses exports
par panel. La personne en est informée à part ; aucune conclusion locale n'est
changée pour faire correspondre les sources. Les réserves nécessaires sur la
fiabilité d'un résultat restent dans le rapport. Les profils étoilés ne sont pas
convertis en diplotypes validés à partir d'un PDF.

Les mentions « Not currently supported » signifient **non évalué**, même si
une carte les étiquette « NOT DETECTED ». Un panel annoncé de plusieurs dizaines
de conditions n'est pas entièrement vérifié si seules quelques cartes sont
lisibles ou si plusieurs pages sont blanches.

Un percentile théorique fondé sur Hardy–Weinberg et indépendance des loci est
explicitement distingué d'une calibration empirique adaptée à la population.
Même une calibration empirique ne constitue pas, seule, une probabilité absolue
de maladie. Un résultat du fournisseur ne figure dans le rapport personnel que
si sa reproduction indépendante et sa qualification ont été vérifiées. Sinon,
il reste dans un audit séparé, même lorsque son percentile dépasse ±1 sigma.
Ce retrait concerne aussi les tableaux et CSV du rapport, ses annexes et ses
pièces jointes. Un calcul partiel n'est pas une reproduction du modèle complet.

La navigation par panels propose un export CSV complet de ses lignes source,
avec protection contre les formules de tableur. Pour un grand modèle, annoncer
le nombre de lignes de l'aperçu, celui des observations dans le CSV et celui
des sites ciblés par le modèle. Ne pas confondre ces trois périmètres.

<p class="caveat">Un catalogue de milliers d'associations GWAS ne teste pas des milliers de maladies causales : les études, phénotypes, loci et observations sont des unités différentes. Les associations gardent leurs publications et ne sont pas additionnées pour fabriquer un score.</p>

</section>

<section id="vocabulaire" data-document-only markdown="1">

## Des mots qui empêchent les mauvaises interprétations

| Terme | Signification dans le rapport | Ce qu'il ne permet pas de conclure seul |
|---|---|---|
| Variant observé | Allèle explicitement appelé dans le VCF | Appel vrai, maladie ou caractère causal |
| Hétérozygote | Deux allèles différents à un site diploïde | Deux copies du gène atteintes |
| Homozygote | Même allèle appelé sur les deux copies au même site | Pathogénicité du variant |
| Deux variants dans un gène | Deux changements distincts ont été observés | Présence sur les deux copies du gène |
| En cis | Les variants sont sur une même copie | Atteinte de la seconde copie |
| En trans | Les variants sont sur des copies différentes | Deux variants causaux pour une même maladie |
| Hétérozygotie composite | Deux variants différents en trans | Diagnostic sans pertinence clinique des deux variants |
| P/LP | Variant pathogène ou probablement pathogène dans une interprétation donnée | Pénétrance complète ou maladie chez tout porteur |
| VUS | Signification clinique incertaine | Indication de traitement ou de PGT |
| Candidat pLoF | Perte de fonction prédite à vérifier | Classification P/LP ou mécanisme approprié à toute maladie du gène |
| Région non évaluable | Données insuffisantes pour un appel fiable | Génotype de référence |

Les recommandations ACMG/AMP précisent qu'un VUS ne doit pas servir de base à une
décision clinique. Un candidat est donc décrit comme un élément à interpréter,
sans le transformer en maladie, portage causal ou indication de FIV. <a data-cite="acmg"></a>

</section>

<section id="fiche" data-document-only markdown="1">

## La fiche obligatoire de chaque variant intéressant

La fiche commence par le résultat **chez cette personne** : allèle appelé,
génotype, zygosité et statut de validation. Elle explique ensuite les maladies
pertinentes pour le **variant exact**, en distinguant les associations générales
du gène. Une page ClinVar agrégeant plusieurs maladies ne prouve pas que son
variant cause toutes ces maladies selon tous les modes de transmission.

| Bloc | Contenu minimal |
|---|---|
| Identité | Assemblage, chromosome, position, REF/ALT, transcript et version, nomenclature HGVS si validée |
| Appel | Génotype, profondeur DP, qualité GQ, équilibre allélique AD, FILTER, métriques complémentaires, fichiers utilisés |
| Effet prédit | Faux-sens, décalage du cadre, codon stop, épissage, indel complexe ; transcript pertinent et limites |
| Preuve | ClinVar par variant/maladie, statut de revue, désaccords, date, fréquence populationnelle et association gène-maladie |
| Maladie | Nom compréhensible, mécanisme, transmission, âge habituel, organes touchés, manifestations et variabilité |
| Chez un porteur | Effets attendus d'un seul allèle, pénétrance lorsqu'elle est connue, absence de conclusion lorsqu'elle ne l'est pas |
| Pour le couple | Variant(s) de chaque parent, maladie commune réellement pertinente, phase et calcul autorisé ou non |
| Suite immédiate | Revue de l'appel, confirmation clinique, complément de couverture, relecture de la classification |
| Diagnostic orienté | Consultation et examens seulement si le génotype validé, les symptômes ou l'histoire familiale le justifient |
| Prise en charge | Traitement établi, mesures de suivi ou absence de traitement spécifique, conditionnés au diagnostic |
| Sources | Liens directs vers GeneReviews, Orphanet, ClinVar, recommandations et publications originales |

Les symptômes sont décrits en phrases concrètes, avec leurs manifestations
légères et sévères. Le rapport doit aussi dire qu'une personne ne présente pas
nécessairement ces symptômes parce qu'un variant incertain a été trouvé.

</section>

<section id="preuves" markdown="1">

## Classer les preuves avant d'afficher les maladies

Une conséquence informatique délétère demande une revue du mécanisme : perte de
fonction effectivement causale pour cette maladie, transcript exprimé dans le
tissu pertinent, exon important, effet réel sur l'épissage et possibilité de
dégradation de l'ARN. Les variantes en fin de gène, les isoformes alternatives,
les indels voisins et les régions paralogues peuvent modifier l'interprétation.
La règle « plusieurs pertes de fonction pathogènes existent ailleurs dans ce
gène » est un filtre de recherche, pas une classification clinique. <a data-cite="pvs1"></a>

Le rapport sépare cinq niveaux : résultat P/LP appelé avec une qualité suffisante
à confirmer cliniquement ; candidat moléculaire plausible ; VUS ; appel complexe
ou techniquement douteux ; résultat non évalué. Les résultats exclus restent
inspectables avec leur motif, mais ne sont pas additionnés aux portages établis.

Une fréquence absente d'une table est notée **non renseignée**. Elle n'est pas
remplacée par zéro. Une fréquence dans une population ne suffit pas à prédire
l'effet chez une personne. Les seuils employés et leurs unités — fréquence
d'allèle, fréquence de portage ou risque par grossesse — sont explicités.

<p class="caveat">Une association gène-maladie peut être dominée par un mécanisme précis : gain de fonction, effet dominant négatif, susceptibilité ou combinaison particulière de variants. Le nom du gène ne détermine pas à lui seul le modèle de risque.</p>

</section>

<section id="couple" markdown="1">

## Écrire une comparaison de couple sans ambiguïté

La synthèse affiche séparément le **nombre de gènes partagés avec des variants
à revoir** et le **nombre de maladies récessives pour lesquelles un portage
causal compatible est identifié chez les deux parents**. Les mutations peuvent
être différentes dans le même gène. Inversement, partager un gène avec des VUS
ne démontre pas un risque récessif établi.

Pour chaque gène commun, une table montre les variants d'un parent puis ceux de
l'autre, leur classement, leur génotype et leur phase. Les observations à
l'intérieur d'une personne ne sont pas confondues avec la transmission de cette
personne à son enfant. Deux génotypes `0/1` ne démontrent ni cis ni trans.

<details class="explanation prose" markdown="1"><summary>Exemple pédagogique : un recouvrement TTN</summary>

Un parent peut avoir un variant faux-sens incertain ; l'autre, un autre faux-sens
et une insertion complexe. Cela donne un gène commun et trois observations au
total, **pas une double atteinte récessive démontrée chez chaque parent**. La
phase des deux observations du second parent peut être inconnue. L'insertion
doit être validée ; son effet dépend aussi du transcript et de l'exon. La
classification des faux-sens reste incertaine. Le rapport explique ces points
avant toute discussion de risque, et n'affiche pas « 25 % de maladie ».

</details>

Une probabilité de transmission biallélique est calculée seulement lorsque les
allèles pertinents et les haplotypes transmissibles sont établis. Pour deux
porteurs hétérozygotes d'une même maladie récessive, le modèle classique donne
une transmission biallélique de 25 % par grossesse ; l'expression clinique
demande encore d'examiner pénétrance, sévérité et compatibilité de la combinaison
d'allèles. Les VUS ne sont pas injectés dans ce calcul.

Les risques liés à l'X ou dominants peuvent exister sans variante commune chez
l'autre parent et ont leur propre section. Les chiffres populationnels obtenus
avant le test du partenaire restent des scénarios historiques : ils ne sont pas
le risque mesuré du couple après ses deux tests.

</section>

<section id="diagnostic" markdown="1">

## Donner des suites concrètes sans prescrire à partir d'une incertitude

Chaque résultat suit un parcours lisible : vérifier l'appel et la classification ;
le faire confirmer avec une méthode adaptée ; confronter le résultat aux symptômes
et à l'histoire familiale ; discuter, si pertinent, des examens orientés ; décrire
la prise en charge d'une maladie **si elle est diagnostiquée**.

La confirmation peut nécessiter un examen du BAM, un séquençage indépendant,
une méthode distinguant un pseudogène, une analyse de délétions/duplications ou
une étude familiale de la phase. Le BAM seul ne résout pas nécessairement les
régions très homologues. Une étude ARN peut être envisagée par le laboratoire
pour un mécanisme d'épissage, lorsqu'un tissu et une méthode pertinents existent.

Une fiche de traitement distingue recommandations établies, prise en charge
symptomatique, essais cliniques et hypothèses. Elle ne transforme pas une voie
métabolique en conseil alimentaire ni un VUS en indication de médicament.
Les posologies, arrêts de traitements et décisions reproductives nécessitent
une interprétation clinique propre à la personne.

</section>

<section id="pharmacogenomique" markdown="1">

## Pharmacogénomique : relier un résultat valide à un médicament pertinent

Le rapport indique les marqueurs recherchés, les allèles explicitement observés,
la couverture du test, les médicaments concernés et la source des recommandations.
Il distingue marqueur isolé, haplotype, diplotype et phénotype métaboliseur. Une
conclusion de métaboliseur « normal » n'est pas produite si des allèles requis
sont absents d'un VCF ne contenant que les variants.

Les recommandations CPIC sont des aides à l'interprétation d'un résultat
pharmacogénétique disponible ; elles ne constituent pas un choix autonome de
traitement. Une fiche précise ce qu'il faudrait confirmer avant une prescription,
et rappelle les autres facteurs pertinents : médicament actuel, indication,
co-prescriptions, fonction rénale/hépatique et contexte clinique. <a data-cite="cpic"></a>

Les gènes complexes tels que CYP2D6, les allèles HLA et certaines variations
structurales exigent des méthodes adaptées. Un proxy de haplotype est identifié
comme tel. Les variants mitochondriaux sont vérifiés sur la référence adéquate
avec leur hétéroplasmie lorsque les données la permettent. Les informations sont
priorisées si un médicament concerné est réellement utilisé ou envisagé.

</section>

<section id="nutrition" markdown="1">

## Nutrition : distinguer un effet interprétable d'un conseil commercial

Les résultats peuvent documenter une composante génétique de la persistance de
la lactase, certains troubles du métabolisme du fer ou de l'alcool. Chaque fiche
doit préciser ce que le marqueur indique dans les populations où il a été validé
et ce qu'il ne prédit pas : symptômes digestifs, carence, intoxication, régime
idéal ou dose de supplément.

La tolérance au lactose se confronte aux symptômes et, si nécessaire, à un examen
adapté. Un variant HFE se confronte au génotype complet et aux marqueurs du fer,
pas à une prescription automatique de saignées ou de restrictions alimentaires.
Un risque lié à l'alcool ne signifie pas que la consommation est médicalement
recommandée chez d'autres génotypes.

Les scores de régimes, performances sportives et associations faibles ne sont
pas mélangés aux résultats validés. Une information éducative peut être incluse
avec sa limite ; elle ne devient pas une recommandation personnalisée en
l'absence de bénéfice démontré et de contexte clinique.

La première lecture n'affiche pas de fiche « rien détecté » ni de VUS comme
information de santé. Elle présente seulement les génotypes pertinents
explicitement observés, avec « votre résultat / ce que cela signifie / ce que
vous pouvez faire ». Les associations de faible effet, comme certains marqueurs
de glycémie, sont étiquetées comme telles et ne justifient pas de supplément.

</section>

<section id="qualite" data-document-only markdown="1">

## Qualité et concordance : des contrôles indépendants, pas des garanties

Les premiers contrôles portent sur l'identité des fichiers, l'échantillon unique,
la référence génomique, les sommes de contrôle, le nombre de variants, DP/GQ,
l'équilibre allélique, les filtres et la disponibilité de la phase. Les
profondeurs calculées aux sites variantés ne sont pas une couverture de tout le
génome. Une référence absente du VCF ne prouve pas qu'elle a été correctement
séquencée. <a data-cite="gatk"></a>

| Contrôle complémentaire | Utilité | Conditions et limites |
|---|---|---|
| ABO | Concordance avec un groupe sanguin connu | Haplotype, allèles rares et FUT1/FUT2 ; confirmer par sérologie |
| RhD | Concordance immunohématologique | RHD comporte des délétions et allèles hybrides ; VCF SNP seul souvent insuffisant |
| Chromosomes sexuels | Cohérence du contenu chromosomique et de la ploïdie | Pas une détermination de genre ; karyotype/mosaïcisme demandent des données adaptées |
| Pigmentation | Information de concordance indicative | Caractère polygénique ; un SNP ne vérifie pas l'identité |
| Empreinte génétique | Échantillons dupliqués ou inversés | Panel validé, taux d'appels et métriques calibrées |
| Parenté génétique | Vérifier relations annoncées | IBD/kinship sur marqueurs indépendants ; la simple intersection de variants ne suffit pas |
| Contamination | Détecter mélange d'échantillons | Méthodes dédiées, souvent lectures ou fréquences de génotypes calibrées |

Une discordance est investiguée ; elle n'est ni cachée ni attribuée d'emblée à
une maladie. Les réponses connues — groupe sanguin, histoire familiale,
traitements — sont datées et leur source déclarée. Elles restent facultatives.

</section>

<section id="ascendance" markdown="1">

## Ascendance génétique et références populationnelles

L'ethnicité est une identité et un contexte social ; elle ne se déduit pas d'une
poignée de SNP. Une estimation d'ascendance génétique demande un panel de
référence documenté, des marqueurs indépendants, une méthode calibrée et une
incertitude explicite. Elle n'est pas une attribution de nationalité. Sans ce
panel et ces contrôles, le résultat est **non calculé**, et non « français » ou
« européen » présumé. <a data-cite="ancestry"></a>

L'utilisation de gnomAD NFE dans une estimation de risque est un choix de
référence, pas une estimation d'ascendance de la personne. Les fréquences
globales, les fréquences par population et les bornes entre populations ont
des dénominateurs différents. Une comparaison française doit distinguer une
cohorte française spécifique d'un proxy européen non finlandais.

La méthode sans BAM compare des marqueurs autosomiques explicites à des
références publiques GRCh37. Le biais d'un VCF de variants est reproduit sur
les références d'entraînement et de validation : un site absent n'est jamais
transformé en 0/0. La PCA, les contrôles indépendants, les sous-panels et les
contrôles par chromosomes sont enregistrés. Le résultat décrit une affinité,
pas une décomposition en pourcentages ancestraux. Le graphique doit utiliser
exactement le modèle PCA ayant produit le résultat, vérifié par empreinte.

</section>

<section id="scores-polygeniques" markdown="1">

## Scores polygéniques : classement qualifié, pas probabilité de maladie

Les poids et métadonnées proviennent de PGS Catalog ; les distributions de
référence sont identifiées par modèle, variantes réellement utilisées et
population. Les allèles, le brin, la référence génomique et la ploïdie sont
harmonisés. Les génotypes inconnus produisent des bornes, pas des dosages nuls.
<a data-cite="pgs-catalog"></a> <a data-cite="1000g"></a>

Le filtre éditorial est **strictement au-delà de ±1 sigma** : `abs(z) > 1`,
ou un intervalle entièrement au-dessus de +1 ou au-dessous de −1. Il ne suffit
pas à rendre un score valide. Un sous-score dont les sommes ont été vérifiées
contre un logiciel indépendant n'est pas automatiquement un substitut validé
du PGS publié. Ces résultats de recherche restent dans un audit privé séparé,
hors du rapport personnel, tant que leur reproduction et leur qualification
prédictive ne sont pas validées.

Un percentile signifie « position du score dans la population de référence ».
Ce n'est ni une probabilité de développer une maladie ni un risque à vie.
L'âge, le contexte clinique et les autres facteurs exigent un modèle de risque
absolu distinct, calibré et validé. Les captures d'un autre rapport ne sont
jamais une source de poids, de génotypes ou de résultats personnels.

</section>

<section id="html" markdown="1">

## Présentation HTML et accès à l'information

L'HTML doit fonctionner hors ligne, sans envoyer de génotypes à un serveur.
La synthèse et les limites importantes restent visibles. Des filtres par
catégorie et une recherche par gène aident à parcourir les fiches ; masquer une
fiche par un filtre ne change pas les calculs. Les longues séquences d'indels
doivent se replier sans casser la mise en page.

La vue document contient les tableaux exhaustifs, la provenance et les sources.
Une vue de présentation peut résumer les mêmes conclusions sans les modifier.
La page doit rester utilisable au clavier, sur mobile et à l'impression.
Les codes couleur ont aussi des libellés textuels. Les liens médicaux sont
spécifiques au contenu et distinguent source originale et ressource grand public.

</section>

<section id="reproductibilite" data-document-only markdown="1">

## Reprendre et préserver une analyse

Le [runbook de reproduction](REPRODUCTION.md) relie les fonctions de calcul,
les annotations de recherche, le rendu et la sauvegarde. Un agent commence par
les entrées figées, leurs empreintes et les références/version/licences utilisées,
pas par des résultats recopiés d'une autre version. Un BAM est un complément
possible pour valider les lectures, pas un prérequis du chemin VCF actuel.

Les profils réutilisable et personnel privé non commercial restent séparés.
Une nouvelle annotation prédictive ne change pas les filtres, la phase ou le
diagnostic. Les résultats non qualifiés et comparaisons fournisseur/VCF restent
dans leurs audits privés, pas dans les panels ou pièces jointes du rapport.

Chaque export annonce son périmètre : allèles uniques, associations ou lignes
score/transcrit/observation ne sont pas des compteurs interchangeables. Avant
une annonce publique, les nombres de conditions, gènes et variants doivent être
recalculés depuis le snapshot exact et les filtres du code exécuté. Distinguer
catalogue inventorié et appels réellement évaluables : un site absent du VCF
reste inconnu, et le catalogue ne constitue pas un test validé pour chaque
maladie. Les futurs runs nécessitent un manifeste complémentaire des références
et du code ; les JSON actuels du moteur ne figent pas toutes ces empreintes.
Ne pas reconstruire rétrospectivement une preuve d'exécution manquante.

Un rapport peut se lire seul tout en ayant des CSV annexes. La fonction
`prepare_report_delivery()` conditionne chaque HTML choisi et ses références
statiques dans un ZIP vérifié, sans ajouter les audits voisins ni envoyer de
message. Extraire les rapports dans des dossiers distincts ; autoriser
explicitement les rapports et leur destinataire avant partage.

Une sauvegarde privée vérifiée préserve les sources locales ; elle ne démontre
pas que l'analyse ou l'interface a été reproduite sur un autre hôte. La version
du guide reste distincte de celle des rapports patients : réviser ce document
ne recalcule ni ne réétiquette ces rapports.

</section>

<section id="validation" data-document-only markdown="1">

## Vérifications avant de remettre le rapport

La génération utilise les résultats structurés et leurs fichiers source. Les
chiffres affichés sont recalculés : variants, gènes distincts et individus ont
leurs propres compteurs. Une table de contrôle relie chaque chiffre à son
expression et à la révision d'analyse. Les variantes du même site multi-gène
ne sont pas comptées plusieurs fois comme observations indépendantes.

La relecture vérifie les affirmations cliniques, les coordonnées, les transcripts,
les références, la conservation de l'incertitude, les compteurs et les scénarios
de transmission. Les liens internes, filtres, détails dépliables, modes de
lecture, navigation clavier, mobile et absence d'appels réseau sont testés.

Un historique enregistre les reclassifications, nouvelles données et corrections
d'interprétation. Le guide et les fonctions peuvent être partagés ; les génomes
et rapports nominaux suivent une gestion privée distincte. Toute version du
rapport indique la date des données, la date de l'interprétation et les éléments
non analysés.

</section>
