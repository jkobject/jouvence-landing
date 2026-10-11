# Objectif éditorial

Public : personnes ayant réalisé un séquençage et couples préparant une grossesse,
avec une seconde lecture possible par un professionnel de santé.

Question : que montrent réellement les données, quelles maladies sont pertinentes
pour chaque résultat et que faut-il vérifier avant d'agir ?

Message central : un rapport utile relie chaque variant à ses preuves, au mode de
transmission et à une suite concrète. Chaque association utile explique le résultat,
la physiologie, l'ampleur publiée, l'action possible et la limite personnelle ;
le rapport ne confond pas observation, prédiction, portage et diagnostic.

Les rapports personnels et les données du couple restent privés. Ce guide et
les fonctions de génération sont génériques. L'HTML est autonome, en français,
lisible hors ligne, avec détails dépliables et filtres qui ne modifient pas les
conclusions. La présentation privilégie la synthèse ; les tableaux exhaustifs
restent dans la vue document.

Un changement de format réutilise les analyses figées du run concerné. Les nouvelles
extractions pharmacogénomiques et de traits sont identifiées séparément. Les
corrections d'interprétation clinique sont tracées, sans modifier les VCF source.

## Révision de lecture inspirée des exemples fournis

La première vue ne montre que les résultats utiles : observations P/LP à
confirmer, candidats pLoF stricts explicitement séparés et marqueurs ciblés
avec génotype et qualité exploitables. Chaque carte répond à trois questions :
votre résultat, sa signification et la prochaine étape. La signification précise
le phénotype réellement étudié et l'ampleur publiée : OR/HR/β, intervalle de
confiance disponible, unité, population et comparateur. L'action peut être
« aucune action spécifique justifiée » ; un chiffre non documenté reste absent.
Les VUS, appels rejetés,
sites absents et contrôles exhaustifs restent dans une annexe fermée par défaut.
Les limites qui changent l'interprétation restent visibles sur la carte.

L'ascendance est calculée sur un panel public figé, avec validation indépendante
et contrôles de sensibilité reproduisant les données variants-only. La PCA
illustre une proximité de références, pas des proportions ancestrales ou une
nationalité. Les images d'exemple ne sont jamais une source de génotypes.

Les scores polygéniques requièrent poids, harmonisation et référence de
calibration documentés. Le guide explique la somme poids × dosage, les moments
théoriques sous Hardy–Weinberg, les covariances LD omises sous indépendance et
la comparaison empirique selon l'ascendance. Normalisation et calibration d'un
risque absolu restent deux étapes distinctes. Seuls des résultats qualifiés au-delà de ±1 sigma,
y compris une borne intégralement hors de cet intervalle, peuvent entrer dans
la synthèse. Un percentile classe un score ; il ne mesure pas une probabilité
de maladie. Les scores incomplets ou non qualifiés restent hors de cette
synthèse clinique. Un volet exploratoire distinct peut être activé explicitement,
avec un sous-panel observé, un modèle HWE conditionnel, un z-score et un rang
discret Monte-Carlo calculés localement. Il ne copie aucun résultat fournisseur,
n'impute aucun inconnu et ne revendique pas le PGS complet. La première lecture
montre `abs(z)>1`, avec les autres calculs repliés et CSV par modèle/référence ;
une méthode partagée évite de répéter les avertissements. Ce seuil éditorial n'est pas diagnostique. La validation arithmétique d'un sous-score ne démontre pas
sa performance prédictive clinique.

## Panels et audit séparé des rapports du fournisseur

Les rapports PDF/HTML sont importés comme données non exécutées, avec provenance
par document et par page, puis confrontés aux génotypes explicitement présents
dans les VCF. Imports et comparaisons restent dans un audit technique séparé.
Le rapport personnel ne réaffiche aucun résultat du fournisseur non reproduit
indépendamment. Le retrait porte aussi sur les tableaux, CSV et pièces jointes
de l'HTML, pas seulement sur les badges. Les réserves cliniques et techniques
qui changent la fiabilité d'une observation réellement extraite des données
locales restent visibles ; une observation n'est pas un diagnostic.
L'audit distingue le périmètre annoncé par un fournisseur du contenu réellement
extractible : une page blanche ou un test non pris en charge n'est pas un négatif.

L'enrichissement GWAS sépare associations, rsID, publications, mesures étudiées
et appels explicites. Il ne transforme pas des milliers d'associations en
milliers de maladies détectées. Un pQTL concerne une mesure de protéine ; le
chronotype une préférence d'horaires de sommeil. Le nom du gène ne permet pas
de passer de l'un à l'autre. Les plus petites valeurs P ne définissent pas
l'utilité clinique ; les lignes corrélées ne s'additionnent pas pour créer un
risque ou un percentile. Chaque effet cité garde son passage primaire justificatif
et le niveau de lecture réel, notamment méthodes/suppléments pour les unités.
Seuls les scores dont la reproduction et la qualification ont été vérifiées
peuvent figurer dans le volet qualifié. Le volet exploratoire séparé présente
uniquement les indices locaux conditionnels explicitement calculés, pas les
valeurs importées ni les anciens audits bruts. Les blocages restent dans
l'audit séparé. Les gros modèles disposent d'une table limitée et
d'un CSV complet des observations, sans alourdir l'HTML par des copies du
catalogue entier.
