# Objectif éditorial

Public : personnes ayant réalisé un séquençage et couples préparant une grossesse,
avec une seconde lecture possible par un professionnel de santé.

Question : que montrent réellement les données, quelles maladies sont pertinentes
pour chaque résultat et que faut-il vérifier avant d'agir ?

Message central : un rapport utile relie chaque variant à ses preuves, au mode de
transmission et à une suite concrète ; il ne confond pas observation, prédiction,
portage et diagnostic.

Les rapports personnels et les données du couple restent privés. Ce guide et
les fonctions de génération sont génériques. L'HTML est autonome, en français,
lisible hors ligne, avec détails dépliables et filtres qui ne modifient pas les
conclusions. La présentation privilégie la synthèse ; les tableaux exhaustifs
restent dans la vue document.

Le changement de format réutilise l'analyse figée du 9 octobre 2026. Les nouvelles
extractions pharmacogénomiques et de traits sont identifiées séparément. Les
corrections d'interprétation clinique sont tracées, sans modifier les VCF source.

## Révision de lecture inspirée des exemples fournis

La première vue ne montre que les résultats utiles : observations P/LP à
confirmer, candidats pLoF stricts explicitement séparés et marqueurs ciblés
avec génotype et qualité exploitables. Chaque carte répond à trois questions :
votre résultat, sa signification et la prochaine étape. Les VUS, appels rejetés,
sites absents et contrôles exhaustifs restent dans une annexe fermée par défaut.
Les limites qui changent l'interprétation restent visibles sur la carte.

L'ascendance est calculée sur un panel public figé, avec validation indépendante
et contrôles de sensibilité reproduisant les données variants-only. La PCA
illustre une proximité de références, pas des proportions ancestrales ou une
nationalité. Les images d'exemple ne sont jamais une source de génotypes.

Les scores polygéniques requièrent poids, harmonisation et référence de
calibration documentés. Seuls des résultats qualifiés au-delà de ±1 sigma,
y compris une borne intégralement hors de cet intervalle, peuvent entrer dans
la synthèse. Un percentile classe un score ; il ne mesure pas une probabilité
de maladie. Les scores incomplets ou non qualifiés restent dans l'audit, sans
classement inventé. La validation arithmétique d'un sous-score ne démontre pas
sa performance prédictive clinique.

## Panels et audit séparé des rapports du fournisseur

Les rapports PDF/HTML sont importés comme données non exécutées. Leurs résultats
restent attribués au fournisseur et confrontés aux génotypes explicitement
présents dans les VCF, avec provenance par document et par page. Cette comparaison
reste dans un audit technique séparé, communiqué à part. Le rapport personnel
traduit les explications et offre un export CSV par panel, sans tableau de
discordances ni mélange des génotypes de deux sources. Les réserves cliniques
et techniques qui changent la fiabilité d'une observation locale restent visibles.
Elle distingue le périmètre annoncé par un fournisseur du contenu réellement
extractible : une page blanche ou un test non pris en charge n'est pas un négatif.

L'enrichissement GWAS sépare associations, rsID, publications, mesures étudiées
et appels explicites. Il ne transforme pas des milliers d'associations en
milliers de maladies détectées. Les scores externes hors ±1 sigma peuvent être
inspectés dans une rubrique attribuée et non reproduite ; ils ne deviennent pas
des scores locaux qualifiés. Les gros modèles disposent d'une table limitée et
d'un CSV complet des observations, sans alourdir l'HTML par des copies du
catalogue entier.
