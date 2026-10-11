# Consignes pour les agents

## Manière de travailler

- Préférer des fonctions Python simples à des scripts ponctuels. Documenter les entrées, sorties, limites et blocs importants ; éviter les couches et sous-fonctions inutiles.
- Utiliser `uv`, pas `pip` ni un environnement `venv` créé manuellement.
- Si l'objectif est ambigu, demander avant d'implémenter. Si une documentation ou une référence est inaccessible, expliquer le blocage et demander les informations manquantes ; ne pas inventer les données.
- Préserver les modifications existantes. Si l'extension devient difficile à comprendre, proposer un refactor ciblé plutôt qu'ajouter des couches.
- Modifier les fichiers avec `apply_patch`, vérifier les changements et leurs effets. Pour une tâche demandant une livraison de code ou de documentation, le travail n'est pas terminé avant le commit et le push des seuls fichiers concernés ; signaler tout blocage, ne jamais forcer le push.

## Avant un travail génomique

Lire dans cet ordre :

1. [Runbook de reprise](docs/rapport_genomique/REPRODUCTION.md).
2. [API et périmètres actuels](tools/jouvence_genomics/README.md).
3. [Guide de qualité des rapports](docs/rapport_genomique/report.md).
4. [Conditions des ressources](tools/jouvence_genomics/resource_use_review.json), pour tout ajout de référence.

Vérifier ensuite les signatures et tests des fonctions concernées. Un document décrit une méthode, pas la preuve qu'un run particulier l'a exécutée. Ne pas relancer une analyse coûteuse simplement pour changer sa présentation.

## Données et conclusions

- Garder génomes, génotypes, rapports nominatifs, références téléchargées, logs et manifestes privés hors Git. Ne jamais utiliser `git add .` pour une livraison génomique.
- Ne transmettre aucune donnée génétique ni position personnelle à un service distant sans autorisation explicite. Une sauvegarde privée ou un partage autorisé a son propre périmètre ; ce n'est pas une autorisation de publication ou de télémétrie.
- Un site absent d'un VCF de variants reste inconnu, jamais `0/0`. Vérifier assemblage, REF/ALT, ploïdie, filtres, phase et empreintes ; ne pas modifier les seuils cliniques pour augmenter les résultats.
- Prédiction moléculaire, classification du variant et diagnostic sont distincts. Les VUS, associations GWAS et annotations AlphaMissense/REVEL/CADD ne deviennent pas des maladies, indications de PGT ou prescriptions.
- Le profil `personal_noncommercial` est explicite, privé et séparé du profil réutilisable ; il ne donne aucun droit commercial et ne valide pas des termes inconnus.
- Respecter les exclusions du rapport définies dans le runbook. Les audits de recherche complets restent séparés ; une limite éditoriale comme ±1 sigma ne valide pas un score.
- Sauvegarder uniquement les sources autorisées vers la destination privée vérifiée. Conserver les originaux sur le Mac, vérifier les reçus et l'intégrité, ne jamais annoncer un transfert vérifié à partir d'un processus seulement lancé.
