"""Rapports génomiques français, autonomes et traçables.

Les fonctions réutilisent les JSON d'analyse figés. Elles enrichissent la lecture
avec des fiches cliniques et un panel ciblé de traits, sans transformer un VUS en
diagnostic. Le rendu utilise la fonction ``render_report`` du guide HTML local.
Les données personnelles et les HTML nominaux doivent rester hors du dépôt.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import shutil
import sys
from base64 import b64encode
from collections import defaultdict
from datetime import datetime, timezone
from html import escape
from importlib.metadata import version
from pathlib import Path
from urllib.parse import urlsplit

from .traits import audit_traits
from .catalogs import review_stars
from .vcf import normalize_allele
from .patient_view import (select_patient_signals, render_patient_sections,
                          render_population_results, ancestry_plot_spec)
from .report_panels import build_report_panels
from .panel_view import render_panel_dashboard


HERE = Path(__file__).parent
TIERS = {
    "plp": "P/LP : résultat à confirmer cliniquement",
    "plof": "Candidat de perte de fonction à confirmer",
    "vus": "Signification clinique incertaine",
    "mechanism": "Effet et mécanisme de maladie à préciser",
    "complex": "Appel complexe : haplotype et lectures à revoir",
    "technical": "Qualité insuffisante : résultat non retenu",
    "unassigned": "Association ou attribution non établie",
}
CONSEQUENCES = {
    "missense_variant": "changement d'acide aminé (faux-sens)",
    "frameshift_variant": "décalage du cadre de lecture prédit",
    "splice_acceptor_variant": "atteinte prédite du site accepteur d'épissage",
    "splice_donor_variant": "atteinte prédite du site donneur d'épissage",
    "inframe_insertion": "insertion conservant le cadre de lecture",
    "inframe_deletion": "délétion conservant le cadre de lecture",
    "inframe_delins": "remplacement conservant le cadre de lecture",
    "stop_gained": "codon d'arrêt prématuré prédit",
    "start_lost": "perte prédite du codon de départ",
    "stop_lost": "perte prédite du codon d'arrêt",
}
GENERAL_SOURCES = {
    "acmg": {"label": "ACMG/AMP : interprétation des variants, 2015", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4544753/"},
    "pvs1": {"label": "ClinGen : interprétation des pertes de fonction PVS1", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6185798/"},
    "cpic": {"label": "CPIC : recommandations pharmacogénétiques", "url": "https://cpicpgx.org/guidelines/"},
    "gatk": {"label": "Broad GATK : qualité et filtrage des variants germinaux", "url": "https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants"},
    "ancestry": {"label": "NHGRI : descripteurs populationnels, ascendance et ethnicité", "url": "https://www.genome.gov/about-genomics/policy-issues/population-descriptors-in-genomics"},
    "pgs-catalog": {"label": "PGS Catalog : poids et métadonnées publics des scores", "url": "https://www.pgscatalog.org/downloads/"},
    "1000g": {"label": "1000 Genomes : populations et variations de référence", "url": "https://www.internationalgenome.org/data-portal/population"},
    "clinvar-review": {"label": "ClinVar : niveaux de revue et critères des assertions", "url": "https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/"},
}


def collect_variant_records(result: dict) -> list[dict]:
    """Rassembler les observations uniques, avec un statut de lecture explicite.

    Les associations de gène ne sont pas des diagnostics. RNF213 reçoit un
    libellé spécifique au mécanisme et le marqueur sanguin XG sort du décompte
    clinique. Les annotations multi-gènes PRSS1/TRB restent inspectables.
    """

    records = []
    for item in result["all_reviewed_clinvar_matches"]:
        gene = item["clinvar"]["gene"]
        if not gene:
            continue  # Le phénotype XG est un trait sanguin, pas une maladie retenue.
        tier = "technical" if item["quality_flags"] else "plp"
        if tier == "plp" and item["clinvar"].get("review_status") and review_stars(item["clinvar"]["review_status"]) == 0:
            tier = "unassigned"  # Une ancienne étoile mal décodée ne promeut pas une assertion sans critères.
        if gene == "TRB":
            tier = "unassigned"
        records.append({**item, "gene": gene, "report_tier": tier})
    novel = result["novel_variant_analysis"]
    for field, tier in (("strict_plof_candidates", "plof"),
                        ("complex_plof_candidates", "complex"),
                        ("protein_altering_vus", "vus"),
                        ("technical_review_candidates", "technical")):
        for item in novel.get(field, []):
            report_tier = "mechanism" if item["gene"] == "RNF213" and tier != "technical" else tier
            records.append({**item, "report_tier": report_tier})
    unique = {}
    for item in records:
        key = (*normalize_allele(item["chrom"], item["pos"], item["ref"], item["alt"]), item["gene"])
        unique.setdefault(key, item)
    return sorted(unique.values(), key=lambda x: (x["report_tier"] != "plp", x["gene"], x["pos"]))


def render_variant_catalog(records: list[dict], catalog: dict, references: dict) -> str:
    """Créer une fiche par gène avec toutes ses mutations et son contexte clinique.

    Chaque mutation décrit son appel et sa suite de validation. Les descriptions
    de maladie, symptômes et prise en charge sont partagées dans la même fiche
    pour éviter de répéter des dizaines de fois une maladie à appels complexes.
    """

    groups = defaultdict(list)
    for item in records:
        groups[item["gene"]].append(item)
    missing = set(groups) - set(catalog) - {"TRB"}
    if missing:
        raise ValueError(f"Fiches cliniques manquantes : {sorted(missing)}")
    html = ['<div data-filter-scope="variants" class="prose">', _filter_controls()]
    for gene, items in groups.items():
        tiers = sorted({item["report_tier"] for item in items})
        card = catalog.get(gene, {
            "title": "Attribution multi-gène du locus PRSS1/TRB à revoir",
            "interpretation": "Cette ligne provient du même enregistrement ClinVar que PRSS1. Elle ne démontre pas une maladie TRB et n'est pas une seconde mutation indépendante.",
            "sources": [{"label": "ClinVar : enregistrement du locus", "url": "https://www.ncbi.nlm.nih.gov/clinvar/variation/11877/"}],
        })
        title = card.get("title") or "; ".join(x["name"] for x in card.get("conditions", []))
        opened = " open" if "plp" in tiers else ""
        html.append(f'<details id="gene-{escape(gene)}" class="gene-card" data-gene-card '
                    f'data-tiers="{escape(" ".join(tiers))}" data-search="{escape(gene + " " + title, quote=True)}"{opened}>'
                    f'<summary><strong>{escape(gene)}</strong> · {len(items)} observation(s) '
                    f'<span class="status">{escape(TIERS[tiers[0]])}</span></summary>')
        html.append(f'<p class="inline-note">{escape(title)}</p>')
        for item in items:
            tier = item["report_tier"]
            ref, alt = item["ref"], item["alt"]
            long = max(len(ref), len(alt)) > 35
            label = f"chr{str(item['chrom']).removeprefix('chr')}:{item['pos']}"
            if long:
                label += f" · {len(ref)} base(s) de référence / {len(alt)} base(s) alternatives"
            else:
                label += f" {ref}>{alt}"
            clinical = item.get("clinvar", {})
            classification = clinical.get("clinical_significance") or ", ".join(item.get("clinvar_significances", []))
            phase = (f"bloc {item['phase_block']}, haplotype(s) {item['alt_haplotypes']}"
                     if item.get("phased") else "non résolue dans les données")
            consequence = CONSEQUENCES.get(item.get("consequence"), "correspondance exacte au catalogue ClinVar")
            html.append(f'<div class="variant-row"><strong><code>{escape(label)}</code></strong> '
                        f'<span class="status status-{tier}">{escape(TIERS[tier])}</span>')
            if long:
                html.append(f'<details class="sequence-details"><summary>Voir la séquence VCF complète</summary>'
                            f'<code>{escape(ref)}&gt;{escape(alt)}</code></details>')
            balance = item.get("allele_balance")
            balance_text = "non renseigné" if balance is None else f"{100 * balance:.1f} %"
            html.append('<div class="fact-grid">'
                        f'<span>Source / personne : {escape(str(item.get("source_file", "non renseignée")))}</span>'
                        f'<span>Génotype appelé : <strong>{escape(str(item.get("genotype")))}</strong></span>'
                        f'<span>Phase : {escape(phase)}</span>'
                        f'<span>DP : {item.get("DP") if item.get("DP") is not None else "non renseignée"} · '
                        f'GQ : {item.get("GQ") if item.get("GQ") is not None else "non renseignée"}</span>'
                        f'<span>Fraction des lectures ALT : {balance_text}</span>'
                        f'<span>FILTER : {escape(str(item.get("FILTER")))}</span>'
                        f'<span>Transcript : {escape(str(item.get("transcript") or "à documenter"))}</span></div>')
            html.append(f'<p><strong>Effet annoté :</strong> {escape(consequence)} '
                        f'{escape(str(item.get("protein_change") or ""))}.</p>')
            if classification:
                html.append(f'<p><strong>Classement :</strong> {escape(_french_classification(classification))}.</p>')
            stars = clinical.get("stars", item.get("clinvar_review_stars"))
            frozen_stars = stars
            if clinical.get("review_status"):
                stars = review_stars(clinical["review_status"])
            if stars is not None:
                html.append(f'<p><strong>Revue ClinVar :</strong> {stars} étoile(s) dans la référence figée ; '
                            'ce niveau décrit la revue du catalogue, pas la validation de cet échantillon.</p>')
                if stars != frozen_stars:
                    html.append(f'<p class="inline-note">Correction de lecture du statut « {escape(clinical["review_status"])} » : '
                                f'{stars} étoile(s), contre {frozen_stars} dans l’ancien décodage figé. '
                                'Le fichier source est conservé inchangé. <a data-cite="clinvar-review"></a></p>')
            if "gnomad_global_af" in item or "gnomad_nfe_af" in item:
                frequencies = []
                for key, label in (("gnomad_global_af", "globale"), ("gnomad_nfe_af", "NFE")):
                    value = item.get(key)
                    frequencies.append(f'{label} : {value:.6g}' if value is not None else f'{label} : non renseignée')
                html.append('<p><strong>Fréquence allélique gnomAD :</strong> ' + "; ".join(frequencies) +
                            '. Une fréquence non renseignée ne prouve ni absence dans la population ni rareté.</p>')
            if clinical.get("variation_id") or item.get("clinvar_variation_ids"):
                ids = [clinical["variation_id"]] if clinical.get("variation_id") else item["clinvar_variation_ids"]
                html.append("<p>" + " · ".join(f'<a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/{escape(str(x))}/">ClinVar {escape(str(x))}</a>' for x in ids) + "</p>")
            notes = item.get("quality_flags", []) + item.get("complex_region_flags", [])
            if notes:
                html.append(f'<p class="caveat">{escape(_french_flags(notes))}</p>')
            next_step = ("Faire revoir et confirmer l'appel avant toute décision clinique ou reproductive."
                         if tier == "plp" else
                         "Vérifier les lectures, le transcript et l'interprétation du variant ; aucune maladie ni indication de traitement n'est déduite de ce seul signal.")
            html.append(f'<p><strong>Suite pour cette observation :</strong> {next_step}</p>')
            if gene == "RNF213" and item.get("consequence") in {
                "splice_acceptor_variant", "splice_donor_variant", "frameshift_variant", "stop_gained"
            }:
                html.append('<p class="caveat">La prédiction d’épissage/LoF est plausible, mais ni le transcript réellement produit ni la dégradation de l’ARN ne sont démontrés. Ce signal n’est pas un portage récessif établi.</p>')
            if gene in {"TUBB8", "PRSS1"} and tier == "technical" and clinical:
                html.append('<p class="caveat">Pourquoi ce résultat peut différer d’un rapport commercial : '
                            'le variant est bien trouvé dans notre catalogue, mais l’appel de cet échantillon '
                            'n’est pas retenu pour qualité insuffisante. Le classement ClinVar ne valide pas les lectures. '
                            'Les métriques ci-dessus et les lectures d’un BAM, ou une confirmation indépendante, '
                            'permettront de trancher.</p>')
            html.append("</div>")
        html.append('<div class="clinical-context">')
        if card.get("evidence"):
            html.append(f'<p><strong>Solidité et portée de l’association génique :</strong> {escape(card["evidence"])}</p>')
        if card.get("conditions"):
            html.append('<h3>Maladies et transmission</h3><ul>' + "".join(
                f'<li><strong>{escape(c["name"])}</strong> — {escape(c["inheritance"])}</li>' for c in card["conditions"]) + "</ul>")
        elif card.get("inheritance"):
            html.append(f'<p><strong>Transmission :</strong> {escape(card["inheritance"])}</p>')
        for heading, value in (("Manifestations et variabilité", card.get("phenotype")),
                               ("Ce que cela signifie ici", card.get("interpretation") or card.get("caution")),
                               ("Étapes diagnostiques possibles", card.get("diagnostic_steps") or card.get("diagnosis")),
                               ("Prise en charge si une maladie est diagnostiquée", card.get("management")),
                               ("Précisions propres au variant", card.get("variant_specific_notes"))):
            if heading == "Précisions propres au variant" and gene == "MVK" and not any(
                item.get("clinvar", {}).get("variation_id") == "97608" for item in items
            ):
                continue
            if value:
                html.append(f"<h3>{heading}</h3>" + _paragraphs(value))
        html.append('<p class="caveat">La description de cette maladie est une information de contexte. Elle ne signifie pas que la personne présente ce phénotype, notamment pour un VUS ou un appel complexe.</p>')
        html.append(_source_links(card.get("sources", []), references))
        html.append("</div></details>")
    html.append("</div>")
    return "\n".join(html)


def render_trait_sections(traits: dict, references: dict) -> str:
    """Décrire marqueurs observés et sites inconnus dans leurs rubriques cliniques."""

    sections = []
    for category, title, introductory in (
        ("pharmacogenomique", "Pharmacogénomique", "Ces observations ciblées peuvent orienter une confirmation avant un médicament concerné. Elles ne donnent ni un diplotype complet, ni une dose."),
        ("nutrition", "Nutrition et métabolisme", "Les marqueurs se confrontent aux symptômes et à la biologie. Ils ne définissent pas un régime idéal ni des compléments à prendre."),
        ("coherence", "Traits et contrôles de concordance", "Ces traits peuvent être comparés à une information connue. Ils ne certifient pas l'identité de l'échantillon ni un groupe sanguin pour une transfusion."),
        ("predisposition", "Autres prédispositions : limites de l'interprétation", "Un marqueur isolé de prédisposition n'est ni une maladie future annoncée ni un diagnostic."),
        ("performance", "Performance physique", "Un marqueur isolé ne classe pas une personne en endurance ou puissance et ne dicte aucun entraînement."),
    ):
        markers = [m for m in traits["markers"] if m["category"] == category]
        called = [m for m in markers if m["observations"]]
        sections.append(f'<section id="{category}" markdown="1"><h2>{title}</h2>'
                        f'<p class="takeaway">{len(called)} site(s) avec un enregistrement explicite sur {len(markers)} sites du panel.</p>'
                        f'<p class="caveat">{introductory}</p><div class="prose">')
        profiles = defaultdict(list)
        for marker in markers:
            profiles[marker["profile"]].append(marker)
        for items in profiles.values():
            profile = items[0]
            present = [m for m in items if m["observations"]]
            interpretable = [m for m in present if m["status"] in {"called", "called_needs_validation"}
                             and m["effect_allele_copies"] is not None]
            label = "observation provisoire ou descriptive" if present else "génotype indéterminé"
            profile_genes = " / ".join(dict.fromkeys(m["gene"] for m in items))
            sections.append(f'<details class="marker-card"{" open" if present else ""}><summary>{escape(profile_genes)} — {label}</summary>'
                            f'<p>{escape(profile["explanation_fr"])}</p>'
                            '<div class="table-wrap"><table><thead><tr><th>Marqueur</th><th>Appel observé</th><th>Qualité</th><th>Allèle ciblé</th></tr></thead><tbody>')
            for marker in items:
                observed = marker["genotype"] or "non appelé / non évaluable"
                quality = (f"DP {marker['dp']:g}, GQ {marker['gq']:g}, FILTER {marker['filter']}"
                           if marker.get("dp") is not None and marker.get("gq") is not None else "données insuffisantes")
                issues = " ; ".join(marker.get("quality_issues_fr", []))
                depths = ",".join(str(x) for x in marker["ad"]) if marker.get("ad") is not None else "non renseigné"
                status = "appel descriptif retenu" if marker["quality_pass"] else "indéterminé / validation requise"
                sections.append(f'<tr><td><strong>{escape(marker["gene"])}</strong> · {escape(marker["rsid"])}<br>{escape(marker["label"])}<br>'
                                f'<small>chr{marker["chrom"]}:{marker["pos"]} {escape(marker["ref"])}&gt;{escape(marker["alt"])}</small></td>'
                                f'<td>{escape(observed)}<br>{"site observé" if marker["observations"] else "absence ≠ référence"}</td>'
                                f'<td>{escape(quality)}<br>AD : {escape(depths)}<br>{status}<br><small>{escape(issues)}</small></td><td>{escape(marker["effect_allele"])} · '
                                f'{marker["effect_allele_copies"] if marker["effect_allele_copies"] is not None else "?"} copie(s) appelée(s)</td></tr>')
                if marker.get("marker_caveat_fr"):
                    sections.append(f'<tr><td colspan="4">{escape(marker["marker_caveat_fr"])}</td></tr>')
            sections.append('</tbody></table></div>')
            if interpretable:
                copies = "; ".join(f"{m['rsid']} : {m['effect_allele_copies']} copie(s) de {m['effect_allele']}" for m in interpretable)
                sections.append(f'<p><strong>Observation de ce fichier :</strong> {escape(copies)}. '
                                'Si elle est confirmée, elle peut être discutée dans le contexte expliqué ci-dessus.</p>')
                issues = sorted({issue for m in present for issue in m["quality_issues_fr"]})
                if issues:
                    sections.append(f'<p class="caveat">Conclusion clinique provisoire : {escape(" ; ".join(issues))}.</p>')
            sections.append(f'<p><strong>Suite possible :</strong> {escape(profile["next_steps_fr"])}</p>')
            sources = [{"title": url.split("/")[2] + " — ressource du marqueur", "url": url}
                       for url in profile["sources"] if "/refsnp/" not in url and "/snp/rs" not in url]
            sections.append(_source_links(sources, references))
            sections.append('</details>')
        sections.append('</div></section>')
    return "\n".join(sections)


def render_couple_section(couple: dict, catalog: dict, references: dict) -> str:
    """Expliquer les recouvrements incertains sans les appeler portages communs."""

    shared = couple["shared_recessive_genes"]
    established = [g for g in shared if g["disease_risks"]]
    rows = []
    for group in shared:
        unique = {}
        for pair in group["variant_pairs"]:
            for role, name in (("first", couple["first_name"]), ("second", couple["second_name"])):
                item = pair[role]
                key = (role, *normalize_allele(item["chrom"], item["pos"], item["ref"], item["alt"]))
                unique.setdefault(key, (name, item))
        for name, item in unique.values():
            change = f"{item['ref']}>{item['alt']}" if len(item["alt"]) < 30 else f"insertion / indel de longueur nette {len(item['alt']) - len(item['ref'])} bases"
            classification = ", ".join(item.get("clinvar_significances", [])) or item.get("evidence_tier", "incertain")
            rows.append(f'<tr><td>{escape(name)}</td><td>{escape(group["gene"])}</td>'
                        f'<td>chr{str(item["chrom"]).removeprefix("chr")}:{item["pos"]}<br>{escape(change)}</td>'
                        f'<td>{escape(str(item["genotype"]))}</td><td>{escape(CONSEQUENCES.get(item.get("consequence"), "effet à préciser"))}</td>'
                        f'<td>{escape(_french_classification(classification))}</td><td>{"documentée" if item.get("phased") else "non résolue"}</td></tr>')
    headline = ("Aucun portage récessif causal commun identifié dans cette analyse." if not established else
                f"{len(established)} gène(s) avec une paire récessive compatible à discuter en génétique clinique.")
    content = [f'<section id="couple-synthese"><h2>Vos résultats pour le couple</h2>'
               f'<p class="takeaway">{headline}</p>'
               '<p>Les résultats utiles de chaque personne sont expliqués plus bas. Un gène partagé avec des variants incertains '
               'ne compte pas comme une maladie récessive commune.</p>'
               '<p class="caveat">Ce résultat ne signifie pas un risque nul : la couverture et les variants non évalués laissent un risque résiduel non quantifié.</p>'
               f'<details class="audit-details"{" open" if established else ""}><summary>Comprendre les recouvrements et leur niveau de preuve</summary>'
               f'<p>{len(shared)} gène(s) commun(s) parmi les observations à revoir ; '
               f'{len(established)} gène(s) avec paire causale P/LP compatible identifiée.</p>',
               '<div class="table-wrap"><table><thead><tr><th>Personne</th><th>Gène</th><th>Variant GRCh37</th><th>GT</th><th>Effet prédit</th><th>Preuve</th><th>Phase</th></tr></thead><tbody>' + "".join(rows) + '</tbody></table></div>']
    for group in established:
        for risk in group["disease_risks"]:
            probability = risk["conditional_biallelic_probability"]
            value = f"{100 * probability:.1f} % par grossesse" if probability is not None else "indéterminée : phase ou haplotypes insuffisants"
            content.append(f'<p><strong>{escape(group["gene"])} · {escape(risk["disease_id"])} :</strong> '
                           f'probabilité conditionnelle de transmission biallélique : {value}. '
                           'Cette valeur suppose les appels et classifications valides ; elle ne signifie pas une probabilité garantie de maladie clinique.</p>')
    if any(group["gene"] == "TTN" for group in shared):
        ttn = next(group for group in shared if group["gene"] == "TTN")
        parent_counts = {
            role: len({normalize_allele(pair[role]["chrom"], pair[role]["pos"],
                                       pair[role]["ref"], pair[role]["alt"])
                       for pair in ttn["variant_pairs"]})
            for role in ("first", "second")
        }
        content.append('<div class="important-note"><strong>TTN : pas de « double récessif » démontré.</strong> '
                       f'{parent_counts["first"]} observation(s) chez {escape(couple["first_name"])} et '
                       f'{parent_counts["second"]} chez {escape(couple["second_name"])}. '
                       f'Les {len(ttn["variant_pairs"])} combinaisons de recherche réutilisent ces observations ; '
                       'elles ne sont pas des mutations supplémentaires. Un génotype 0/1 signifie un appel hétérozygote, '
                       'et ne démontre pas une double atteinte. La phase et la causalité doivent être établies. '
                       'Une paire de VUS ne justifie aucun chiffre de risque de maladie de 25 %.</div>')
        content.append('<div class="prose">' + _paragraphs(catalog["TTN"]["interpretation"]) +
                       _source_links(catalog["TTN"]["sources"], references) + '</div>')
    content.append('</details>')
    locus_rows = []
    for name, checks in couple.get("partner_locus_checks", {}).items():
        for check in checks:
            status = ("Site absent du VCF de variants : génotype et couverture indéterminés"
                      if check["status"] == "absent_from_variant_vcf_not_callable_reference"
                      else "Appel explicite : vérifier les détails dans la source de comparaison")
            locus_rows.append(f'<tr><td>{escape(name)}</td><td>{escape(check["gene"] or "trait non clinique")}</td>'
                              f'<td><code>{escape(check["variant"])}</code></td><td>{status}</td></tr>')
    if locus_rows:
        content.append('</section><section id="couple-sites" data-document-only><h2>Sites recherchés chez l’autre partenaire</h2>'
                       '<p>Ce contrôle porte sur les allèles précis, pas sur toutes les mutations possibles du gène. '
                       'Une absence de ligne ne démontre pas un résultat négatif fiable.</p>'
                       '<div class="table-wrap"><table><thead><tr><th>Personne vérifiée</th><th>Gène</th><th>Site GRCh37</th><th>Résultat</th></tr></thead><tbody>'
                       + "".join(locus_rows) + '</tbody></table></div>')
    content.append('</section><section id="couple-suites"><h2>Les vérifications qui feraient avancer le dossier</h2>'
                   '<p class="takeaway">Valider les appels et leur pertinence clinique avant de choisir une stratégie reproductive.</p>'
                   '<ol><li>Confirmer les résultats P/LP avec un laboratoire clinique et compléter l’évaluation fiable des gènes pertinents chez l’autre partenaire.</li>'
                   '<li>Examiner les lectures des indels candidats ; reconstituer ensemble les variants proches et établir la phase lorsque plusieurs observations persistent.</li>'
                   '<li>Pour les gènes partagés, vérifier les isoformes et exons pertinents, les preuves de chaque variant et le contexte clinique familial.</li>'
                   '<li>Compléter les CNV/SV, répétitions et régions non évaluables. Une absence dans un VCF de variants n’est pas un résultat négatif callable.</li></ol>'
                   '<p class="caveat">Aucune indication de PGT ou de FIV n’est déduite automatiquement de la présence de VUS. La discussion reproductive dépend des résultats confirmés et de leurs implications.</p></section>')
    return "\n".join(content)


def build_genomic_reports(
    first_result: Path,
    second_result: Path,
    couple_result: Path,
    output_directory: Path,
    *,
    skill_scripts: Path,
    refresh_traits: bool = True,
    report_date: str | None = None,
    ancestry_results: dict[str, Path] | None = None,
    polygenic_results: dict[str, Path] | None = None,
    ancestry_reference_projection: Path | None = None,
    supplementary_results: dict[str, dict[str, Path | list[Path]]] | None = None,
) -> dict[str, Path]:
    """Générer deux rapports individuels et un rapport de couple français.

    Parameters
    ----------
    first_result, second_result, couple_result : Path
        JSON issus de l'analyse à lire, sans les écraser.
    output_directory : Path
        Dossier privé, idéalement sous ``analysis/`` exclu de Git.
    skill_scripts : Path
        Dossier des fonctions ``report_builder.py`` et ``report_shell.py``.
    refresh_traits : bool
        Relire les VCF pour le panel ciblé ; False réutilise le cache local
        seulement si les empreintes des entrées, du panel et du décodeur sont identiques.
    report_date : str | None
        Date de revue au format ISO ; par défaut, date courante en UTC.
    ancestry_results, polygenic_results : dict[str, Path] | None
        JSON validés par personne, indexés par ``sample_name``. Les données
        d'ascendance et de scores sont calculées séparément du rendu ; aucun
        pourcentage n'est déduit d'une capture ou fabriqué par la présentation.
    ancestry_reference_projection : Path | None
        Nuage public des références, issu du même modèle que les résultats
        d'ascendance. Nécessaire pour leur représentation graphique.
    supplementary_results : dict[str, dict[str, Path | list[Path]]] | None
        Audits privés, indexés par personne puis par type : associations,
        exact_polygenic, dante_html, dante_comparison, pdf_reports (liste),
        pdf_comparison, ancestry_sensitivity. Les audits VCF doivent porter les mêmes empreintes.
        Les documents externes restent attribués au fournisseur, sans les
        transformer en résultats locaux ou compléter les génotypes absents.

    Returns
    -------
    dict[str, Path]
        Chemins des trois HTML autonomes et de leurs sources privées.
    """

    output_directory.mkdir(parents=True, exist_ok=True)
    report_date = report_date or datetime.now(timezone.utc).date().isoformat()
    catalog = json.loads((HERE / "clinical_extended_cards.json").read_text())
    catalog.update(json.loads((HERE / "clinical_priority_cards.json").read_text()))
    outputs = {}
    source_paths = [first_result, second_result]
    sources = [json.loads(path.read_text()) for path in source_paths]
    frozen_hashes = {path.name + ":" + path.parent.name: hashlib.sha256(path.read_bytes()).hexdigest()
                     for path in [first_result, second_result, couple_result]}
    trait_signature = hashlib.sha256((HERE / "traits.py").read_bytes() +
                                     (HERE / "trait_markers.json").read_bytes() +
                                     (HERE / "vcf.py").read_bytes()).hexdigest()
    for index, result in enumerate(sources, 1):
        name = result["sample_name"]
        folder = output_directory / f"personne_{index}"
        folder.mkdir(exist_ok=True)
        cache = folder / "traits.json"
        identities = [entry["sha256"] for entry in result["provenance"]["inputs"]]
        for entry in result["provenance"]["inputs"]:
            with Path(entry["path"]).open("rb") as handle:
                current_hash = hashlib.file_digest(handle, "sha256").hexdigest()
            if current_hash != entry["sha256"]:
                raise ValueError(f"VCF modifié depuis l'analyse : {Path(entry['path']).name}")
        if not refresh_traits and cache.exists():
            traits = json.loads(cache.read_text())
            if (traits["input_analysis_sha256"] != identities or
                traits.get("trait_definition_sha256") != trait_signature):
                raise ValueError("Cache de traits périmé : relancer avec refresh_traits=True")
        else:
            traits = audit_traits([Path(entry["path"]) for entry in result["provenance"]["inputs"]])
            traits["input_analysis_sha256"] = identities
            traits["trait_definition_sha256"] = trait_signature
            cache.write_text(json.dumps(traits, ensure_ascii=False, indent=2), encoding="utf-8")
        references = dict(GENERAL_SOURCES)
        records = collect_variant_records(result)
        counts = {tier: sum(item["report_tier"] == tier for item in records) for tier in TIERS}
        extras = {}
        for key, paths in (("ancestry", ancestry_results), ("polygenic", polygenic_results)):
            if paths and name in paths:
                path = Path(paths[name])
                value = json.loads(path.read_text())
                if value.get("sample_name") != name:
                    raise ValueError(f"Résultat complémentaire attribué à une autre personne : {key}")
                extra_inputs = value.get("inputs", value.get("input_files"))
                if not extra_inputs or [entry["sha256"] for entry in extra_inputs] != identities:
                    raise ValueError(f"Empreintes d'entrées incompatibles dans le résultat {key}")
                if value.get("analysis_result_sha256") and value["analysis_result_sha256"] != hashlib.sha256(source_paths[index-1].read_bytes()).hexdigest():
                    raise ValueError("Résultat d'ascendance lié à une autre révision d'analyse")
                extras[key] = value
        for key, paths in (supplementary_results or {}).get(name, {}).items():
            if key not in {"associations", "exact_polygenic", "dante_html", "dante_comparison", "pdf_reports", "pdf_comparison", "ancestry_sensitivity"}:
                raise ValueError(f"Résultat complémentaire non pris en charge : {key}")
            if key == "pdf_reports":
                if not isinstance(paths, list):
                    raise ValueError("pdf_reports doit être une liste de chemins JSON")
                value = [json.loads(Path(path).read_text()) for path in paths]
                for document in value:
                    original = Path(document["path"])
                    with original.open("rb") as handle:
                        current = hashlib.file_digest(handle, "sha256").hexdigest()
                    if current != document["sha256"]:
                        raise ValueError("PDF modifié depuis l'extraction")
            else:
                value = json.loads(Path(paths).read_text())
                if key in {"associations", "exact_polygenic", "dante_comparison", "pdf_comparison", "ancestry_sensitivity"}:
                    if [entry["sha256"] for entry in value.get("input_files", [])] != identities:
                        raise ValueError(f"Empreintes d'entrées incompatibles dans {key}")
                    if key in {"associations", "exact_polygenic", "ancestry_sensitivity"} and value.get("sample_name") != name:
                        raise ValueError(f"Résultat complémentaire attribué à une autre personne : {key}")
                    if key == "ancestry_sensitivity":
                        for configuration in value["configurations"]:
                            if (configuration["sample_name"] != name or
                                [entry["sha256"] for entry in configuration["inputs"]] != identities):
                                raise ValueError("Configuration de sensibilité attribuée à d'autres entrées")
            extras[key] = value
        if extras.get("dante_html") or extras.get("dante_comparison"):
            if not extras.get("dante_html") or not extras.get("dante_comparison"):
                raise ValueError("Le rapport HTML externe doit accompagner son audit de concordance")
            documents = {(doc["name"], doc["sha256"]) for doc in extras["dante_html"]["documents"]}
            if any((row["source_document"], row["source_sha256"]) not in documents
                   for row in extras["dante_comparison"]["comparisons"]):
                raise ValueError("Comparaison liée à une autre révision de document HTML")
        if extras.get("pdf_comparison"):
            documents = {(doc["filename"], doc["sha256"]) for doc in extras.get("pdf_reports", [])}
            if any((row["source_document"], row["source_sha256"]) not in documents
                   for row in extras["pdf_comparison"]["comparisons"]):
                raise ValueError("Comparaison liée à une autre révision de PDF")
        plots = {}
        if extras.get("ancestry", {}).get("status") == "reference_affinity_supported":
            if ancestry_reference_projection is None:
                raise ValueError("Le nuage des références d'ascendance est requis")
            projection = json.loads(Path(ancestry_reference_projection).read_text())
            plots["ancestry-space"] = ancestry_plot_spec(extras["ancestry"], projection)
        signals = select_patient_signals(records, traits)
        for score in extras.get("exact_polygenic", {}).get("scores", []):
            if score.get("trace_csv_path"):
                source_csv = Path(score["trace_csv_path"])
                with source_csv.open("rb") as handle:
                    fingerprint = hashlib.file_digest(handle, "sha256").hexdigest()
                if fingerprint != score["trace_csv_sha256"]:
                    raise ValueError("CSV de score modifié depuis son audit")
                pgs_id = score["pgs_id"]
                if len(pgs_id) != 9 or not pgs_id.startswith("PGS") or not pgs_id[3:].isdigit():
                    raise ValueError("Identifiant PGS invalide pour l'export")
                shutil.copyfile(source_csv, folder / (pgs_id + ".csv"))
        for domain, export in extras.get("associations", {}).get("domain_exports", {}).items():
            from .association_catalog import TRAIT_DOMAINS
            if domain not in TRAIT_DOMAINS:
                raise ValueError("Domaine GWAS inconnu pour l'export")
            source_csv = Path(export["csv_path"])
            with source_csv.open("rb") as handle:
                fingerprint = hashlib.file_digest(handle, "sha256").hexdigest()
            if fingerprint != export["sha256"]:
                raise ValueError("CSV d'associations modifié depuis son audit")
            shutil.copyfile(source_csv, folder / ("GWAS-" + domain + ".csv"))
        description = (f'{len(signals["clinical"])} résultat(s) de portage ou santé à confirmer, '
                       f'{signals["profiles"]} sujet(s) de médicaments, nutrition ou traits avec une observation utile.')
        sections = ['<section id="synthese"><h2>Vos informations clés</h2>'
                    f'<p class="takeaway">{escape(description)}</p>'
                    '<p>Commencez par les cartes : votre résultat, ce qu’il signifie et la prochaine étape. '
                    'Les détails complets sont accessibles sans encombrer cette première lecture.</p>'
                    '<p class="caveat">Ce rapport n’est pas un diagnostic. Les observations provisoires demandent une confirmation ; '
                    'les traitements et antécédents n’ont pas été renseignés.</p></section>',
                    render_panel_dashboard(build_report_panels(records, traits, catalog,
                        polygenic=extras.get("exact_polygenic") or extras.get("polygenic"), supplementary=extras)),
                    _external_coverage_section(extras, records),
                    render_patient_sections(signals, catalog, references),
                    render_population_results(extras.get("ancestry"), extras.get("polygenic"), references),
                    _ancestry_sensitivity_section(extras.get("ancestry_sensitivity"))]
        sections.append('<section id="annexes" data-document-only><h2>Détails et transparence de l’analyse</h2>'
                        f'<details class="audit-details"><summary>Ouvrir l’annexe : {len(records)} observations, qualité, VUS et sites non évalués</summary>'
                        '<p class="caveat">Les VUS et appels techniques sont conservés pour audit, pas présentés comme des maladies. '
                        'Une absence dans un VCF ne démontre pas un résultat négatif fiable.</p>'
                        + _quality_section(result, traits, extras.get("ancestry"))
                        + _structural_section(result)
                        + '<section id="variants"><h2>Catalogue détaillé de chaque observation</h2>'
                        + render_variant_catalog(records, catalog, references) + '</section>'
                        + render_trait_sections(traits, references)
                        + _population_audit(extras, references)
                        + '</details></section>')
        sections.append(_provenance_section(result, source_paths[index - 1], frozen_hashes, traits, report_date))
        md = "\n".join(sections)
        (folder / "report.md").write_text(md, encoding="utf-8")
        _save_variant_table(records, folder / "variants.csv")
        metadata = {"title": f"Vos résultats — {name}", "subtitle": "Ce qui compte pour vous, puis les détails pour comprendre et vérifier.",
                    "author": "Jouvence · analyse exploratoire et revue documentaire", "date": report_date, "version": "v4.1 · panels et sources distinctes",
                    "references": references, "analysis_input_sha256": frozen_hashes,
                    "clinical_cards_sha256": hashlib.sha256(json.dumps(catalog, sort_keys=True).encode()).hexdigest(),
                    "counts": counts, "patient_signals": {"clinical": len(signals["clinical"]), "validation": len(signals["validation"]), "profiles": signals["profiles"]},
                    "additional_results_sha256": {key: hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest() for key, value in extras.items()},
                    "clinical_history": "non renseignée"}
        structural_data = json.dumps(result.get("structural_candidates", []), ensure_ascii=False, indent=2).encode()
        for key, value in extras.items():
            (folder / (key + ".json")).write_text(json.dumps(value, ensure_ascii=False), encoding="utf-8")
        # Les audits volumineux et comparaisons fournisseur/VCF restent séparés.
        # Les panels présentent chaque source sans afficher ces comparaisons.
        embedded_extras = {key: value for key, value in extras.items()
                           if key not in {"associations", "exact_polygenic", "pdf_reports", "ancestry_sensitivity",
                                          "dante_comparison", "pdf_comparison"}}
        attachments = {"variants.csv": (folder / "variants.csv").read_bytes(), "traits.json": cache.read_bytes(),
                       "structural_candidates.json": structural_data,
                       **{key + ".json": json.dumps(value, ensure_ascii=False, indent=2).encode() for key, value in embedded_extras.items()}}
        outputs[f"personne_{index}"] = write_french_html(folder, folder / "rapport.html", metadata, skill_scripts,
                                                        attachments=attachments, plots=plots)
    couple = json.loads(couple_result.read_text())
    folder = output_directory / "couple"
    folder.mkdir(exist_ok=True)
    references = dict(GENERAL_SOURCES)
    md = render_couple_section(couple, catalog, references)
    for gene in sorted({group["gene"] for group in couple["shared_recessive_genes"]}):
        records = []
        for result in sources:
            for item in collect_variant_records(result):
                if item["gene"] == gene:
                    records.append({**item, "source_file": result["sample_name"] + " · " + item["source_file"]})
        md += f'<section id="couple-contexte-{escape(gene)}" data-document-only><h2>{escape(gene)} : détails de la comparaison</h2>' \
              + '<details class="audit-details"><summary>Afficher les observations et les maladies associées, sans diagnostic</summary>' \
              + render_variant_catalog(records, catalog, references) + '</details></section>'
    md += _couple_individual_context(sources, catalog, references)
    md += '<section id="couple-limites"><h2>Ce qui demeure non évalué</h2><p class="takeaway">Une absence de paire causale identifiée ne quantifie pas le risque résiduel.</p>'
    md += '<p class="prose">La couverture hors des sites variantés, les CNV/SV de la seconde personne, les expansions répétées et la phase distante demeurent incomplets. '
    md += 'Le volet de nouveaux variants a utilisé un transcript choisi par gène et un filtre de recherche, sans classification ACMG complète. '
    md += 'Les risques dominants ou liés à l’X exigent une interprétation variant-maladie spécifique, indépendante d’un gène commun chez les partenaires.</p></section>'
    md += f'<section id="couple-provenance" data-document-only><h2>Provenance et lecture de la comparaison</h2><p>Source : {escape(couple_result.parent.name)}. '
    md += f'Les références de chaque analyse figurent dans sa provenance ; les explications ont été révisées le {escape(report_date)}.</p>'
    md += '<p>Les combinaisons de comparaison réutilisent les observations de chaque parent ; elles ne sont pas des événements supplémentaires.</p></section>'
    (folder / "report.md").write_text(md, encoding="utf-8")
    outputs["couple"] = write_french_html(folder, folder / "rapport.html", {
        "title": f"Rapport de couple — {couple['first_name']} et {couple['second_name']}",
        "subtitle": "Une comparaison par maladie et niveau de preuve, avec les limites de phase explicites.",
        "author": "Jouvence · analyse exploratoire et revue documentaire", "date": report_date, "version": "v4.1 · panels et sources distinctes",
        "references": references, "analysis_input_sha256": frozen_hashes,
    }, skill_scripts)
    return outputs


def _ancestry_sensitivity_section(audit: dict | None) -> str:
    """Montrer la sensibilité du modèle, jamais des proportions d'ascendance."""
    if not audit:
        return ""
    panel = audit["configurations"][0]["panel"]
    fine = panel.get("fine_population_counts", {})
    html = ['<section id="ascendance-sensibilite"><h2>Ascendance : ce que change la référence enrichie</h2>',
            f'<p>Le panel enrichi contient {panel["reference_people"]} références non apparentées, '
            f'{len(fine)} populations documentées et {panel["selected_markers"]} SNP retenus. '
            f'Les petits groupes French ({fine.get("French", 0)}) et Basque ({fine.get("Basque", 0)}) décrivent des références, '
            'pas un échantillonnage représentatif de la France.</p>',
            '<div class="table-wrap"><table><thead><tr><th>Modèle testé</th><th>Contrôle indépendant des références</th>'
            '<th>Stabilité de votre classement</th><th>Contrôle par chromosomes</th><th>Conclusion autorisée</th></tr></thead><tbody>']
    for result in audit["configurations"]:
        prediction = result["provisional_model_assignment"]
        assignments = result["subsample_80_percent_assignments"]
        stable_count = assignments.get(prediction, 0)
        total = sum(assignments.values())
        folds = result["validation"]["folds"]
        conclusion = ("Affinité descriptive soutenue : " + result["closest_reference_group_label_fr"]
                      if result["status"] == "reference_affinity_supported" else "Affinité non résolue par ces contrôles")
        chromosome_checks = "; ".join(f"{key.replace('_', ' ')} : {value}" for key, value in result["chromosome_split_assignments"].items())
        html.append(f'<tr><td>{result["pca_components"]} axes PCA</td>'
                    f'<td>{100*min(f["variant_only_balanced_accuracy"] for f in folds):.1f}–'
                    f'{100*max(f["variant_only_balanced_accuracy"] for f in folds):.1f} % d’exactitude équilibrée</td>'
                    f'<td>{stable_count}/{total} sous-panels classés {escape(prediction)}</td>'
                    f'<td>{escape(chromosome_checks)}</td><td>{escape(conclusion)}</td></tr>')
    html.append('</tbody></table></div><p class="caveat">Une fréquence d’assignation à EUR ou MID n’est pas un pourcentage '
                'd’origine européenne ou moyen-orientale. Un échec de stabilité est conservé ; aucun seuil n’est abaissé '
                'pour obtenir une population souhaitée. Le nombre d’axes est une analyse de sensibilité, pas un diagnostic.</p>'
                '<p>Un classement qui change selon les chromosomes indique ici une instabilité du modèle et des données disponibles, '
                'pas la découverte d’une autre origine. La réussite des contrôles sur les références ne suffit pas si les contrôles '
                'propres à votre profil échouent.</p>'
                '<p>Les références françaises plus denses comme FranceGenRef et PREGO nécessitent une procédure d’accès ; '
                'aucun génome n’a été envoyé à un serveur d’ascendance ou d’imputation.</p>'
                '<p><a id="ancestry-sensitivity-download" href="#ancestry-sensitivity-download" '
                'data-companion-file="ancestry_sensitivity.json" download>Télécharger les contrôles complets (JSON)</a></p></section>')
    return "\n".join(html)


def _external_coverage_section(extras: dict, records: list[dict]) -> str:
    """Décrire les documents importés et leurs limites, sans comparaison de GT."""
    if not extras.get("dante_html") and not extras.get("pdf_reports"):
        return ""
    html = ['<section id="sources-fournisseur"><h2>Documents du fournisseur et périmètre des résultats</h2>',
            '<p class="takeaway">Les résultats importés restent attribués à leur source ; ils ne deviennent pas des calculs locaux validés.</p>',
            '<div class="prose"><h3>Ce que les documents renseignent</h3>']
    reports = extras.get("dante_html", {}).get("panels", [])
    if reports:
        available = [row for row in reports if row["pgs_id"]]
        empty = len(reports) - len(available)
        html.append(f'<p>Les {len(extras["dante_html"]["documents"])} HTML contiennent {len(available)} identifiants de scores et {empty} panels sans valeur ni modèle. '
                    'Leurs résultats viennent du tableau des méthodes : les compteurs animés sauvegardés peuvent rester à « 1 ». '
                    'Les scores hors ±1 σ figurent comme déclarations du fournisseur, jamais comme résultats locaux validés.</p>')
        theoretical = sum(row.get("theoretical_reference") is True for row in available)
        if theoretical:
            html.append(f'<p>Pour {theoretical} modèles, le fournisseur indique une référence théorique : fréquences publiées, '
                        'équilibre de Hardy–Weinberg et indépendance des loci. Cette approximation ne remplace pas une distribution '
                        'empirique adaptée à l’ascendance et ne tient pas nécessairement compte des corrélations entre variants.</p>')
        html.append('<p>Un percentile 94 signifie un rang de score, pas 94 % de probabilité de maladie. '
                    'Lorsque la méthode n’est pas attestée dans le document, elle reste non documentée dans cet import.</p>')
        unavailable = [row["title"] for row in reports if not row["pgs_id"]]
        if unavailable:
            html.append('<p>Panels sans résultat dans les sauvegardes : ' + escape(", ".join(unavailable)) + '.</p>')
        html.append('<p>Un classement aperçu sur une capture, sans identifiant de modèle ni poids '
                    '(cancer colorectal, par exemple), ne suffit pas à reproduire un score. '
                    'Un modèle public différent doit être étiqueté comme tel.</p>')
    pdfs = extras.get("pdf_reports", [])
    if pdfs:
        profiles = sum(1 for doc in pdfs for row in doc["findings"] if row.get("section_kind") not in {"carrier_panel", "not_currently_supported"})
        unsupported = sum(row.get("status") == "not_evaluated_by_provider" for doc in pdfs for row in doc["findings"])
        html.append(f'<p>Les PDF décrivent {profiles} profils de sommeil, sport, nutrition et médicaments. '
                    f'Le PDF pharmacogénomique présente aussi {unsupported} loci dans « Not currently supported » : '
                    'leur libellé « NOT DETECTED » signifie ici non évalué, pas un test négatif. '
                    'UGT1A1*28 et plusieurs loci HLA/CNV font notamment partie de ces limites.</p>')
        carriers = [row for doc in pdfs for row in doc["findings"] if row.get("section_kind") == "carrier_panel"]
        if carriers:
            html.append(f'<p>{len(carriers)} cartes de portage sont lisibles dans les PDF ; le nombre annoncé en couverture '
                        'n’atteste pas leur extraction complète. Le syndrome de l’X fragile et les délétions de SMN1 '
                        'demandent des tests spécialisés. Le badge « CLEAR » ne démontre pas leur exclusion par ce VCF.</p>')
    if any(row["gene"] == "TUBB8" and row["report_tier"] == "technical" for row in records):
        html.append('<h3>La maturation des ovocytes : signal trouvé, mais à confirmer</h3>'
                '<p>Un appel TUBB8 est visible dans le VCF, mais un ou plusieurs '
                'contrôles de qualité ne sont pas satisfaits : il n’est pas retenu comme génotype fiable avec nos filtres. '
                'Le panel technique et la fiche du gène expliquent les métriques, le phénotype et la confirmation nécessaire. '
                    'Cette observation ne diagnostique pas une infertilité.</p>')
    if extras.get("associations"):
        summary = extras["associations"]["summary"]
        html.append(f'<h3>Le catalogue exploratoire a été élargi</h3><p>{summary["selected_associations"]:,} associations '
                    f'publiées sélectionnées, {summary["selected_rsid_count"]:,} rsID ciblés, '
                    f'{summary["reported_trait_count"]:,} mesures ou phénotypes source et {summary["study_count"]:,} études. '
                    'Ces unités ne sont pas interchangeables : ce ne sont pas autant de maladies détectées chez vous. '
                    'Les associations présentes figurent par domaine ; les inconnus et exclusions demeurent dans l’audit.</p>')
    html.append('<details class="audit-details"><summary>Audits sources complets dans ce dossier privé</summary><ul>')
    for key in ("associations", "exact_polygenic", "pdf_reports", "dante_html"):
        if key in extras:
            html.append(f'<li><a id="audit-download-{key}" href="#audit-download-{key}" '
                        f'data-companion-file="{key}.json" download>{key}.json</a></li>')
    html.append('</ul><p>Conserver ces fichiers avec l’HTML pour l’audit complet. Les CSV des panels restent directement embarqués.</p>'
                '</details></div></section>')
    return "\n".join(html)


def write_french_html(source: Path, output: Path, metadata: dict, skill_scripts: Path,
                      attachments: dict[str, bytes] | None = None, plots: dict | None = None) -> Path:
    """Utiliser le renderer HTML documenté, puis franciser son interface.

    Le runtime partagé garde ses modes document/présentation. Les ajouts locaux
    ne font aucun appel réseau et ne changent jamais les données affichées.
    """

    sys.path.insert(0, str(skill_scripts))
    from report_builder import render_report
    base = skill_scripts.parent / "assets"
    (source / "report.css").write_text((base / "report.css").read_text() + "\n" +
        (HERE / "report_assets.css").read_text() + "\n" + (HERE / "panel_assets.css").read_text())
    runtime = (base / "report.js").read_text().replace(
        "Full screen is blocked here. Open the file in its own browser tab and press F.",
        "Le plein écran est bloqué ici. Ouvrir le fichier dans un onglet dédié puis appuyer sur F.")
    (source / "report.js").write_text(runtime + "\n" + (HERE / "report_assets.js").read_text() +
        "\n" + (HERE / "panel_assets.js").read_text())
    metadata["build_dependency_versions"] = {name: version(name) for name in ("Markdown", "beautifulsoup4")}
    metadata["presentation_build_utc"] = datetime.now(timezone.utc).isoformat()
    render_report(source, output, metadata=metadata, plots=plots or {}, attachments=attachments)
    html = output.read_text()
    replacements = {'<html lang="en">': '<html lang="fr">', 'Skip to report': 'Aller au rapport',
                    'CI <span>Research / Report</span>': 'Jouvence <span>Génomique / Rapport</span>',
                    '>Contents<': '>Sommaire<', 'Table of contents': 'Sommaire', 'Reading mode': 'Mode de lecture',
                    'Scientific report': 'Rapport documenté', '<h2>Bibliography</h2>': '<h2>Sources et lectures utiles</h2>',
                    'aria-label="Reference ': 'aria-label="Référence ',
                    'Previous slide': 'Diapositive précédente', 'Next slide': 'Diapositive suivante',
                    '← Previous': '← Précédent', 'Next →': 'Suivant →', '>Presentation<': '>Présentation<', '>Full screen<': '>Plein écran<',
                    'Enable JavaScript for interactive figures and presentation mode. The narrative and embedded source downloads remain available.':
                    'Activer JavaScript pour les filtres et la présentation. Le texte et les téléchargements intégrés restent disponibles.'}
    for old, new in replacements.items():
        html = html.replace(old, new)
    from bs4 import BeautifulSoup
    document = BeautifulSoup(html, "html.parser")
    # Le renderer commun embarque automatiquement chaque href local. Les CSV
    # de centaines de milliers de lignes restent volontairement des pièces
    # compagnons : leur lien est rétabli seulement après ce rendu documenté.
    for link in document.select("a[data-companion-file]"):
        filename = link["data-companion-file"]
        if Path(filename).name != filename or ".." in filename or not filename.endswith((".csv", ".json")):
            raise ValueError("Nom de fichier compagnon non sûr")
        if not (source / filename).is_file():
            raise ValueError(f"Fichier compagnon absent : {filename}")
        link["href"] = filename
        link["download"] = filename
    bibliography = document.select_one("#bibliography > ol")
    if bibliography:
        disclosure = document.new_tag("details", attrs={"class": "audit-details bibliography-details"})
        summary = document.new_tag("summary")
        summary.string = f"Consulter les {len(bibliography.find_all('li', recursive=False))} références"
        bibliography.wrap(disclosure)
        disclosure.insert(0, summary)
    html = str(document)
    icon = HERE.parents[1] / "assets/brand/jouvence-icon-192.png"
    if icon.exists():
        favicon = b64encode(icon.read_bytes()).decode("ascii")
        html = html.replace("</head>", f'<link rel="icon" type="image/png" href="data:image/png;base64,{favicon}"></head>')
    output.write_text(html, encoding="utf-8")
    return output


def _source_links(sources: list, references: dict) -> str:
    links = []
    for source in sources:
        if isinstance(source, str):
            source = {"url": source, "label": urlsplit(source).netloc + " — " + urlsplit(source).path.rstrip("/").split("/")[-1]}
        url = source["url"]
        key = "src-" + hashlib.sha256(url.encode()).hexdigest()[:12]
        label = source.get("title") or source.get("label") or url
        references[key] = {"url": url, "label": label}
        links.append(f'<li><a href="{escape(url, quote=True)}">{escape(label)}</a> <a data-cite="{key}"></a></li>')
    return '<ul class="sources-list">' + "".join(links) + '</ul>'


def _paragraphs(value) -> str:
    if isinstance(value, str):
        return f"<p>{escape(value)}</p>"
    return "<ul>" + "".join(f"<li>{escape(str(item))}</li>" for item in value) + "</ul>"


def _french_classification(value: str) -> str:
    return {"Pathogenic": "pathogène", "Likely pathogenic": "probablement pathogène",
            "Pathogenic/Likely pathogenic": "pathogène / probablement pathogène",
            "Uncertain significance": "signification incertaine",
            "Conflicting classifications of pathogenicity": "classifications de pathogénicité contradictoires",
            "vus": "signification incertaine", "complex_plof": "prédiction de perte de fonction sur un appel complexe",
            "novel_plof": "candidat de perte de fonction", "clinvar_plp": "P/LP dans ClinVar"}.get(value, value)


def _french_flags(flags: list[str]) -> str:
    labels = {"indel_ge_50bp_requires_read_validation": "indel ≥50 bases : validation des lectures nécessaire",
              "nearby_indel_calls_require_joint_haplotype_review": "indels voisins : reconstruire le haplotype conjointement"}
    return " ; ".join(labels.get(flag, flag.replace("raw GATK filter", "filtre GATK du fichier brut").replace("(gt ", "(> ").replace("(lt ", "(< ")) for flag in flags)


def _filter_controls() -> str:
    options = '<option value="all">Tous les niveaux de preuve</option>' + "".join(
        f'<option value="{key}">{escape(value)}</option>' for key, value in TIERS.items())
    return ('<div class="filter-tools"><label>Rechercher un gène ou une maladie '
            '<input data-gene-search type="search" placeholder="Ex. TTN, myopathie…" aria-label="Rechercher un gène ou une maladie"></label>'
            f'<label>Niveau de preuve <select data-tier-select>{options}</select></label>'
            '<button data-reset-filters>Réinitialiser</button><button data-open-cards>Tout ouvrir</button>'
            '<button data-close-cards>Tout replier</button><span class="filter-count" data-filter-count aria-live="polite"></span></div>')


def _quality_section(result: dict, traits: dict, ancestry: dict | None = None) -> str:
    rows = []
    for qc in result["qc"]:
        counts = qc["counts"]
        rows.append(f'<tr><td>{escape(Path(qc["path"]).name)}</td><td>{counts["records"]:,}</td>'
                    f'<td>{qc["depth"]["median"]}</td><td>{qc["genotype_quality"]["median"]}</td>'
                    f'<td>{counts.get("phase_block", 0):,}</td><td>{counts.get("filter_pass", 0):,} PASS / '
                    f'{counts.get("filter_not_applied", 0):,} sans filtre déclaré</td></tr>')
    unavailable = "".join(f'<li><strong>{escape(item["item"])}</strong> : {escape(item["reason_fr"])}</li>'
                          for item in traits["not_evaluated"] if not (ancestry and "Ascendance" in item["item"]))
    population_note = ("L'affinité aux références d'ascendance a été calculée séparément et ses contrôles figurent dans ce rapport. "
                       if ancestry else "L'ascendance n'a pas été calculée pour ce rendu. ")
    return ('<section id="qualite"><h2>Qualité, couverture et concordance</h2>'
            '<p class="takeaway">DP et GQ décrivent les sites appelés ; ils ne mesurent pas la couverture complète du génome.</p>'
            '<div class="table-wrap"><table><thead><tr><th>Fichier</th><th>Appels</th><th>DP médiane</th><th>GQ médiane</th><th>Sites avec bloc de phase</th><th>Filtres</th></tr></thead><tbody>'
            + "".join(rows) + '</tbody></table></div>'
            '<p class="caveat">Une absence de ligne dans ces VCF de variants demeure non évaluable. FILTER « . » signifie qu’aucun filtre n’est déclaré, même avec GQ élevé.</p>'
            '<div class="prose"><h3>Informations utiles qui manquent encore</h3><ul>' + unavailable + '</ul>'
            '<p>Les groupes sanguins sérologiques, l’histoire familiale et les traitements déclarés ne sont pas renseignés. '
            + population_note + 'Les pourcentages de mélange ancestral, la nationalité, la parenté IBD et la contamination '
            'ne sont pas déduits de ce panel.</p></div></section>')


def _population_audit(extras: dict, references: dict) -> str:
    """Rendre visibles les méthodes et exclusions, sans les confondre avec risques."""

    text = []
    ancestry = extras.get("ancestry")
    if ancestry:
        text.append('<section id="ascendance-audit"><h2>Contrôles de l’ascendance</h2>'
                    f'<p>Méthode : {escape(ancestry["method"])}. '
                    f'{ancestry["marker_counts"]["explicit_sites_accepted"]} marqueurs explicites retenus sur '
                    f'{ancestry["panel"]["selected_markers"]} sites du panel ; les autres restent inconnus.</p>'
                    + _paragraphs(ancestry.get("limitations_fr", []))
                    + _source_links(ancestry.get("sources", []), references)
                    + '<p><a href="ancestry.json">Résultat d’ascendance et contrôles complets</a></p></section>')
    polygenic = extras.get("polygenic")
    if polygenic:
        rows = []
        for score in polygenic.get("scores", []):
            state = {"incomplete_not_calibrated": "Incomplet, sans calibration personnelle",
                     "incomplete_reference_calibrated_bounds_only": "Bornes de recherche seulement",
                     "research_bounds_not_predictively_validated": "Bornes de recherche, non validées pour prédire la maladie"}.get(score["status"], score["status"])
            rows.append(f'<tr><td>{escape(score["trait_fr"])}</td><td>{escape(score["pgs_id"])}</td>'
                        f'<td>{score.get("usable_genotypes", 0)} / {score["n_variants"]}</td>'
                        f'<td>{escape(str(state))}<br>{escape(" ; ".join(score.get("reasons_fr", [])))}</td></tr>')
        text.append('<section id="scores-audit"><h2>Scores testés et raisons de non-affichage</h2>'
                    '<p>Le filtre de synthèse est |z| &gt; 1 dans une référence calibrée. '
                    'Il ne correspond pas à un seuil de signification clinique. Les génotypes absents ne sont pas transformés en référence.</p>'
                    '<p>La calibration d’un sous-ensemble de recherche peut être vérifiée numériquement sans '
                    'être validée pour prédire la maladie comme le PGS publié. Les percentiles personnels ne '
                    'sont pas affichés tant que cette qualification manque ; cela ne signifie pas un score normal.</p>'
                    '<div class="table-wrap"><table><thead><tr><th>Sujet</th><th>Modèle</th><th>Sites utilisables / attendus</th><th>État</th></tr></thead><tbody>'
                    + "".join(rows) + '</tbody></table></div>'
                    + _source_links([source for score in polygenic.get("scores", []) for source in score.get("sources", [])], references)
                    + '<p><a href="polygenic.json">Audit des scores, poids, bornes et calibration</a></p></section>')
    return "\n".join(text)


def _structural_section(result: dict) -> str:
    """Conserver les CNV/SV exploratoires, sans les compter comme maladies."""

    summary = result.get("structural_summary", {})
    candidates = result.get("structural_candidates", [])
    available = summary.get("available", False)
    message = (f'{summary.get("records", 0):,} appels structuraux dans les fichiers fournis ; '
               f'{len(candidates)} chevauchements de recherche avec des gènes de maladie, '
               f'{summary.get("strict_clinvar_matches", 0)} correspondance(s) ClinVar stricte(s).'
               if available else "Aucun fichier CNV/SV exploité pour cette personne : ce volet n’est pas évalué.")
    content = ['<section id="structure" data-document-only><h2>CNV, variants structuraux et répétitions</h2>'
               f'<p class="takeaway">{message}</p>'
               '<p class="caveat">Les appels structuraux disponibles ne comportent pas les preuves suffisantes de support des lectures, '
               'incertitudes de position et qualité du caller. Un chevauchement de gène, ou FILTER PASS seul, '
               'ne permet pas de conclure à une mutation causale ni à un second allèle pathogène.</p>']
    if candidates:
        rows = []
        for item in candidates:
            rows.append(f'<tr><td>chr{escape(str(item["chrom"]))}:{item["start"]}–{item["end"]}</td>'
                        f'<td>{escape(item["interpreted_type"])}<br>{item["length"]:,} bases</td>'
                        f'<td>{escape(", ".join(item["genes"]))}</td>'
                        f'<td>{escape(", ".join(item["exon_genes"]) or "non renseigné")}</td>'
                        f'<td>{escape(item["source_file"])}</td></tr>')
        content.append('<details><summary>Voir tous les chevauchements structuraux à valider</summary>'
                       '<div class="table-wrap"><table><thead><tr><th>Intervalle GRCh37</th><th>Type / taille</th>'
                       '<th>Gènes chevauchés</th><th>Exons chevauchés</th><th>Fichier</th></tr></thead><tbody>'
                       + "".join(rows) + '</tbody></table></div></details>'
                       '<p><a href="structural_candidates.json">Télécharger tous les appels structuraux candidats et leurs annotations</a></p>')
    repeat = result.get("repeat_analysis", {})
    content.append('<div class="prose"><h3>Expansions répétées</h3><p>' + escape(repeat.get("message") or
                   "Consulter les résultats spécialisés ; ce volet n’est pas interprété par le catalogue de petites variations.") + '</p>'
                   '<p><strong>Suite possible :</strong> récupérer les fichiers de support, examiner les lectures et faire confirmer '
                   'les événements pertinents par une méthode adaptée, après priorisation clinique. Les répétitions demandent une analyse spécialisée.</p>'
                   '<p>Ces événements techniques ne sont pas ajoutés aux compteurs des petites variations ni aux portages P/LP.</p></div></section>')
    return "\n".join(content)


def _provenance_section(result: dict, source_path: Path, hashes: dict, traits: dict, report_date: str) -> str:
    rows = "".join(f'<tr><td>{escape(Path(item["path"]).name)}</td><td><code>{item["sha256"]}</code></td></tr>'
                   for item in result["provenance"]["inputs"])
    reference_names = ", ".join(Path(item["path"]).name for item in result["provenance"].get("references", []))
    return ('<section id="provenance" data-document-only><h2>Méthode, provenance et limites</h2>'
            '<details class="audit-details"><summary>Sources, méthodes, empreintes et téléchargements</summary>'
            f'<p>Analyse figée : {escape(result["provenance"]["run_utc"])}. Lecture clinique et présentation : {escape(report_date)}.</p>'
            f'<p>Références de l’analyse : {escape(reference_names)}. '
            'Les fréquences gnomAD globales et NFE sont distinguées. NFE est une référence populationnelle européenne non finlandaise, pas une ascendance calculée ni une cohorte française.</p>'
            '<div class="table-wrap"><table><thead><tr><th>Entrée</th><th>SHA-256 enregistré dans l’analyse</th></tr></thead><tbody>' + rows + '</tbody></table></div>'
            '<p>Les petites variations connues sont confrontées à ClinVar et aux associations GenCC. Le filtre de nouveaux variants privilégie les gènes récessifs '
            'et un transcript par gène : il ne couvre pas tous les variants dominants, liés à l’X, introniques ou régulateurs.</p>'
            '<p>Les liens géniques et prédictions LoF ne sont pas des classifications ACMG complètes. Les appels complexes, les pseudogènes, la phase et la couverture '
            'requièrent des données supplémentaires. Les fréquences non renseignées ne valent pas zéro.</p>'
            f'<p>Panel ciblé de traits : {traits["summary"]["markers_tested"]} sites, coordonnées vérifiées le {escape(traits["panel_verified_at"])}. '
            'Il conserve les sites inconnus et ne calcule ni diplotype complet ni posologie.</p>'
            '<p class="download-row"><a href="variants.csv">Table complète des observations (CSV)</a>'
            '<a href="traits.json">Résultats du panel ciblé (JSON)</a><a href="report.md">Source éditable du rapport</a>'
            '<a href="build.json">Provenance du rendu HTML</a></p>'
            '<p>Aucune donnée génomique n’est envoyée à un service web par ce rapport. Les liens de documentation ouvrent des ressources externes.</p></details></section>')


def _save_variant_table(records: list[dict], destination: Path) -> None:
    fields = ["gene", "chrom", "pos", "ref", "alt", "genotype", "transcript", "consequence", "protein_change", "DP", "GQ", "allele_balance", "report_tier"]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fields, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(records)
    destination.write_text(buffer.getvalue(), encoding="utf-8")


def _couple_individual_context(results: list[dict], catalog: dict, references: dict) -> str:
    genes = sorted({x["gene"] for result in results for x in collect_variant_records(result)
                    if x["report_tier"] == "plp"})
    content = ['<section id="portages-personnels"><h2>Les résultats individuels utiles au couple</h2>'
               f'<p class="takeaway">Gènes concernés par un résultat pathogène du catalogue : {escape(", ".join(genes) or "aucun retenu")}.</p>'
               '<div class="finding-grid">']
    for result in results:
        for item in collect_variant_records(result):
            if item["report_tier"] != "plp":
                continue
            gene = item["gene"]
            card = catalog[gene]
            title = card.get("title") or card["conditions"][0]["name"]
            content.append(f'<article class="finding-card"><p class="finding-domain">{escape(result["sample_name"])} · {escape(gene)}</p>'
                           f'<h3>{escape(title)}</h3><dl class="plain-result"><dt>Votre résultat</dt>'
                           f'<dd>Un variant classé pathogène ou probablement pathogène dans le catalogue, appelé {escape(item["genotype"])} ; à confirmer cliniquement.</dd>'
                           '<dt>Ce que cela signifie pour le couple</dt><dd>' + escape(card["interpretation"]) + '</dd>'
                           '<dt>La prochaine étape</dt><dd>Confirmer le résultat et compléter le test fiable du gène chez l’autre partenaire, '
                           'plutôt que considérer l’absence de ligne VCF comme un résultat négatif.</dd></dl>'
                           '<details class="explanation prose"><summary>Symptômes, examens et prise en charge</summary>')
            content.append(_paragraphs(card["phenotype"]))
            content.append('<h4>Étapes possibles</h4>' + _paragraphs(card.get("diagnostic_steps") or card.get("diagnosis", [])))
            content.append('<h4>Prise en charge si diagnostic</h4>' + _paragraphs(card["management"]) + '</details>')
            content.append(_source_links(card["sources"], references))
            content.append('</article>')
    content.append('</div></section>')
    return "\n".join(content)
