"""Import local traçable de PDF Dante, sans interprétation médicale automatique.

Utiliser ``uv run --with pypdf==6.20.0 --with pdfplumber==0.11.10``.
Les sorties contiennent des données personnelles et doivent rester hors Git.
Les valeurs sont attribuées au fournisseur : aucun allèle de référence manquant,
score, barre colorée ou diplotype n'est reconstruit à partir d'un PDF.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections import Counter
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from urllib.parse import unquote


GENE_HEADING = re.compile(r"^([A-Z][A-Z0-9]*(?:[-/][A-Z0-9]+)*(?:\*\d+(?::\d+)*)?)(?:\s+(RS\d+(?:\s+RS\d+)*))?$", re.I)


def extract_dante_pdf_document(path: Path | str, output_directory: Path | str | None = None) -> dict:
    """Extraire texte, tableaux, liens, cartes et provenance de toutes les pages.

    Parameters
    ----------
    path : Path or str
        PDF local autorisé. Il n'est jamais envoyé à un service web.
    output_directory : Path or str, optional
        Dossier privé où écrire le JSON complet et le texte page par page.
        Sans dossier, la fonction retourne les données sans écrire de fichier.

    Returns
    -------
    dict
        Pages 1-based, texte de deux moteurs, lignes positionnées, tableaux,
        URLs annotées, cartes fournisseur et empreinte SHA-256 du PDF original.
        Les génotypes sont littéraux et leur orientation reste non vérifiée.
    """

    import pdfplumber
    from pypdf import PdfReader

    path = Path(path)
    reader = PdfReader(path)
    pages = []
    with pdfplumber.open(path) as pdf:
        if len(pdf.pages) != len(reader.pages):
            raise ValueError("Les deux lecteurs ne donnent pas le même nombre de pages.")
        for index, page in enumerate(pdf.pages):
            lines = [{"text": item["text"], "top": item["top"], "bottom": item["bottom"]}
                     for item in page.extract_text_lines()]
            links = []
            for item in page.hyperlinks:
                if not item.get("uri"):
                    continue
                bbox = (max(0, item["x0"]), max(0, item["top"]),
                        min(page.width, item["x1"]), min(page.height, item["bottom"]))
                label = page.crop(bbox).extract_text() if bbox[2] > bbox[0] and bbox[3] > bbox[1] else None
                links.append({"url": item["uri"], "title": label, "top": item["top"],
                              "bottom": item["bottom"], "bbox": list(bbox)})
            pages.append({"page": index + 1, "text": page.extract_text() or "",
                          "pypdf_text": reader.pages[index].extract_text() or "", "lines": lines,
                          "tables_raw": page.extract_tables(), "links": links})
    document = {"schema_version": 1, "source_kind": "provider_pdf_not_locally_validated",
                "path": str(path.resolve()), "filename": path.name,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "page_count": len(pages),
                "extracted_at": datetime.now(timezone.utc).isoformat(),
                "dependencies": {name: version(name) for name in ("pypdf", "pdfplumber")},
                "pages": pages, "findings": parse_dante_findings(pages),
                "references": extract_document_references(pages),
                "methods_and_limitations": [{"page": page["page"], "text": page["text"]} for page in pages
                    if re.search(r"Limitations|Disclaimers|Common alleles only|Not included|Not medical advice|Not currently supported", page["text"], re.I)],
                "score_and_method_mentions": [{"page": page["page"], "text": line["text"]}
                    for page in pages for line in page["lines"]
                    if re.search(r"single[- ](?:nucleotide|variant)|activity score|polygenic|percentile|haplotype|unphased|cis|CPIC convention", line["text"], re.I)],
                "extraction_limits_fr": [
                    "Résultats et recommandations attribués au fournisseur, sans validation locale ni diagnostic.",
                    "Un diplotype étoilé rapporté n'est pas transformé en génotypes nucléotidiques ou en rsID absents du PDF.",
                    "Les nombres de copies nucléotidiques comptent seulement les allèles littéralement écrits ; brin et assembly non vérifiés.",
                    "Les tableaux bruts et textes complets sont conservés ; aucune valeur n'est inférée depuis la couleur ou la longueur d'une barre.",
                ]}
    if output_directory is not None:
        folder = Path(output_directory)
        folder.mkdir(parents=True, exist_ok=True)
        (folder / f"{path.stem}.extraction.json").write_text(json.dumps(document, ensure_ascii=False, indent=2), encoding="utf-8")
        text = "\n\n".join(f"PAGE {page['page']}\n{page['text']}" for page in pages)
        (folder / f"{path.stem}.text.txt").write_text(text, encoding="utf-8")
    return document


def parse_dante_findings(pages: list[dict]) -> list[dict]:
    """Découper les cartes Dante par titres/gènes et marqueurs textuels explicites.

    Chaque carte conserve son bloc intégral, résultat, description et actions.
    Les liens sont rattachés par position verticale, pas par proximité médicale.
    Une page sans carte reste disponible dans le texte complet du document.
    """

    findings = []
    for page in pages:
        positioned_lines = page.get("lines", [])
        # pypdf respecte ici l'ordre logique titre/gène/statut ; pdfplumber
        # fusionne parfois le titre à gauche et le statut à droite d'une carte.
        text = page.get("pypdf_text") or page["text"]
        lines = [{"text": text, "top": index, "bottom": index + 1}
                 for index, text in enumerate(text.splitlines())]
        headings = []
        for index, line in enumerate(lines[:-1]):
            match = GENE_HEADING.fullmatch(line["text"].strip())
            following = lines[index + 1]["text"].strip()
            if match and line["text"] == line["text"].upper() and following == following.upper() and len(following) > 3:
                if any(item["text"].strip() == "YOUR FINDING" for item in lines[index + 1:]):
                    headings.append((index, match))
        for order, (index, match) in enumerate(headings):
            end = headings[order + 1][0] - 1 if order + 1 < len(headings) else len(lines)
            block = lines[index + 1:end]
            finding_at = next((i for i, item in enumerate(block) if item["text"].strip() == "YOUR FINDING"), None)
            action_at = next((i for i, item in enumerate(block) if item["text"].strip() == "WHAT YOU CAN DO"), None)
            if finding_at is None or action_at is None or finding_at >= action_at:
                continue
            source_at = next((i for i, item in enumerate(block) if item["text"].strip().startswith("SOURCES")), len(block))
            # Sans SOURCES, ne pas inclure le pied de page dans les actions.
            footer_at = next((i for i, item in enumerate(block) if item["text"].startswith(("Analysis ", "Confidential "))), len(block))
            action_end = min(source_at, footer_at)
            result = "\n".join(item["text"] for item in block[finding_at + 1:action_at])
            gene = match.group(1)
            rsids = [value.lower() for value in re.findall(r"rs\d+", (match.group(2) or "") + " " + result, re.I)]
            if re.fullmatch(r"RS\d+", gene, re.I):
                rsids.append(gene.lower())
                gene = None  # Un rsID seul n'est pas un nom de gène.
            rsids = list(dict.fromkeys(rsids))
            position_candidates = [item["top"] for item in positioned_lines if item["text"].strip() == lines[index]["text"].strip()]
            start_top = min(position_candidates) if position_candidates else lines[index]["top"]
            next_text = lines[headings[order + 1][0]]["text"].strip() if order + 1 < len(headings) else None
            next_positions = [item["top"] for item in positioned_lines if item["text"].strip() == next_text]
            end_top = min(next_positions) if next_positions else (
                lines[headings[order + 1][0]]["top"] if order + 1 < len(headings) else float("inf"))
            sources = [{key: value for key, value in link.items() if key in {"title", "url"}}
                       for link in page.get("links", []) if start_top <= link.get("top", -1) < end_top]
            item = {"page": page["page"], "title": lines[index - 1]["text"] if index else gene,
                    "gene_reported": gene, "rsids_reported": rsids,
                    "result_label_reported": block[0]["text"], "finding_reported": result,
                    "description_reported": "\n".join(value["text"] for value in block[1:finding_at]),
                    "actions_reported": "\n".join(value["text"] for value in block[action_at + 1:action_end]),
                    "source_labels_reported": block[source_at]["text"] if source_at < len(block) else None,
                    "sources": sources, "references": extract_document_references([{"page": page["page"], "links": sources, "text": result}]),
                    "bbox_top": start_top, "raw_block": "\n".join(value["text"] for value in lines[max(index - 1, 0):end]),
                    "status": "reported_by_provider_not_locally_validated", **parse_reported_genotype(result)}
            findings.append(item)
        unsupported = [index for index, line in enumerate(lines) if line["text"].strip() == "NOT DETECTED"]
        for order, index in enumerate(unsupported):
            end = unsupported[order + 1] - 2 if order + 1 < len(unsupported) else len(lines)
            block = [line["text"] for line in lines[index + 1:end]]
            block = block[:next((i for i, value in enumerate(block) if value.startswith(("PHARMACOGENOMICS", "v3d", "Confidential"))), len(block))]
            title = lines[index - 1]["text"] if index else "Locus non évalué"
            findings.append({"page": page["page"], "section_kind": "not_currently_supported",
                             "title": title, "gene_reported": title, "rsids_reported": list(dict.fromkeys(
                                 value.lower() for value in re.findall(r"rs\d+", title + " " + " ".join(block), re.I))),
                             "result_label_reported": "NOT DETECTED", "finding_reported": "\n".join(block),
                             "genotype_reported": None, "allele_counts_reported": None,
                             "status": "not_evaluated_by_provider",
                             "interpretation_limit_fr": "Le libellé NOT DETECTED est dans la rubrique Not currently supported : il ne démontre pas une absence de variant."})
        # Le PDF de portage fournit des statuts CLEAR sans génotype : ne pas en créer.
        carrier = [(index, re.fullmatch(r"([A-Z0-9-]+)\s*·\s*(MONDO:\d+)", line["text"]))
                   for index, line in enumerate(lines)]
        carrier = [(index, match) for index, match in carrier if match]
        for order, (index, match) in enumerate(carrier):
            end = carrier[order + 1][0] if order + 1 < len(carrier) else len(lines)
            block = [line["text"] for line in lines[index + 1:end]]
            if len(block) < 3:
                continue
            findings.append({"page": page["page"], "section_kind": "carrier_panel", "title": block[0],
                             "gene_reported": match.group(1), "disease_id_reported": match.group(2),
                             "inheritance_reported": block[1], "result_label_reported": block[2],
                             "finding_reported": "\n".join(block[3:]), "rsids_reported": [],
                             "genotype_reported": None, "allele_counts_reported": None,
                             "status": "reported_by_provider_not_locally_validated"})
    return findings


def parse_reported_genotype(text: str) -> dict:
    """Conserver un génotype littéral ; distinguer SNP, répétition et diplotype."""

    genotype = None
    labeled = re.search(r"\bGenotype:\s*([^\n.]+)", text, re.I)
    if labeled:
        genotype = labeled.group(1).split(",")[0].strip()
    else:
        match = re.match(r"\s*((?:[ACGT]+|\d+)(?:[/|](?:[ACGT]+|\d+))+)\b", text)
        if match:
            genotype = match.group(1)
    alleles = None
    kind = None
    if genotype:
        if re.fullmatch(r"[ACGT]+(?:[/|][ACGT]+)+", genotype):
            alleles = genotype.replace("|", "/").split("/")
            kind = "nucleotide_alleles_reported_strand_unverified"
        elif "*" in genotype:
            kind = "star_diplotype_reported_not_reconstructed"
        elif all(re.fullmatch(r"(?:rs\d+\s+)?(?:reference|variant)\s*\(([ACGT]+)\)", value, re.I)
                 for value in genotype.split("/")):
            alleles = [re.search(r"\(([ACGT]+)\)", value, re.I).group(1).upper() for value in genotype.split("/")]
            kind = "nucleotide_alleles_reported_via_labels_strand_unverified"
        elif re.fullmatch(r"\d+(?:[/|]\d+)+", genotype):
            kind = "repeat_or_numeric_alleles_reported"
        else:
            kind = "provider_genotype_label_without_nucleotide_call"
    activity = re.search(r"activity score\s*([0-9]+(?:\.[0-9]+)?)", text, re.I)
    return {"genotype_reported": genotype, "genotype_kind": kind,
            "alleles_reported": alleles, "allele_counts_reported": dict(Counter(alleles)) if alleles else None,
            "activity_score_reported": activity.group(1) if activity else None,
            "genotype_interpretation_limit_fr": "Lecture littérale du PDF, sans vérification du brin, de l'assembly, de la couverture ni de la phase."}


def extract_document_references(pages: list[dict]) -> list[dict]:
    """Extraire les URLs et identifiants explicitement présents, sans les inventer."""

    references = []
    seen = set()
    for page in pages:
        for link in page.get("links", []):
            signature = (page["page"], link["url"])
            if signature in seen:
                continue
            seen.add(signature)
            url = unquote(link["url"])
            pmid = re.search(r"(?:pubmed\.ncbi\.nlm\.nih\.gov/|ncbi\.nlm\.nih\.gov/pubmed/)(\d+)", url)
            doi = re.search(r"10\.\d{4,9}/[^\s\"<>]+", url)
            references.append({"page": page["page"], "title": link.get("title"), "url": link["url"],
                               "pmid": pmid.group(1) if pmid else None,
                               "doi": doi.group(0).rstrip(".,;") if doi else None})
        for match in re.finditer(r"\bPMID\s*[:#]?\s*(\d+)|\b(10\.\d{4,9}/[^\s\"<>]+)", page.get("text", ""), re.I):
            references.append({"page": page["page"], "pmid": match.group(1),
                               "doi": match.group(2).rstrip(".,;") if match.group(2) else None,
                               "url": None, "literal": match.group(0)})
    return references


def compare_dante_pdf_genotypes(
    documents: list[dict], vcf_paths: list[Path | str], coordinate_cache: Path | str,
    output_path: Path | str | None = None,
) -> dict:
    """Comparer les allèles littéraux des PDF à des appels VCF explicites.

    Les coordonnées GRCh37 viennent du panel local vérifié ou des placements
    chromosomiques GRCh37.p13 de dbSNP (cache public, contenant seulement la
    réponse du référentiel). Seuls les rsID publics sont interrogés ; aucun
    génotype ou fichier personnel ne sort de la machine. Le moteur de lecture
    VCF existant ``dante_html.compare_reported_genotypes`` est réutilisé.

    Le PDF ne documente pas ici l'assembly ou le brin de ses allèles : une
    comparaison est donc littérale, et non une validation clinique. Un site
    absent reste inconnu, jamais homozygote de référence. Les appels et leurs
    DP/GQ/AD/FILTER sont conservés, sans déclarer une qualité médicale acquise.
    """

    import time
    import urllib.error
    import urllib.request

    from .dante_html import compare_reported_genotypes
    from .traits import load_trait_markers

    cache = Path(coordinate_cache)
    cache.mkdir(parents=True, exist_ok=True)
    known = {marker["rsid"]: marker for marker in load_trait_markers()["markers"]}
    requested = {rsid for document in documents for finding in document["findings"]
                 if finding.get("alleles_reported") for rsid in finding.get("rsids_reported", [])}
    coordinates, lookup_issues = {}, {}
    for rsid in sorted(requested):
        if not re.fullmatch(r"rs\d+", rsid):
            raise ValueError(f"Identifiant dbSNP invalide : {rsid}")
        if rsid in known:
            marker = known[rsid]
            coordinates[rsid] = {"chrom": marker["chrom"], "pos": marker["pos"], "ref": marker["ref"],
                                 "alts": [marker["alt"]], "assembly": "GRCh37", "strand": "+",
                                 "coordinate_source": "verified_local_trait_panel", "sources": marker["sources"]}
            continue
        path = cache / f"{rsid}.refsnp.json"
        url = f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{rsid[2:]}"
        if path.exists():
            public = json.loads(path.read_text())
        else:
            public = None
            for attempt in range(3):
                try:
                    with urllib.request.urlopen(url, timeout=30) as response:
                        public = json.load(response)
                    break
                except urllib.error.HTTPError as error:
                    if error.code != 429 or attempt == 2:
                        lookup_issues[rsid] = f"dbSNP HTTP {error.code}; aucune coordonnée supposée"
                        break
                    time.sleep(attempt + 1)
                except urllib.error.URLError as error:
                    lookup_issues[rsid] = f"dbSNP inaccessible : {error.reason}; aucune coordonnée supposée"
                    break
            if public is None:
                continue
            path.write_text(json.dumps(public, ensure_ascii=False, indent=2), encoding="utf-8")
            time.sleep(.55)  # cadence de consultation publique, pas d'envoi de données privées
        if str(public.get("refsnp_id")) != rsid[2:]:
            lookup_issues[rsid] = "Identité du rsID non attestée dans la réponse/copie publique ; site non comparé"
            continue
        matches = []
        for placement in public.get("primary_snapshot_data", {}).get("placements_with_allele", []):
            annotation = placement.get("placement_annot", {})
            seq_id = placement.get("seq_id", "")
            chromosome = re.fullmatch(r"NC_0*(\d+)\.\d+", seq_id)
            traits = annotation.get("seq_id_traits_by_assembly", [])
            if annotation.get("is_mismatch") is True or not chromosome or not 1 <= int(chromosome.group(1)) <= 24 or not any(
                trait.get("assembly_name") == "GRCh37.p13" and trait.get("is_chromosome") is True
                and trait.get("is_alt") is False and trait.get("is_patch") is False for trait in traits
            ):
                continue
            alleles = [item["allele"]["spdi"] for item in placement["alleles"] if "spdi" in item.get("allele", {})]
            positions = {(item["position"], item["deleted_sequence"]) for item in alleles}
            if len(positions) != 1 or any(item["deleted_sequence"] not in "ACGT" or len(item["deleted_sequence"]) != 1
                                         or item["inserted_sequence"] not in "ACGT" or len(item["inserted_sequence"]) != 1 for item in alleles):
                continue  # Une représentation complexe n'est pas ramenée artificiellement à un SNP.
            position, reference = next(iter(positions))
            number = int(chromosome.group(1))
            matches.append({"chrom": str(number) if number <= 22 else "X" if number == 23 else "Y",
                            "pos": position + 1, "ref": reference,
                            "alts": sorted({item["inserted_sequence"] for item in alleles if item["inserted_sequence"] != reference}),
                            "assembly": "GRCh37", "strand": "+", "seq_id": seq_id,
                            "coordinate_source": "NCBI_dbSNP_GRCh37_p13_chromosomal_SPDI", "source_url": url,
                            "public_reference_sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
        if len(matches) == 1:
            coordinates[rsid] = matches[0]
        else:
            lookup_issues[rsid] = f"{len(matches)} placements chromosomiques SNP GRCh37.p13 non ambigus ; site non comparé"

    panels, excluded = [], []
    for document in documents:
        for ordinal, finding in enumerate(document["findings"], 1):
            if not finding.get("alleles_reported"):
                continue
            ids = finding.get("rsids_reported", [])
            base = {"source_document": document["filename"], "source_sha256": document["sha256"],
                    "source_page": finding["page"], "title": finding["title"],
                    "gene": finding.get("gene_reported"), "reported_genotype": "/".join(finding["alleles_reported"]),
                    "reported_result": finding["result_label_reported"], "source_genotype_literal": finding["genotype_reported"]}
            if len(ids) != 1 or ids[0] not in coordinates:
                excluded.append({**base, "rsids_reported": ids, "comparison_status": "coordinate_unverified",
                                 "reason_fr": lookup_issues.get(ids[0]) if len(ids) == 1 else "Aucun rsID unique explicitement rattaché au génotype."})
                continue
            rsid = ids[0]
            locus = coordinates[rsid]
            variant = {**base, "rsid": rsid, "chrom": locus["chrom"], "pos": locus["pos"],
                       "reference_coordinate_evidence": locus, "provider_strand": "not_documented"}
            panels.append({"id": f"{Path(document['filename']).stem}-{ordinal}", "pgs_id": None,
                           "assembly": "GRCh37", "assembly_source": "public_rsID_mapping_not_PDF_declaration",
                           "variants": [variant]})
    result = compare_reported_genotypes({"panels": panels}, list(map(Path, vcf_paths)))
    result["comparisons"].extend(excluded)
    for row in result["comparisons"]:
        observed = row.get("vcf_genotype")
        if row["comparison_status"] == "discordant" and observed:
            complemented = sorted(allele.translate(str.maketrans("ACGT", "TGCA"))
                                  for allele in row["reported_genotype"].split("/"))
            row["strand_complement_compatible"] = complemented == sorted(observed.split("/"))
        row["quality_interpretation"] = "comparison_only_not_medically_validated"
    statuses = Counter(row["comparison_status"] for row in result["comparisons"])
    result.update(counts=dict(statuses), coordinate_lookup_issues=lookup_issues,
                  source_documents=[{key: document[key] for key in ("filename", "sha256", "page_count")} for document in documents],
                  interpretation_fr="Comparaison littérale aux appels explicites. Le brin fournisseur n'est pas attesté, les absences restent inconnues et une concordance ne valide ni qualité clinique ni interprétation du trait.")
    if output_path is not None:
        target = Path(output_path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return result
