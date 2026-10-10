"""Targeted quality audit for called clinical and research candidate alleles."""

from pathlib import Path

from .vcf import normalize_allele, open_vcf, parse_info


def audit_candidate_quality(paths: list[Path], candidates: list[dict]) -> None:
    """Attach INFO quality metrics and flag raw calls using GATK hard filters.

    DP/GQ/allele-balance checks have already been applied by the matcher.
    For FILTER='.' calls, available annotations are checked against the Broad
    recommendations: SNP QD<2, FS>60, SOR>3, MQ<40, MQRankSum<-12.5,
    ReadPosRankSum<-8; indel QD<2, FS>200, SOR>10, ReadPosRankSum<-20.
    Missing metrics are recorded, and are never silently treated as passing.
    The function changes candidate dictionaries in place, never the VCF.
    """

    selected = {}
    for item in candidates:
        key = normalize_allele(item["chrom"], item["pos"], item["ref"], item["alt"])
        selected.setdefault(key, []).append(item)
    for path in paths:
        with open_vcf(path) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip().split("\t")
                chrom, pos, _, ref, raw_alt, _, filt, raw_info = fields[:8]
                for alt in raw_alt.split(","):
                    key = normalize_allele(chrom, pos, ref, alt)
                    if key not in selected:
                        continue
                    info = parse_info(raw_info)
                    snp = len(ref) == len(alt) == 1
                    rules = {"QD": ("lt", 2), "FS": ("gt", 60 if snp else 200),
                             "SOR": ("gt", 3 if snp else 10),
                             "ReadPosRankSum": ("lt", -8 if snp else -20)}
                    if snp:
                        rules.update({"MQ": ("lt", 40), "MQRankSum": ("lt", -12.5)})
                    metrics = {name: float(info[name]) for name in rules
                               if info.get(name) not in {None, ".", ""}}
                    for item in selected[key]:
                        if item["source_file"] != path.name:
                            continue
                        item["call_quality_annotations"] = metrics
                        item["missing_hard_filter_annotations"] = sorted(set(rules) - metrics.keys())
                        flags = item.setdefault("quality_flags", [])
                        if filt == ".":
                            for name, value in metrics.items():
                                direction, threshold = rules[name]
                                failed = value < threshold if direction == "lt" else value > threshold
                                if failed:
                                    flag = f"raw GATK filter: {name}={value:g} ({direction} {threshold})"
                                    if flag not in flags:
                                        flags.append(flag)
