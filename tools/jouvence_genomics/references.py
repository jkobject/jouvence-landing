"""Download public references used by the no-BAM research workflow."""

from __future__ import annotations

import shutil
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


GENCODE_TRANSCRIPT_FASTA = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/"
    "gencode.v19.pc_transcripts.fa.gz"
)
TOGOVAR_GNOMAD_ROOT = (
    "https://grch37.togovar.org/downloads/release/20241203/grch37/"
    "frequency/vcf/gnomad_exomes"
)
CONTIGS = tuple(str(value) for value in range(1, 23)) + ("X", "Y")


def prepare_no_bam_references(
    reference_directory: Path,
    *,
    workers: int = 4,
) -> dict[str, object]:
    """Download the compact public references needed for VCF-only annotation.

    Existing non-empty files are reused, making the function safe to call for
    each new participant. Personal VCF data is never uploaded.
    """

    reference_directory.mkdir(parents=True, exist_ok=True)
    transcript_path = reference_directory / "gencode.v19.pc_transcripts.fa.gz"
    jobs = [(GENCODE_TRANSCRIPT_FASTA, transcript_path)]
    frequency_root = reference_directory / "gnomad_exomes_grch37"
    for population in ("all", "nfe"):
        destination = frequency_root / population
        destination.mkdir(parents=True, exist_ok=True)
        for contig in CONTIGS:
            name = f"gnomad_exomes_v2.1.1.{population}.chr{contig}.grch37.vcf.gz"
            jobs.append((f"{TOGOVAR_GNOMAD_ROOT}/{population}/{name}", destination / name))

    with ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(lambda job: _download_if_missing(*job), jobs))
    return {
        "transcript_fasta": transcript_path,
        "gnomad_frequency_directory": frequency_root,
        "downloaded_or_reused_files": len(jobs),
    }


def _download_if_missing(url: str, destination: Path) -> None:
    """Download atomically so interrupted files are never treated as complete."""

    if destination.exists() and destination.stat().st_size > 0:
        return
    temporary = destination.with_suffix(destination.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "Jouvence-research/1"})
    with urllib.request.urlopen(request, timeout=120) as response, temporary.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    temporary.replace(destination)
