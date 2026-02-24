from __future__ import annotations

from collections import Counter
from pathlib import Path
from dataclasses import asdict
from typing import Optional, Union

import click

from .config import (
    CACHE_DIR,
    DEFAULT_FEATURES,
    DEFAULT_FLANK,
    DEFAULT_CACHE_TTL_HOURS,
    DEFAULT_RETRIES,
    DEFAULT_TIMEOUT,
    OUTPUT_DIR,
    FeatureScanOptions,
)
from .models.data_schemas import SequenceRecordBundle
from .modules.coordinate_resolver import CoordinateResolver
from .modules.feature_scanner import FeatureScanner
from .modules.output_generator import write_outputs
from .modules.sequence_fetcher import SequenceFetcher
from .utils.api_client import ApiClient


def _parse_features(features_csv: str) -> list[str]:
    items = [item.strip() for item in features_csv.split(",") if item.strip()]
    return items if items else list(DEFAULT_FEATURES)


def run_pipeline(
    uniprot_id: str,
    outdir: Path = OUTPUT_DIR,
    flank: int = DEFAULT_FLANK,
    flank_mode: str = "genomic",
    assembly: str = "auto",
    mask: str = "soft",
    features: Union[list[str], str] = DEFAULT_FEATURES,
    maf_threshold: float = 0.01,
    gc_window: int = 50,
    gc_step: int = 10,
    gc_min: float = 30.0,
    gc_max: float = 70.0,
    homopolymer_at: int = 5,
    homopolymer_gc: int = 4,
    timeout: float = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_RETRIES,
    cache: str = "on",
    cache_ttl_hours: int = DEFAULT_CACHE_TTL_HOURS,
    offline: bool = False,
    write_metadata_json: bool = True,
    ) -> tuple[Path, Optional[Path], dict]:
    selected_features = _parse_features(",".join(features) if isinstance(features, list) else features)
    feature_options = FeatureScanOptions(
        maf_threshold=maf_threshold,
        gc_window=gc_window,
        gc_step=gc_step,
        gc_min=gc_min,
        gc_max=gc_max,
        homopolymer_at=homopolymer_at,
        homopolymer_gc=homopolymer_gc,
    )

    cache_enabled = cache == "on"
    api = ApiClient(
        timeout=timeout,
        retries=retries,
        cache_enabled=cache_enabled,
        cache_path=str(CACHE_DIR),
        ttl_hours=cache_ttl_hours,
        offline=offline,
    )

    resolver = CoordinateResolver(api)
    resolver_result = resolver.resolve(
        uniprot_id=uniprot_id,
        flank_bp=flank,
        flank_mode=flank_mode,
        assembly_preference=assembly,
    )
    coordinates = resolver_result.coordinates

    fetcher = SequenceFetcher(api)
    sequence, fetch_warnings = fetcher.fetch(coordinates=coordinates, mask=mask)

    scanner = FeatureScanner(api)
    detected_features, scan_warnings = scanner.scan(
        coordinates=coordinates,
        full_sequence=sequence,
        requested_features=selected_features,
        options=feature_options,
    )

    warnings = [*resolver_result.warnings, *fetch_warnings, *scan_warnings]
    feature_counts = Counter(f.feature_type for f in detected_features)

    metadata = {
        "uniprot_id": uniprot_id,
        "ensembl_gene_id": coordinates.ensembl_gene_id,
        "coordinate_source": coordinates.coordinate_source,
        "ncbi_accession": coordinates.ncbi_accession,
        "assembly": coordinates.assembly_name,
        "region": f"{coordinates.seq_region_name}:{coordinates.ext_start_1based}-{coordinates.ext_end_1based}:{coordinates.strand}",
        "flank_bp": flank,
        "flank_mode": flank_mode,
        "mask": mask,
        "api_cache": cache_enabled,
        "options": asdict(feature_options),
        "warnings": warnings,
        "feature_counts": dict(feature_counts),
    }

    bundle = SequenceRecordBundle(
        coordinates=coordinates,
        full_sequence=sequence,
        features=detected_features,
        metadata=metadata,
    )
    gb_path, metadata_path = write_outputs(
        bundle=bundle,
        outdir=outdir,
        write_metadata_json=write_metadata_json,
    )

    return gb_path, metadata_path, {
        "uniprot_id": uniprot_id,
        "outdir": str(outdir),
        "gb_path": str(gb_path),
        "metadata_path": str(metadata_path) if metadata_path else None,
        "n_features": len(detected_features),
        "feature_counts": dict(feature_counts),
        "warnings": warnings,
        "coordinates": coordinates.model_dump(),
        "coordinate_source": coordinates.coordinate_source,
        "ncbi_accession": coordinates.ncbi_accession,
    }


@click.command()
@click.argument("uniprot_id")
@click.option("--outdir", default=str(OUTPUT_DIR), type=click.Path(file_okay=False, path_type=Path), help="output directory")
@click.option("--flank", default=DEFAULT_FLANK, type=int, help="flanking length in bp")
@click.option(
    "--flank-mode",
    type=click.Choice(["genomic", "strand_relative"], case_sensitive=False),
    default="genomic",
    help="flank expansion mode",
)
@click.option(
    "--assembly",
    default="auto",
    help="assembly preference (GRCh38/GRCh37/auto)",
)
@click.option(
    "--mask",
    type=click.Choice(["none", "soft", "hard"]),
    default="soft",
    help="Ensembl masking mode",
)
@click.option(
    "--features",
    default=",".join(DEFAULT_FEATURES),
    help="comma separated feature types",
)
@click.option("--maf-threshold", default=0.01, type=float)
@click.option("--gc-window", default=50, type=int)
@click.option("--gc-step", default=10, type=int)
@click.option("--gc-min", default=30.0, type=float)
@click.option("--gc-max", default=70.0, type=float)
@click.option("--homopolymer-at", default=5, type=int)
@click.option("--homopolymer-gc", default=4, type=int)
@click.option("--timeout", default=DEFAULT_TIMEOUT, type=float)
@click.option("--retries", default=DEFAULT_RETRIES, type=int)
@click.option("--cache", type=click.Choice(["on", "off"]), default="on")
@click.option("--cache-ttl-hours", default=DEFAULT_CACHE_TTL_HOURS, type=int)
@click.option("--offline", is_flag=True, default=False)
@click.option("--debug", is_flag=True, default=False)
@click.option("--write-metadata-json", is_flag=True, default=True)
def cli(
    uniprot_id: str,
    outdir: Path,
    flank: int,
    flank_mode: str,
    assembly: str,
    mask: str,
    features: str,
    maf_threshold: float,
    gc_window: int,
    gc_step: int,
    gc_min: float,
    gc_max: float,
    homopolymer_at: int,
    homopolymer_gc: int,
    timeout: float,
    retries: int,
    cache: str,
    cache_ttl_hours: int,
    offline: bool,
    debug: bool,
    write_metadata_json: bool,
):
    del debug
    gb_path, metadata_path, summary = run_pipeline(
        uniprot_id=uniprot_id,
        outdir=outdir,
        flank=flank,
        flank_mode=flank_mode,
        assembly=assembly,
        mask=mask,
        features=features,
        maf_threshold=maf_threshold,
        gc_window=gc_window,
        gc_step=gc_step,
        gc_min=gc_min,
        gc_max=gc_max,
        homopolymer_at=homopolymer_at,
        homopolymer_gc=homopolymer_gc,
        timeout=timeout,
        retries=retries,
        cache=cache,
        cache_ttl_hours=cache_ttl_hours,
        offline=offline,
        write_metadata_json=write_metadata_json,
    )
    click.echo(f"GenBank: {gb_path}")
    if metadata_path:
        click.echo(f"Metadata: {metadata_path}")
    click.echo(f"Features: {summary['n_features']}")
    for feature_name, count in summary["feature_counts"].items():
        click.echo(f"- {feature_name}: {count}")


if __name__ == "__main__":
    cli()
