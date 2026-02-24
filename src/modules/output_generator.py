from __future__ import annotations

from typing import Optional
import re

import json
from datetime import datetime, date
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from ..config import GENBANK_FEATURE_MAP, OUTPUT_FILE_SUFFIX
from ..models.data_schemas import SequenceRecordBundle


def _flatten_qualifier_value(value):
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return [str(v) for v in value]
    return [str(value)]


def _as_location(start: int, end: int, strand: Optional[int]):
    if start < 0:
        start = 0
    if end < start:
        end = start
    if strand == 0:
        strand = None
    return FeatureLocation(start, end, strand=strand)


def _feature_qualifiers(feature):
    qualifiers: dict[str, list[str]] = {}
    qualifiers["note"] = _flatten_qualifier_value(feature.description)
    qualifiers["label"] = [feature.feature_type]
    qualifiers["source"] = [feature.source]
    if feature.score is not None:
        qualifiers["score"] = [f"{feature.score}"]
    for key, value in feature.attributes.items():
        if key == "id":
            qualifiers.setdefault("db_xref", []).append(str(value))
        else:
            qualifiers[key] = _flatten_qualifier_value(value)
    return qualifiers


def _build_record(bundle: SequenceRecordBundle) -> SeqRecord:
    seq = bundle.full_sequence.upper()
    coords = bundle.coordinates
    source_label = "ENSEMBL" if coords.coordinate_source == "ensembl" else "NCBI"

    source_db_xrefs = [f"UniProtKB:{coords.uniprot_id}"]
    if coords.coordinate_source == "ensembl" and coords.ensembl_gene_id:
        source_db_xrefs.append(f"Ensembl:{coords.ensembl_gene_id}")
    elif coords.coordinate_source != "ensembl":
        source_db_xrefs.append(f"NCBI_nuccore:{coords.ncbi_accession}")
        if coords.ensembl_gene_id and re.match(r"^ENSG|^ENS\w+G\d+", coords.ensembl_gene_id):
            source_db_xrefs.append(f"Ensembl:{coords.ensembl_gene_id}")

    record = SeqRecord(
        Seq(seq),
        id=coords.uniprot_id,
        name=coords.uniprot_id,
        description=(
            f"UniPcrTemplate gDNA region for {coords.uniprot_id} "
            f"({coords.assembly_name} {coords.seq_region_name}:{coords.ext_start_1based}-{coords.ext_end_1based}, strand={coords.strand}, source={source_label})"
        ),
    )
    record.annotations["molecule_type"] = "DNA"
    record.annotations["organism"] = coords.species
    record.annotations["taxonomy"] = [source_label]
    record.annotations["data_file_division"] = "UNC"
    record.annotations["date"] = date.today().strftime("%d-%b-%Y").upper()

    source_feature = SeqFeature(
        _as_location(0, len(seq), 1),
        type="source",
        qualifiers={
            "organism": [coords.species],
            "db_xref": source_db_xrefs,
            "note": [
                f"extracted with ±{coords.ext_end_1based - coords.ext_start_1based + 1} bp flank",
                f"original genomic: {coords.assembly_name} {coords.seq_region_name}:{coords.ext_start_1based}-{coords.ext_end_1based}",
            ],
        },
    )
    record.features.append(source_feature)

    gene_start = max(0, coords.gene_start_1based - coords.ext_start_1based)
    gene_end = max(gene_start, coords.gene_end_1based - coords.ext_start_1based + 1)
    gene_feature = SeqFeature(
        _as_location(gene_start, gene_end, coords.strand),
        type="gene",
        qualifiers={
            "gene": [coords.display_name or coords.ensembl_gene_id],
            "db_xref": (
                [f"Ensembl:{coords.ensembl_gene_id}"]
                if coords.coordinate_source == "ensembl" and coords.ensembl_gene_id
                else [f"NCBI:{coords.ncbi_accession}"] if coords.ncbi_accession else []
            ),
            "note": [f"gene span from {source_label} lookup"],
        },
    )
    record.features.append(gene_feature)

    for feature in sorted(
        bundle.features,
        key=lambda item: (item.start, -(item.end - item.start), item.feature_type),
    ):
        feature_type = feature.feature_type
        seq_feature_type = GENBANK_FEATURE_MAP.get(feature_type, "misc_feature")
        qualifiers = _feature_qualifiers(feature)
        rec_feature = SeqFeature(
            _as_location(feature.start, feature.end, feature.strand),
            type=seq_feature_type,
            qualifiers=qualifiers,
        )
        record.features.append(rec_feature)

    return record


def _feature_counts(features):
    counts: dict[str, int] = {}
    for feature in features:
        counts[feature.feature_type] = counts.get(feature.feature_type, 0) + 1
    return counts


def output_paths(outdir: Path, uniprot_id: str, assembly: str, chr_name: str, ext_start: int, ext_end: int) -> tuple[Path, Path]:
    safe_chr = "".join(ch if ch.isalnum() else "_" for ch in chr_name)
    safe_asm = "".join(ch if ch.isalnum() else "_" for ch in assembly)
    base = f"{uniprot_id}.{safe_asm}.{safe_chr}_{ext_start}_{ext_end}"
    gb_path = outdir / f"{base}{OUTPUT_FILE_SUFFIX}"
    json_path = outdir / f"{base}.metadata.json"
    return gb_path, json_path


def write_outputs(
    bundle: SequenceRecordBundle,
    outdir: Path,
    write_metadata_json: bool = True,
) -> tuple[Path, Optional[Path]]:
    outdir.mkdir(parents=True, exist_ok=True)
    gb_path, json_path = output_paths(
        outdir,
        bundle.coordinates.uniprot_id,
        bundle.coordinates.assembly_name,
        bundle.coordinates.seq_region_name,
        bundle.coordinates.ext_start_1based,
        bundle.coordinates.ext_end_1based,
    )
    record = _build_record(bundle)
    SeqIO.write(record, gb_path, "genbank")

    if write_metadata_json:
        metadata = dict(bundle.metadata)
        metadata.setdefault("run_timestamp", datetime.utcnow().isoformat() + "Z")
        metadata.setdefault("feature_counts", _feature_counts(bundle.features))
        json_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    else:
        json_path = None
    return gb_path, json_path
