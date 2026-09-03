"""
Annotation file validator and converter (GFF3 / GTF).
"""

import re
import sys
from pathlib import Path

from .MPRUN import logger


def detect_format(filepath):
    """Detect whether file is GFF3 or GTF."""
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                if "##gff-version" in line:
                    return "gff3"
                continue
            cols = line.split("\t")
            if len(cols) != 9:
                continue
            attrs = cols[8]
            if attrs.startswith("ID=") or "Parent=" in attrs:
                return "gff3"
            if "gene_id" in attrs:
                return "gtf"
    return "unknown"


def parse_gff3_attributes(attr_str):
    """Parse GFF3 attributes into dict."""
    attrs = {}
    for pair in attr_str.split(";"):
        pair = pair.strip()
        if not pair:
            continue
        if "=" in pair:
            key, val = pair.split("=", 1)
            attrs[key] = val
    return attrs


def gff3_to_gtf_attr(attrs, feature_type="CDS"):
    """Convert GFF3 attributes to GTF format."""
    gtf_attrs = {}
    if "ID" in attrs:
        gtf_attrs["gene_id"] = attrs["ID"]
    elif "Name" in attrs:
        gtf_attrs["gene_id"] = attrs["Name"]
    else:
        gtf_attrs["gene_id"] = "unknown"

    if "Parent" in attrs:
        gtf_attrs["transcript_id"] = attrs["Parent"]
    elif "ID" in attrs:
        gtf_attrs["transcript_id"] = attrs["ID"]
    else:
        gtf_attrs["transcript_id"] = gtf_attrs["gene_id"]

    for key, val in attrs.items():
        if key not in ("ID", "Name", "Parent"):
            gtf_attrs[key] = val

    return "; ".join(f'{k} "{v}"' for k, v in gtf_attrs.items())


def validate(filepath):
    """Validate annotation file. Returns (is_valid, message, records, fmt)."""
    filepath = Path(filepath)
    if not filepath.exists():
        return False, f"File not found: {filepath}", [], "unknown"

    fmt = detect_format(filepath)
    if fmt == "unknown":
        return False, "Could not detect file format. Expected GFF3 or GTF.", [], "unknown"

    records = []
    cds_count = 0
    gene_ids = set()
    total_cds_length = 0
    errors = []
    line_num = 0

    with open(filepath, "r") as f:
        for line in f:
            line_num += 1
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            cols = line.split("\t")
            if len(cols) != 9:
                errors.append(f"Line {line_num}: expected 9 columns, got {len(cols)}")
                continue

            seqid, source, ftype, start, end, score, strand, phase, attrs_str = cols

            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                errors.append(f"Line {line_num}: non-integer coordinates")
                continue

            if start_i > end_i:
                errors.append(f"Line {line_num}: start > end ({start_i} > {end_i})")
                continue
            if start_i < 1:
                errors.append(f"Line {line_num}: start < 1")
                continue

            length = end_i - start_i + 1

            if fmt == "gff3":
                attrs = parse_gff3_attributes(attrs_str)
                has_id = "ID" in attrs or "Name" in attrs
            else:
                has_id = "gene_id" in attrs_str
                attrs = {}

            if not has_id:
                errors.append(f"Line {line_num}: missing gene identifier")
                continue

            if ftype == "CDS":
                cds_count += 1
                total_cds_length += length
                if fmt == "gff3":
                    gid = attrs.get("ID", attrs.get("Name", "unknown"))
                else:
                    gid_match = re.search(r'gene_id\s+"([^"]+)"', attrs_str)
                    gid = gid_match.group(1) if gid_match else "unknown"
                gene_ids.add(gid)

            records.append({
                "seqid": seqid, "type": ftype, "start": start_i, "end": end_i,
                "length": length, "strand": strand, "attributes": attrs_str,
                "parsed_attrs": attrs
            })

    if not records:
        return False, "No valid annotation records found", [], fmt

    if cds_count == 0:
        errors.append("No CDS features found. Viral expression quantification requires CDS annotations.")

    if errors:
        return False, f"Validation failed ({len(errors)} errors):\n  " + "\n  ".join(errors[:10]), records, fmt

    msg = (
        f"Format: {fmt.upper()} | Records: {len(records)} | "
        f"CDS: {cds_count} | Genes: {len(gene_ids)} | Total CDS length: {total_cds_length} bp"
    )
    return True, msg, records, fmt


def convert_gff3_to_gtf(records, output_path):
    """Convert validated GFF3 records to GTF."""
    with open(output_path, "w") as out:
        out.write("##gff-version 2.5\n")
        out.write("# Converted from GFF3 by viralexpr\n")
        for rec in records:
            if rec["type"] != "CDS":
                continue
            attrs = rec["parsed_attrs"]
            gtf_attr = gff3_to_gtf_attr(attrs, rec["type"])
            out.write(
                f"{rec['seqid']}\t{rec.get('source', '.')}\tCDS\t"
                f"{rec['start']}\t{rec['end']}\t.\t{rec['strand']}\t.\t"
                f"{gtf_attr}\n"
            )
    logger.info("Converted GFF3 to GTF: {}".format(output_path))


def ensure_valid_gtf(input_path, output_path=None):
    """Validate annotation and return path to valid GTF."""
    input_path = Path(input_path)

    is_valid, msg, records, fmt = validate(input_path)
    if is_valid:
        logger.info(msg)
    else:
        logger.error(msg)
        sys.exit(1)

    if fmt == "gff3":
        if output_path is None:
            output_path = input_path.with_suffix(".gtf")
        convert_gff3_to_gtf(records, output_path)
        return str(output_path)
    else:
        return str(input_path)