"""
Tools for parsing gtf files for peak annotation.  For make new annotations
files:
    genes = get_gtf_genes(snap.genome.GRCh38.annotation)
    exons, tss = parse_gtf_features(snap.genome.GRCh38.annotation, genes)
    promoters = tss.slop(l=2000, r=100, genome='hg38')
Convert to dataframe and save as csv.
"""

from typing import Dict, List, Union
from pathlib import Path

import pybedtools
from pybedtools import BedTool


def parse_gtf_attribute(attribute_str: str) -> Dict[str, str]:
    """Parse GTF attribute string into a dictionary"""
    attributes = {}
    # Split by semicolon and process each attribute
    for attribute in attribute_str.strip(";").split(";"):
        attribute = attribute.strip()
        if "=" in attribute:
            key, value = attribute.split("=", 1)
        else:
            # Handle space-separated format
            parts = attribute.strip().split(" ", 1)
            if len(parts) != 2:
                continue
            key, value = parts
        # Clean up quotes if present
        value = value.strip("'")
        attributes[key] = value
    return attributes


def get_gtf_genes(gtf_path: Union[Path, str]) -> BedTool:
    import pandas as pd
    """Extract genes from GTF file and convert to appropriate format"""
    # Read GTF file
    gtf = pd.read_csv(gtf_path,
                      sep="\t",
                      comment="#",
                      header=None,
                      names=["chrom", "source", "feature", "start", "end",
                             "score", "strand", "frame", "attributes"])

    # Filter for genes
    genes = gtf[gtf["feature"] == "gene"].copy()

    # Parse attributes and extract gene name
    genes["gene_name"] = genes["attributes"].apply(
        lambda x: parse_gtf_attribute(x).get(
            "gene_name", parse_gtf_attribute(x).get("gene_id", "Unknown")
        )
    )

    # Convert to BED format (0-based start)
    genes_bed = pd.DataFrame({
        "chrom": genes["chrom"],
        "start": genes["start"] - 1,  # Convert to 0-based
        "end": genes["end"],
        "name": genes["gene_name"],
        "score": ".",
        "strand": genes["strand"]
    })

    return BedTool.from_dataframe(genes_bed).sort()


def parse_gtf_features(
    gtf_path: Union[Path, str], genes: BedTool
) -> List[BedTool]:
    """Extract genes, exons, and TSS from GTF file"""
    gtf = BedTool(gtf_path)

    # Extract exons
    exons = gtf.filter(lambda x: x[2] == "exon").sort().saveas()

    # Create TSS regions (gene start positions)
    tss = (
        genes.each(lambda f: create_tss_feature(f))
        .sort().saveas()
    )
    return exons, tss


def create_tss_feature(feature: BedTool) -> BedTool:
    """Create a TSS feature from a gene feature"""
    if feature.strand == '+':
        start = feature.start
        end = start + 1
    else:
        end = feature.end
        start = end - 1

    return pybedtools.create_interval_from_list([
        feature.chrom,
        str(start),
        str(end),
        feature.name,
        feature.score,
        feature.strand
    ])
