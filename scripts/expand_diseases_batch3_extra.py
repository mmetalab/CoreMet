#!/usr/bin/env python3
"""
Batch 3 supplement: 3 additional diseases to complete the 30-disease batch.
"""

import csv
import json
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DISEASE_DIR = PROJECT_ROOT / "data" / "mpidatabase" / "disease_mpi"

BATCH3_EXTRA = {

    "gastroparesis": {
        "label": "Gastroparesis",
        "category": "Gastrointestinal",
        "tissue": "Stomach",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000243", "Pyruvic acid"),
        ],
        "proteins": [
            ("P01308", "Insulin", "INS"),
            ("P06213", "Insulin receptor", "INSR"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P40926", "Malate dehydrogenase, mitochondrial", "MDH2"),
        ],
    },

    "polycythemia_vera": {
        "label": "Polycythemia Vera",
        "category": "Hematological",
        "tissue": "Blood / Bone marrow",
        "metabolites": [
            ("HMDB0000289", "Uric acid"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000177", "L-Histidine"),
        ],
        "proteins": [
            ("O60674", "Tyrosine-protein kinase JAK2", "JAK2"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P40926", "Malate dehydrogenase, mitochondrial", "MDH2"),
            ("P60174", "Triosephosphate isomerase", "TPI1"),
        ],
    },

    "thalassemia": {
        "label": "Beta-Thalassemia",
        "category": "Hematological",
        "tissue": "Blood",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000929", "L-Tryptophan"),
        ],
        "proteins": [
            ("P68871", "Hemoglobin subunit beta", "HBB"),
            ("P69905", "Hemoglobin subunit alpha", "HBA1"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P40926", "Malate dehydrogenase, mitochondrial", "MDH2"),
            ("P02768", "Serum albumin", "ALB"),
        ],
    },
}


def write_metabolites_csv(disease_dir, metabolites):
    path = disease_dir / "metabolites.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["hmdb_id", "name"])
        for hmdb_id, name in metabolites:
            w.writerow([hmdb_id, name])
    return len(metabolites)


def write_proteins_csv(disease_dir, proteins):
    path = disease_dir / "proteins.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["uniprot_id", "name", "gene"])
        for uid, name, gene in proteins:
            w.writerow([uid, name, gene])
    return len(proteins)


def main():
    sys.path.insert(0, str(PROJECT_ROOT))
    from scripts.generate_disease_predictions import (
        load_mpi_database,
        generate_predictions,
        compute_network_stats,
        compute_pathway_enrichment,
    )

    logger.info("Loading MPI database for cross-referencing...")
    known_pairs, hmdb_to_name, uniprot_to_info, pair_pathways = load_mpi_database()

    created = 0
    for folder_name, info in sorted(BATCH3_EXTRA.items()):
        disease_dir = DISEASE_DIR / folder_name
        if disease_dir.exists():
            pred_file = disease_dir / "predictions.csv"
            if pred_file.exists() and pred_file.stat().st_size > 100:
                logger.info(f"SKIP {folder_name}: already exists")
                continue

        disease_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Creating {folder_name}: {info['label']} ({info['category']})")

        n_met = write_metabolites_csv(disease_dir, info["metabolites"])
        n_prot = write_proteins_csv(disease_dir, info["proteins"])
        logger.info(f"  Wrote {n_met} metabolites, {n_prot} proteins")

        pred_df = generate_predictions(disease_dir, known_pairs, hmdb_to_name, uniprot_to_info)
        compute_network_stats(pred_df, disease_dir)
        compute_pathway_enrichment(pred_df, pair_pathways, disease_dir)

        created += 1

    logger.info(f"\nBATCH 3 SUPPLEMENT COMPLETE: {created} created")
    logger.info(f"Total diseases now: {len(list(DISEASE_DIR.iterdir()))}")


if __name__ == "__main__":
    main()
