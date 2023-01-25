import argparse
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
from loguru import logger

from xchemalign import constants
from xchemalign.data import (
    Dataset,
    DatasetID,
    LigandBindingEvent,
    LigandBindingEvents,
    LigandID,
    SystemData,
)


def get_closest_lig(structure, coord):
    coord_gemmi = gemmi.Position(*coord)

    distances = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.name == "LIG":
                    poss = []
                    for atom in residue:
                        pos = atom.pos
                        poss.append([pos.x, pos.y, pos.z])

                    arr = np.array(poss)
                    mean = np.mean(arr, axis=0)
                    mean_pos = gemmi.Position(*mean)
                    distance = coord_gemmi.dist(mean_pos)
                    distances[residue.seqid.num] = distance

    if len(distances) == 0:
        return None

    return min(distances, key=lambda x: distances[x])


def make_data_json_from_pandda_dir(pandda_dir: Path, output_dir: Path):

    logger.info(f"PanDDA directory is: {pandda_dir}")
    logger.info(f"Output dir is: {output_dir}")

    # Get the PanDDA dirs
    analyses_dir: Path = pandda_dir / constants.PANDDA_ANALYSES_DIR
    processed_datasets_dir: Path = (
        pandda_dir / constants.PANDDA_PROCESSED_DATASETS_DIR
    )

    # Get the event table
    event_table_path: Path = (
        analyses_dir / constants.PANDDA_EVENTS_INSPECT_TABLE_PATH
    )
    event_table = pd.read_csv(event_table_path)

    # Iterate the event table, pulling out associated ligands and
    # their event maps
    initial_datasets: dict = {}
    for idx, row in event_table.iterrows():
        dtag = row["dtag"]
        event_id = row["event_idx"]
        x = row["x"]
        y = row["y"]
        z = row["z"]
        bdc = row["1-BDC"]
        ligand_confidence = row["Ligand Confidence"]

        logger.debug(f"Processing {dtag} {event_id}")

        if ligand_confidence != "High":
            logger.debug("No high confidence ligand!")
            continue

        # Get the structure
        processed_dataset_dir = processed_datasets_dir / dtag
        final_structure_dir_path = (
            processed_dataset_dir / constants.PANDDA_FINAL_STRUCTURE_PDB_DIR
        )
        final_structure_path = (
            final_structure_dir_path
            / constants.PANDDA_FINAL_STRUCTURE_PDB_TEMPLATE.format(dtag=dtag)
        )
        structure = gemmi.read_structure(str(final_structure_path))

        # Identify the closest ligand to the event
        residue_num = get_closest_lig(structure, (x, y, z))

        if not residue_num:
            continue

        # Get the event map
        xmap_path = (
            processed_dataset_dir
            / constants.PANDDA_EVENT_MAP_TEMPLATE.format(
                dtag=dtag, event_id=event_id, bdc=bdc
            )
        )

        # If dtag not already processed, add
        if dtag not in initial_datasets:
            initial_datasets[dtag] = {}

        initial_datasets[dtag][event_id] = LigandBindingEvent(
            id=event_id, residue=residue_num, xmap=str(xmap_path)
        )

    # Get the datasets
    datasets = []
    dataset_ids = []
    for dtag, events in initial_datasets.items():
        processed_dataset_dir = processed_datasets_dir / dtag
        final_structure_dir_path = (
            processed_dataset_dir / constants.PANDDA_FINAL_STRUCTURE_PDB_DIR
        )
        final_structure_path = (
            final_structure_dir_path
            / constants.PANDDA_FINAL_STRUCTURE_PDB_TEMPLATE.format(dtag=dtag)
        )
        event_ids = [
            LigandID(dtag=dtag, id=event_id) for event_id in events.keys()
        ]
        ligand_binding_events = [event for event in events.values()]
        dataset = Dataset(
            dtag=dtag,
            pdb=str(final_structure_path),
            ligand_binding_events=LigandBindingEvents(
                ligand_ids=event_ids,
                ligand_binding_events=ligand_binding_events,
            ),
        )
        dataset_ids.append(DatasetID(dtag=dtag))
        datasets.append(dataset)

    system_data: SystemData = SystemData(
        dataset_ids=dataset_ids, datasets=datasets
    )

    logger.info(f"Logging {len(system_data.datasets)} datasets")
    logger.info(f"Saveing output json to {output_dir}/data.json")
    with open(output_dir / "data.json", "w") as f:
        f.write(system_data.json())


def make_data_json(data_dir: Path, output_dir: Path):

    # Get the PanDDA dir

    #
    paths = data_dir.glob("*.pdb")
    datasets = []
    dataset_ids = []
    for path in paths:
        dtag = path.stem

        # Extract the

        dataset = Dataset(dtag=dtag, pdb=str(path))
        datasets.append(dataset)
        dataset_ids.append(DatasetID(dtag=dtag))

    system_data: SystemData = SystemData(
        dataset_ids=dataset_ids, datasets=datasets
    )

    logger.info(f"Logging {len(system_data.datasets)} datasets")
    logger.info(f"Saveing output json to {output_dir}/data.json")
    with open(output_dir / "data.json", "w") as f:
        f.write(system_data.json())


# Press the green button in the gutter to run the script.
if __name__ == "__main__":
    # Parse args
    logger.info("Parsing args...")

    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--data_dir",
        type=Path,
        help="",
    )

    parser.add_argument(
        "--output_dir",
        type=Path,
        help="",
    )

    args = parser.parse_args()

    data_dir = args.data_dir
    logger.info(f"Data dir is: {data_dir}")
    output_dir = args.output_dir
    logger.info(f"Output dir is: {output_dir}")

    make_data_json(data_dir, output_dir)
