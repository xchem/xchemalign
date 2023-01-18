import argparse
from pathlib import Path

from loguru import logger

from xchemalign.data import Dataset, SystemData


def make_data_json(data_dir: Path, output_dir: Path):
    paths = data_dir.glob("*.pdb")
    datasets = []
    for path in paths:
        dtag = path.stem
        dataset = Dataset(dtag=dtag, pdb=str(path))
        datasets.append(dataset)

    system_data: SystemData = SystemData(dataset=datasets)

    logger.info(f"Logging {len(system_data.dataset)} datasets")
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
