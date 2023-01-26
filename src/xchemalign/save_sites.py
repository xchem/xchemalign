from pathlib import Path

from xchemalign import constants
from xchemalign.data import Sites


def save_sites(sites: Sites, path: Path):
    with open(path / constants.SITES_FILE_NAME, "w") as f:
        f.write(sites.json())
