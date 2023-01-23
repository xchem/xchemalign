from pathlib import Path

from xchemalign.data import Transforms


def save_transforms(transforms: Transforms, path: Path):
    with open(path, "w") as f:
        f.write(transforms.json())
