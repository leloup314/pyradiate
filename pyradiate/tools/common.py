import yaml
from pathlib import Path

from pyradiate import logger


def load_yaml(yaml_file: Path) -> dict:

    assert yaml_file.exists(), f"YAML file {str(yaml_file)} does not exist"

    with open(yaml_file, "r") as f:
        res = yaml.safe_load(f)

    return res


def save_yaml(yaml_file: Path, data: dict) -> None:

    if yaml_file.exists():
        logger.debug(f"Overwriting YAML {str(yaml_file)}")

    with open(yaml_file, "w") as f:
        yaml.safe_dump(data, f)

    assert yaml_file.exists(), f"YAML file {yaml_file} does not exist"
