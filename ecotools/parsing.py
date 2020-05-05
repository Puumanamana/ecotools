from pathlib import Path
import yaml

SRC_DIR = Path(__file__).resolve().parent

def parse_config():
    with open('{}/settings.yaml'.format(SRC_DIR)) as handle:
        cfg = yaml.load(handle, Loader=yaml.FullLoader)

    return cfg


