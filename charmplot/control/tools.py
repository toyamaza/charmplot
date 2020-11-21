import os
import yaml


dirname = os.path.dirname(__file__)


def parse_yaml(config_name):
    if not config_name.endswith(".yaml"):
        config_name = config_name + ".yaml"
    config_path = config_name
    if not os.path.isfile(config_name):
        config_path = os.path.join(dirname, "../config", config_name)
    with open(config_path, 'r') as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
