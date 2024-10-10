import json
import os

def load_config(env):
    with open('functions_environment/config.json', 'r') as f:
        config_data = json.load(f)

    if env in config_data:
        return config_data[env]
    else:
        raise ValueError(f'Environment "{env}" not found in configuration file.')