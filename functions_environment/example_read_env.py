# These lines need to be added at the beggining of the codes to automatically set the environment
import json
import os
from detect_environment import detect_environment
from load_conf import load_conf

# Automatically detect environment
env = detect_environment()

# Load configuration
config = load_conf(env)

# Access input and output directories
input_dir = config['input_dir']
output_dir = config['output_dir']

print(f'Detected environment: {env}')
print(f'Using input directory: {input_dir}')
print(f'Using output directory: {output_dir}')

# Your main code here
# For example:
# process_data(input_dir, output_dir)


