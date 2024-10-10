% Main script
addpath(genpath([pwd,filesep,'functions_environment']))

% Automatically detect environment
env = detect_environment();

% Load configuration
config = load_conf(env);

% Access input and output directories
input_dir = config.input_dir;
output_dir = config.output_dir;

fprintf('Detected environment: %s\n', env);
fprintf('Using input directory: %s\n', input_dir);
fprintf('Using output directory: %s\n', output_dir);

% Your main code here