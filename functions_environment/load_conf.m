function config = load_config(env)
    % Read JSON file
    fid = fopen('config.json');
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    
    % Parse JSON
    configData = jsondecode(str);
    
    % Select the appropriate configuration based on the environment
    if isfield(configData, env)
        config = configData.(env);
    else
        error('Environment "%s" not found in configuration file.', env);
    end
end
