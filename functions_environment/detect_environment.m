function env = detect_environment()
    % Detect the environment based on unique characteristics    
    % Check hostname
    [status, hostname] = system('hostname');
    hostname = strtrim(hostname); % Remove any trailing newline or spaces

    if status == 0
        if contains(hostname, 'local-hostname') || contains(hostname, 'xxxx')% Replace aaaa and xxxx with your local hostname
            env = 'localUserA';
            return;
        elseif contains(hostname, 'local-hostname') % Replace with your local hostname
                env = 'localUserB';
                return;
        else % Server hostname can vary depending on the node, etc.
            env = 'server';
            return;
        end
    end
end
