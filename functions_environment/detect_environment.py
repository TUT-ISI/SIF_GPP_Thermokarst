import os
import socket

def detect_environment():
    #  Check hostname
    hostname = socket.gethostname()

    if 'Drakarykseni.local' in hostname:  # Replace with your local hostname
        return 'localNeus'
    elif 'local-hostname' in hostname:  # Replace with your server hostname
        return 'localAnttiL'
    else  # Server hostname can vary depending on the node, interactive seesion, etc.
        return 'server'

    raise EnvironmentError('Environment could not be detected.')
