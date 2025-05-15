import json
import os

from pathlib import Path


def load_config(config_path: Path) -> dict[str, object]:
    """
    Load auto-cpo config.

    :param config_path: Path to auto-cpo config file.
    :type config_path: Path
    :return: Parsed auto-cpo config
    :rtype: dict
    """
    config = {}
    with open(config_path, 'r') as f:
        config = json.load(f)

    if 'notification' in config:
        notification_system_config_file = config['notification'].get('system_config_file', None)
        if notification_system_config_file and os.path.exists(notification_system_config_file):
            with open(notification_system_config_file, 'r') as f:
                notification_system_config = json.load(f)
                for k, v in notification_system_config.items():
                    config['notification'][k] = v

    return config
