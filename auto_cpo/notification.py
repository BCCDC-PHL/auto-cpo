import argparse
import datetime
import json
import logging
import os
import uuid

from pathlib import Path

import requests
from requests.auth import HTTPBasicAuth

from jinja2 import Environment, Template
from jinja2 import BaseLoader
from importlib.resources import files

import auto_cpo.parsers as parsers
from auto_cpo.config import load_config



def _get_access_token(email_config: dict):
    """
    Get an access token from the MCMS auth service.

    :param config: A dict containing the MCMS auth service URL, client ID, and client secret.
                   Required keys are: ['auth_url', 'client_id', 'client_secret'].
    :type config: dict
    :return: A dict containing the access token and other info. Keys are: ['access_token', 'token_type', 'expires_in', 'timestamp_token_received'].
    :rtype: dict
    """
    auth_url = email_config['auth_url']
    client_id = email_config['client_id']
    client_secret = email_config['client_secret']

    auth = HTTPBasicAuth(client_id, client_secret)
    headers = {
        "Accept": "application/json",
        "Content-Type": "application/x-www-form-urlencoded",
    }
    data = {
        "client_id": client_id,
        "grant_type": "client_credentials",
    }
    response = requests.post(auth_url, data=data, headers=headers, auth=auth)
    response_json = {}
    if response.status_code == 200:
        response_json = response.json()
        timestamp = datetime.datetime.now().isoformat()
        response_json['timestamp_token_received'] = timestamp
    else:
        logging.error(json.dumps({'event_type': 'email_authentication_failed', 'status_code': response.status_code, 'message': response.text}))
        return None

    access_token = response_json['access_token']

    return access_token


def _prepare_email_body(email_data: dict, notification_config: dict):
    """
    """
    message_id = str(uuid.uuid4())
    sender_email = notification_config['sender_email']
    recipients = notification_config['recipient_email_addresses']
    sequencing_run_id = email_data['sequencing_run_id']
    subject = f"[auto-cpo] Analysis Complete: {sequencing_run_id}"

    template_path = files("auto_cpo.templates").joinpath("notification_email.html")
    template_text = template_path.read_text()

    env = Environment(loader=BaseLoader())
    template = env.from_string(template_text)
    
    body = template.render(email_data)
    
    email_request_body = {
        "messageId": message_id,
        "from": sender_email,
        "email": {
            "to": recipients,
            "subject": subject,
            "bodyType": "html",
            "body": body,
        }
    }

    return email_request_body
    

def _collect_email_data(analysis_dir):
    """
    """
    email_data = {}
    analysis_dir_dirname = os.path.dirname(analysis_dir)
    sequencing_run_id = os.path.basename(analysis_dir_dirname)
    email_data['sequencing_run_id'] = sequencing_run_id
    
    sample_qc_summary_path = Path(os.path.join(analysis_dir, f"{sequencing_run_id}_auto-cpo_sample_qc_summary.csv"))
    sample_qc_summary_by_library_id = parsers.parse_auto_cpo_sample_qc_summary(sample_qc_summary_path)

    libraries_by_library_id = {}
    for library_id, sample_qc in sample_qc_summary_by_library_id.items():
        if library_id not in libraries_by_library_id:
            libraries_by_library_id[library_id] = {'library_id': library_id}

        qc_status_fields = [
            'pre_alignment_estimated_depth_coverage_qc_status',
            'q30_percent_qc_status',
        ]
        qc_checks_passed = [sample_qc[qc_status_field] == "PASS" for qc_status_field in qc_status_fields]
        if all(qc_checks_passed):
            libraries_by_library_id[library_id]['qc_status'] = 'Pass'
        else:
            libraries_by_library_id[library_id]['qc_status'] = 'Fail'

    taxon_abundance_top_5_path = Path(os.path.join(analysis_dir, 'taxon-abundance-v0.1-output', f"{sequencing_run_id}_S_bracken_abundances_top_5.csv"))
    taxon_abundance_top_5_by_library_id = parsers.parse_bracken_abundances_top_5(taxon_abundance_top_5_path)

    for library_id, taxon_abundance_top_5 in taxon_abundance_top_5_by_library_id.items():
        most_abundant_species = taxon_abundance_top_5['abundance_1_name']
        most_abundant_percent = taxon_abundance_top_5['abundance_1_fraction_total_reads']
        if most_abundant_percent:
            most_abundant_percent = round(most_abundant_percent, 1)
        species_percentage = f"{most_abundant_species} ({most_abundant_percent}%)"
        libraries_by_library_id[library_id]['species_percentage'] = species_percentage

    
    for library_id in libraries_by_library_id.keys():
        mlst_sequence_type_path = Path(os.path.join(analysis_dir, 'mlst-nf-v0.1-output', library_id, f"{library_id}_sequence_type.csv"))
        if os.path.exists(mlst_sequence_type_path):
            sequence_type_records = parsers.parse_mlst_sequence_type(mlst_sequence_type_path)
            sequence_type = ' '.join([
                sequence_type_records[0]['scheme'],
                sequence_type_records[0]['sequence_type'],
            ])
            libraries_by_library_id[library_id]['mlst'] = sequence_type

    

    email_data['libraries'] = []
    for library_id, library_data in libraries_by_library_id.items():
        email_data['libraries'].append(library_data)

    return email_data
    

def send_notification_email(analysis_dir: Path, notification_config: dict):
    """
    """
    access_token = _get_access_token(notification_config)
    if not access_token:
        return None

    email_info = _collect_email_data(analysis_dir)
    email_body = _prepare_email_body(email_info, notification_config)
    email_url = notification_config['email_url']
    headers = {
        "Accept": "application/json",
        "Content-Type": "application/json",
        "Authorization": "Bearer " + access_token,
    }
    print(json.dumps(email_body, indent=2))

    response = requests.post(email_url, data=json.dumps(email_body), headers=headers)

    print(json.dumps(response.json(), indent=2))



def main(args):
    config = load_config(args.config)
    send_notification_email(args.analysis_outdir, config['notification'])
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--analysis-outdir')
    parser.add_argument('--config')
    args = parser.parse_args()
    main(args)
