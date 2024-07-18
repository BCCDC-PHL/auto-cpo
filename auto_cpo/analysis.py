import datetime
import json
import logging
import os
import shutil
import subprocess


def build_pipeline_command(config, pipeline):
    """
    Builds the pipeline command to be executed.

    :param config: The config dictionary
    :type config: dict
    :param pipeline: The pipeline dictionary
    :type pipeline: dict
    :return: The pipeline command to be executed, as a list of strings
    :rtype: list
    """

    pipeline_command = [
            'nextflow',
            '-log', pipeline['parameters']['log_path'],
            'run',
            pipeline['name'],
            '-r', pipeline['version'],
            '-profile', 'conda',
            '--cache', config['conda_cache_dir'],
            '-work-dir', pipeline['parameters']['work_dir'],
            '-with-report', pipeline['parameters']['report_path'],
            '-with-trace', pipeline['parameters']['trace_path'],
            '-with-timeline', pipeline['parameters']['timeline_path'],
    ]
    pipeline['parameters'].pop('log_path', None)
    pipeline['parameters'].pop('work_dir', None)
    pipeline['parameters'].pop('report_path', None)
    pipeline['parameters'].pop('trace_path', None)
    pipeline['parameters'].pop('timeline_path', None)

    for flag, value in pipeline['parameters'].items():
        if value is None:
            pipeline_command += ['--' + flag]
        else:
            pipeline_command += ['--' + flag, value]
    
    return pipeline_command


def run_pipeline(config, pipeline, run):
    """
    Run a pipeline.

    :param config: The config dictionary
    :type config: dict
    :param run: The run dictionary
    :type run: dict
    :param pipeline: The pipeline dictionary
    :type pipeline: dict
    :return: None
    :rtype: None
    """

    analysis_tracking = {
        "timestamp_analysis_start": datetime.datetime.now().isoformat()
    }
    pipeline_command = build_pipeline_command(config, pipeline)
    pipeline_command_str = list(map(str, pipeline_command))

    sequencing_run_id = run['sequencing_run_id']
    
    try:
        logging.info(json.dumps({
            "event_type": "analysis_started",
            "sequencing_run_id": sequencing_run_id,
            "pipeline_command": pipeline_command_str
        }))
        analysis_result = subprocess.run(pipeline_command_str, capture_output=True, check=True)
        analysis_tracking["timestamp_analysis_complete"] = datetime.datetime.now().isoformat()
        analysis_complete_path = os.path.join(pipeline['parameters']['outdir'], 'analysis_complete.json')
        with open(analysis_complete_path, 'w') as f:
                json.dump(analysis_tracking, f, indent=2)
                f.write('\n')
        logging.info(json.dumps({
            "event_type": "analysis_complete",
            "sequencing_run_id": sequencing_run_id,
            "pipeline_command": pipeline_command_str,
        }))
    except subprocess.CalledProcessError as e:
        logging.error(json.dumps({
            "event_type": "analysis_failed",
            "sequencing_run_id": sequencing_run_id,
            "pipeline_command": pipeline_command_str
        }))
