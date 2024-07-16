import csv
import datetime
import glob
import json
import logging
import os
import re
import shutil
import subprocess
import uuid

from typing import Iterator, Optional

import auto_cpo.fastq as fastq
import auto_cpo.pre_analysis as pre_analysis
import auto_cpo.analysis as analysis
import auto_cpo.post_analysis as post_analysis


def find_fastq_dirs(config, check_symlinks_complete=True):
    """
    Find all directories in the fastq_by_run_dir that match the expected format for a sequencing run directory.

    :param config: Application config.
    :type config: dict[str, object]
    :param check_symlinks_complete: Whether or not to check for the presence of a `symlinks_complete.json` file in each run directory.
    :type check_symlinks_complete: bool
    :return: A run directory to analyze, or None
    :rtype: Iterator[Optional[dict[str, object]]]
    """
    miseq_run_id_regex = "\\d{6}_M\\d{5}_\\d+_\\d{9}-[A-Z0-9]{5}"
    nextseq_run_id_regex = "\\d{6}_VH\\d{5}_\\d+_[A-Z0-9]{9}"
    gridion_run_id_regex = "\\d{8}_\\d{4}_X[1-5]_[A-Z0-9]+_[a-z0-9]{8}"
    fastq_by_run_dir = config['fastq_by_run_dir']
    subdirs = os.scandir(fastq_by_run_dir)
    if 'analyze_runs_in_reverse_order' in config and config['analyze_runs_in_reverse_order']:
        subdirs = sorted(subdirs, key=lambda x: os.path.basename(x.path), reverse=True)
    for subdir in subdirs:
        run_id = subdir.name
        run_fastq_directory = os.path.abspath(subdir.path)

        matches_miseq_regex = re.match(miseq_run_id_regex, run_id)
        matches_nextseq_regex = re.match(nextseq_run_id_regex, run_id)
        matches_gridion_regex = re.match(gridion_run_id_regex, run_id)

        if check_symlinks_complete:
            ready_to_analyze = os.path.exists(os.path.join(subdir.path, "symlinks_complete.json"))
        else:
            ready_to_analyze = True
        conditions_checked = {
            "is_directory": subdir.is_dir(),
            "matches_illumina_run_id_format": ((matches_miseq_regex is not None) or (matches_nextseq_regex is not None)),
            "ready_to_analyze": ready_to_analyze,
        }
        conditions_met = list(conditions_checked.values())
        
        analysis_parameters = {}
        if all(conditions_met):

            logging.info(json.dumps({"event_type": "fastq_directory_found", "sequencing_run_id": run_id, "fastq_directory_path": os.path.abspath(subdir.path)}))
            analysis_parameters['fastq_input'] = run_fastq_directory
            run = {
                "sequencing_run_id": run_id,
                "fastq_directory": run_fastq_directory,
                "instrument_type": "illumina",
                "analysis_parameters": analysis_parameters
            }
            yield run
        else:
            logging.debug(json.dumps({"event_type": "directory_skipped", "fastq_directory": run_fastq_directory, "conditions_checked": conditions_checked}))
            yield None
    

def scan(config: dict[str, object]) -> Iterator[Optional[dict[str, object]]]:
    """
    Scanning involves looking for all existing runs and storing them to the database,
    then looking for all existing symlinks and storing them to the database.
    At the end of a scan, we should be able to determine which (if any) symlinks need to be created.

    :param config: Application config.
    :type config: dict[str, object]
    :return: A run directory to analyze, or None
    :rtype: Iterator[Optional[dict[str, object]]]
    """
    logging.info(json.dumps({"event_type": "scan_start"}))
    for symlinks_dir in find_fastq_dirs(config):    
        yield symlinks_dir


def check_analysis_dependencies_complete(pipeline: dict[str, object], analysis: dict[str, object], analysis_run_output_dir: str):
    """
    Check that all of the entries in the pipeline's `dependencies` config have completed. If so, return True. Return False otherwise.

    Pipeline completion is determined by the presence of an `analysis_complete.json` file in the analysis output directory.

    :param pipeline:
    :type pipeline: dict[str, object]
    :param analysis:
    :type analysis: dictp[str, object]
    :param analysis_run_output_dir:
    :type analysis_run_output_dir: str
    :return: Whether or not all of the pipelines listed in `dependencies` have completed.
    :rtype: bool
    """
    all_dependencies_complete = False
    dependencies = pipeline['dependencies']
    if dependencies is None:
        return True
    dependencies_complete = []
    dependency_infos = []
    for dependency in dependencies:
        dependency_pipeline_short_name = dependency['pipeline_name'].split('/')[1]
        dependency_pipeline_minor_version = ''.join(dependency['pipeline_version'].rsplit('.', 1)[0])
        dependency_analysis_output_dir_name = '-'.join([dependency_pipeline_short_name, dependency_pipeline_minor_version, 'output'])
        dependency_analysis_complete_path = os.path.join(analysis_run_output_dir, dependency_analysis_output_dir_name, 'analysis_complete.json')
        dependency_analysis_complete = os.path.exists(dependency_analysis_complete_path)
        dependency_info = {
            'pipeline_name': dependency['pipeline_name'],
            'pipeline_version': dependency['pipeline_version'],
            'analysis_complete_path': dependency_analysis_complete_path,
            'analysis_complete': dependency_analysis_complete
        }
        dependency_infos.append(dependency_info)
    dependencies_complete = [dep['analysis_complete'] for dep in dependency_infos]
    logging.info(json.dumps({"event_type": "checked_analysis_dependencies", "all_analysis_dependencies_complete": all(dependencies_complete), "analysis_dependencies": dependency_infos}))
    if all(dependencies_complete):
        all_dependencies_complete = True

    return all_dependencies_complete


def get_library_fastq_paths(fastq_input_dir: str):
    """
    Get the paths to all of the fastq files in a directory.
    param: fastq_input_dir: Path to a directory containing fastq files.
    type: fastq_input_dir: str
    return: Paths to R1 and R2 fastq files, indexed by library ID. Keys of the dict are library IDs, values are dicts with keys: ['ID', 'R1', 'R2'].
    rtype: dict[str, dict[str, str]]
    """
    fastq_paths_by_library_id = {}
    for fastq_file in glob.glob(os.path.join(fastq_input_dir, '*.f*q.gz')):
        fastq_file_basename = os.path.basename(fastq_file)
        fastq_file_abspath = os.path.abspath(fastq_file)
        fastq_file_basename_parts = fastq_file_basename.split('_')
        library_id = fastq_file_basename_parts[0]
        if library_id not in fastq_paths_by_library_id:
            fastq_paths_by_library_id[library_id] = {
                'ID': library_id,
                'R1': None,
                'R2': None,
            }
        if '_R1' in fastq_file_basename:
            fastq_paths_by_library_id[library_id]['R1'] = fastq_file_abspath
        elif '_R2' in fastq_file_basename:
            fastq_paths_by_library_id[library_id]['R2'] = fastq_file_abspath

    return fastq_paths_by_library_id


def analyze_run(config: dict[str, object], run: dict[str, object], analysis_type: str = "short"):
    """
    Initiate an analysis on one directory of fastq files. We assume that the directory of fastq files is named using
    a sequencing run ID.

    Runs the pipeline as defined in the config, with parameters configured for the run to be analyzed. Skips any
    analyses that have already been initiated (whether completed or not).

    Some pipelines may specify that they depend on the outputs of another through their 'dependencies' config.
    For those pipelines, we confirm that all of the upstream analyses that we depend on are complete, or
    the analysis will be skipped.

    :param config:
    :type config: dict[str, object]
    :param run: Dictionary describing the run to be analyzed. Keys: ['sequencing_run_id', 'fastq_directory', 'instrument_type', 'analysis_parameters']
    :type run: dict[str, object] Keys: ['sequencing_run_id', 'fastq_directory', 'instrument_type', 'analysis_parameters']
    :param analysis_type: The type of analysis to perform. Default is 'short', alternative is 'hybrid'.
    :type analysis_type: str
    :return: None
    :rtype: NoneType
    """
    sequencing_run_id = run['sequencing_run_id']
    for pipeline in config['pipelines']:
        try:
            logging.debug(json.dumps({"event_type": "prepare_analysis_started", "sequencing_run_id": sequencing_run_id, "pipeline_name": pipeline['name']}))
            pipeline = pre_analysis.prepare_analysis(config, pipeline, run)
        except Exception as e:
            logging.error(json.dumps({"event_type": "prepare_analysis_failed", "sequencing_run_id": sequencing_run_id, "pipeline_name": pipeline['name'], "error": str(e)}))
            return

        print(json.dumps(pipeline, indent=4))
        logging.debug(json.dumps({"event_type": "prepare_analysis_complete", "sequencing_run_id": sequencing_run_id, "pipeline_name": pipeline.get('name', "unknown")}))

        analysis_dependencies_complete = pre_analysis.check_analysis_dependencies_complete(config, pipeline, run)
        analysis_not_already_started = not os.path.exists(pipeline['parameters']['outdir'])
        conditions_checked = {
            'pipeline_dependencies_met': analysis_dependencies_complete,
            'analysis_not_already_started': analysis_not_already_started,
        }
        conditions_met = list(conditions_checked.values())

        if not all(conditions_met):
            logging.warning(json.dumps({
                "event_type": "analysis_skipped",
                "pipeline_name": pipeline['name'],
                "pipeline_version": pipeline['version'],
                "pipeline_dependencies": pipeline['dependencies'],
                "sequencing_run_id": sequencing_run_id,
                "conditions_checked": conditions_checked,
            }))
            continue

        if pipeline:
            analysis.analyze_run(config, pipeline, run)
            post_analysis.post_analysis(config, pipeline, run)
        else:
            logging.error(json.dumps({"event_type": "analysis_skipped", "sequencing_run_id": sequencing_run_id, "reason": "analysis_preparation_failed"}))
            continue
