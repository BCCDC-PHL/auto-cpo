import csv
import datetime
import glob
import json
import logging
import os
import shutil

from . import parsers


def post_analysis_basic_sequence_qc(config, pipeline, run, analysis_mode):
    """
    Perform post-analysis tasks for the basic-sequence-qc pipeline. This includes generating:
    - A sample QC summary CSV file
    - A routine-assembly samplesheet
    - A taxon-abundance samplesheet

    :param config: The config dictionary
    :type config: dict
    :param pipeline: The pipeline dictionary
    :type pipeline: dict
    :param run: The run dictionary
    :type run: dict
    :return: None
    :rtype: None
    """
    logging.info(json.dumps({
        "event_type": "post_analysis_started",
        "sequencing_run_id": run['sequencing_run_id'],
        "pipeline": pipeline,
        "run": run,
    }))
    sequencing_run_id = run['sequencing_run_id']
    analysis_run_output_dir = os.path.join(config['analysis_output_dir'], sequencing_run_id, analysis_mode)
    pipeline_short_name = pipeline['name'].split('/')[1]
    pipeline_minor_version = ''.join(pipeline['version'].rsplit('.', 1)[0])
    analysis_pipeline_output_dir = pipeline.get('parameters', {}).get('outdir', None)
    logging.debug(json.dumps({
        "event_type": "analysis_pipeline_output_dir",
        "sequencing_run_id": sequencing_run_id,
        "analysis_pipeline_output_dir": analysis_pipeline_output_dir
    }))
    basic_qc_stats_csv_path = None
    basic_qc_stats_csv_path = os.path.join(
        analysis_pipeline_output_dir,
        sequencing_run_id + '_basic_qc_stats.csv'
    )
    basic_sequence_qc = {}
    if not basic_qc_stats_csv_path or not os.path.exists(basic_qc_stats_csv_path):
        logging.warning(json.dumps({
            "event_type": "basic_qc_stats_file_not_found",
            "sequencing_run_id": sequencing_run_id,
            "basic_qc_stats_csv_path": basic_qc_stats_csv_path
        }))
        return None

    logging.debug(json.dumps({
        "event_type": "basic_qc_stats_file_found",
        "sequencing_run_id": sequencing_run_id,
        "basic_qc_stats_csv_path": basic_qc_stats_csv_path
    }))
    basic_qc_stats_by_library = parsers.parse_basic_qc_stats(basic_qc_stats_csv_path)
    
    sample_qc_output_fieldnames = [
        'sequencing_run_id',
        'library_id',
        'total_bases_input',
        'filtered_bases_input',
        'pre_alignment_estimated_depth_coverage_per_mb',
        'pre_alignment_estimated_depth_coverage_qc_status',
        'q30_percent_before_filtering',
        'q30_percent_after_filtering',
        'q30_percent_qc_status',
    ]
    sample_qc_summary_path = os.path.join(
        analysis_run_output_dir,
        sequencing_run_id + "_auto-cpo_sample_qc_summary.csv"
    )
    with open(sample_qc_summary_path, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=sample_qc_output_fieldnames, dialect='unix', quoting=csv.QUOTE_MINIMAL, extrasaction='ignore')
        writer.writeheader()
        for library_id, qc_stats in basic_qc_stats_by_library.items():
            basic_qc_stats_by_library[library_id]['sequencing_run_id'] = sequencing_run_id
            coverage_pass_fail = "FAIL"
            estimated_depth_coverage = qc_stats['pre_alignment_estimated_depth_coverage_per_mb']
            if estimated_depth_coverage > config['qc_filters']['samples']['minimum_pre_alignment_estimated_depth_coverage_per_mb']:
                coverage_pass_fail = "WARN"
            if estimated_depth_coverage > config['qc_filters']['samples']['warning_pre_alignment_estimated_depth_coverage_per_mb']:
                coverage_pass_fail = "PASS"
            basic_qc_stats_by_library[library_id]['pre_alignment_estimated_depth_coverage_qc_status'] = coverage_pass_fail
            q30_pass_fail = "FAIL"
            q30_percent_after_filtering = qc_stats['q30_percent_after_filtering']
            if q30_percent_after_filtering > config['qc_filters']['samples']['minimum_q30_percent_before_filtering']:
                q30_pass_fail = "WARN"
            if q30_percent_after_filtering > config['qc_filters']['samples']['warning_q30_percent_before_filtering']:
                q30_pass_fail = "PASS"
            basic_qc_stats_by_library[library_id]['q30_percent_qc_status'] = q30_pass_fail
            
            writer.writerow(qc_stats)

    logging.info(json.dumps({
        "event_type": "qc_summary_written",
        "sequencing_run_id": sequencing_run_id,
        "qc_summary_path": sample_qc_summary_path
    }))

    samplesheets_dir = os.path.join(analysis_run_output_dir, 'samplesheets')
    os.makedirs(samplesheets_dir, exist_ok=True)

    assembly_samplesheet_path = os.path.join(samplesheets_dir, sequencing_run_id + '_routine-assembly_samplesheet.csv')

    taxon_abundance_samplesheet_path = os.path.join(samplesheets_dir, sequencing_run_id + '_taxon-abundance_samplesheet.csv')

    samplesheet_fieldnames = [
        'ID',
        'R1',
        'R2',
    ]
    fastq_by_run_dir = config['fastq_by_run_dir']
    run_fastq_dir = os.path.join(fastq_by_run_dir, sequencing_run_id)

    assembly_samplesheet_f = open(assembly_samplesheet_path, 'w')
    taxon_abundance_samplesheet_f = open(taxon_abundance_samplesheet_path, 'w')
    
    assembly_writer = csv.DictWriter(assembly_samplesheet_f, fieldnames=samplesheet_fieldnames, dialect='unix', quoting=csv.QUOTE_MINIMAL, extrasaction='ignore')
    taxon_abundance_writer = csv.DictWriter(taxon_abundance_samplesheet_f, fieldnames=samplesheet_fieldnames, dialect='unix', quoting=csv.QUOTE_MINIMAL, extrasaction='ignore')
    assembly_writer.writeheader()
    taxon_abundance_writer.writeheader()
    for library_id, qc_stats in basic_qc_stats_by_library.items():
        if qc_stats['pre_alignment_estimated_depth_coverage_qc_status'] == "FAIL" or qc_stats['q30_percent_qc_status'] == "FAIL":
            continue
        r1_fastq_glob = os.path.join(run_fastq_dir, library_id + '_R1*')
        r2_fastq_glob = os.path.join(run_fastq_dir, library_id + '_R2*')
        r1_fastq_files = glob.glob(r1_fastq_glob)
        r2_fastq_files = glob.glob(r2_fastq_glob)
        if len(r1_fastq_files) == 0 or len(r2_fastq_files) == 0:
            logging.warning(json.dumps({
                "event_type": "fastq_files_not_found",
                "sequencing_run_id": sequencing_run_id,
                "library_id": library_id,
                "r1_fastq_glob": r1_fastq_glob,
                "r2_fastq_glob": r2_fastq_glob
            }))
            continue
        if len(r1_fastq_files) > 1 or len(r2_fastq_files) > 1:
            logging.warning(json.dumps({
                "event_type": "multiple_fastq_files_found",
                "sequencing_run_id": sequencing_run_id,
                "library_id": library_id,
                "r1_fastq_files": r1_fastq_files,
                "r2_fastq_files": r2_fastq_files
            }))
            continue
        r1_fastq_path = os.path.abspath(r1_fastq_files[0])
        r2_fastq_path = os.path.abspath(r2_fastq_files[0])
        assembly_samplesheet = {
            'ID': library_id,
            'R1': r1_fastq_path,
            'R2': r2_fastq_path,
        }
        assembly_writer.writerow(assembly_samplesheet)
        taxon_abundance_writer.writerow(assembly_samplesheet)

    return None


def post_analysis_taxon_abundance(config, pipeline, run):
    """
    Perform post-analysis tasks for the taxon-abundance pipeline.

    :param config: The config dictionary
    :type config: dict
    :param pipeline: The pipeline dictionary
    :type pipeline: dict
    :param run: The run dictionary
    :type run: dict
    :return: None
    :rtype: None
    """
    logging.info(json.dumps({
        "event_type": "post_analysis_started",
        "sequencing_run_id": run['sequencing_run_id'],
        "pipeline": pipeline,
        "run": run,
    }))
    sequencing_run_id = run['sequencing_run_id']
    analysis_run_output_dir = os.path.join(config['analysis_output_dir'], sequencing_run_id)

    return None


def post_analysis_routine_assembly(config, pipeline, run):
    """
    Perform post-analysis tasks for the routine-assembly pipeline.

    :param config: The config dictionary
    :type config: dict
    :param pipeline: The pipeline dictionary
    :type pipeline: dict
    :param run: The run dictionary
    :type run: dict
    :return: None
    :rtype: None
    """

    logging.info(json.dumps({
        "event_type": "post_analysis_started",
        "sequencing_run_id": run['sequencing_run_id'],
        "pipeline": pipeline,
        "run": run,
    }))
    sequencing_run_id = run['sequencing_run_id']
    analysis_run_output_dir = os.path.join(config['analysis_output_dir'], sequencing_run_id)

    assembly_symlinks_dir = os.path.join(config['analysis_output_dir'], sequencing_run_id, 'assemblies')
    os.makedirs(assembly_symlinks_dir, exist_ok=True)

    pipeline_short_name = pipeline['name'].split('/')[1]
    pipeline_minor_version = ''.join(pipeline['version'].rsplit('.', 1)[0])
    analysis_pipeline_output_dir = pipeline.get('parameters', {}).get('outdir', None)
    logging.debug(json.dumps({
        "event_type": "analysis_pipeline_output_dir",
        "sequencing_run_id": sequencing_run_id,
        "analysis_pipeline_output_dir": analysis_pipeline_output_dir
    }))
    assemblies_glob = os.path.join(analysis_pipeline_output_dir, 'BC*', 'BC*_unicycler_short.fa')
    assembly_paths = glob.glob(assemblies_glob)
    if len(assembly_paths) == 0:
        logging.warning(json.dumps({
            "event_type": "no_assemblies_found",
            "sequencing_run_id": sequencing_run_id,
            "assemblies_glob": assemblies_glob
        }))
        return None
    for assembly_path in assembly_paths:
        assembly_basename = os.path.basename(assembly_path)
        assembly_symlink_path = os.path.join(assembly_symlinks_dir, assembly_basename)
        try:
            os.symlink(assembly_path, assembly_symlink_path)
        except OSError as e:
            logging.error(json.dumps({
                "event_type": "symlink_assembly_failed",
                "sequencing_run_id": sequencing_run_id,
                "assembly_path": assembly_path,
                "assembly_symlink_path": assembly_symlink_path
            }))
            continue

    samplesheets_dir = os.path.join(analysis_run_output_dir, 'samplesheets')
    os.makedirs(samplesheets_dir, exist_ok=True)

    mlst_nf_samplesheet_fieldnames = [
        'ID',
        'ASSEMBLY',
    ]

    mlst_nf_samplesheet = []
    for assembly_path in assembly_paths:
        assembly_basename = os.path.basename(assembly_path)
        sample_id = assembly_basename.split('_')[0]
        mlst_nf_samplesheet.append({
            'ID': sample_id,
            'ASSEMBLY': os.path.abspath(assembly_path),
        })

    mlst_nf_samplesheet_path = os.path.join(samplesheets_dir, sequencing_run_id + '_mlst-nf_samplesheet.csv')
    with open(mlst_nf_samplesheet_path, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=mlst_nf_samplesheet_fieldnames, dialect='unix', quoting=csv.QUOTE_MINIMAL, extrasaction='ignore')
        writer.writeheader()
        for record in mlst_nf_samplesheet:
            writer.writerow(record)

    logging.info(json.dumps({"event_type": "mlst_nf_samplesheet_written", "sequencing_run_id": sequencing_run_id, "mlst_nf_samplesheet_path": mlst_nf_samplesheet_path}))

    plasmid_screen_samplesheet = []
    assembly_samplesheet_path = os.path.join(samplesheets_dir, sequencing_run_id + '_routine-assembly_samplesheet.csv')
    with open(assembly_samplesheet_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            plasmid_screen_samplesheet.append(row)

    for plasmid_screen_samplesheet_record in plasmid_screen_samplesheet:
        sample_id = record['ID']
        for mlst_nf_samplesheet_record in mlst_nf_samplesheet:
            if mlst_nf_samplesheet_record['ID'] == sample_id:
                plasmid_screen_samplesheet_record['ASSEMBLY'] = mlst_nf_samplesheet_record['ASSEMBLY']
                break

    plasmid_screen_samplesheet_fieldnames = [
        'ID',
        'R1',
        'R2',
        'ASSEMBLY',
    ]

    plasmid_screen_samplesheet_path = os.path.join(samplesheets_dir, sequencing_run_id + '_plasmid-screen_samplesheet.csv')
    with open(plasmid_screen_samplesheet_path, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=plasmid_screen_samplesheet_fieldnames, dialect='unix', quoting=csv.QUOTE_MINIMAL, extrasaction='ignore')
        writer.writeheader()
        for record in plasmid_screen_samplesheet:
            writer.writerow(record)

    return None


def post_analysis_mlst_nf(config, pipeline, run):
    """
    Perform post-analysis tasks for the mlst-nf pipeline.

    :param config: The config dictionary
    :type config: dict
    :param pipeline: The pipeline dictionary
    :type pipeline: dict
    :param run: The run dictionary
    :type run: dict
    :return: None
    :rtype: None
    """
    logging.info(json.dumps({
        "event_type": "post_analysis_started",
        "sequencing_run_id": run['sequencing_run_id'],
        "pipeline": pipeline,
        "run": run,
    }))
    sequencing_run_id = run['sequencing_run_id']
    analysis_run_output_dir = os.path.join(config['analysis_output_dir'], sequencing_run_id)

    return None


def post_analysis_plasmid_screen(config, pipeline, run, analysis_mode):
    """
    Perform post-analysis tasks for the plasmid-screen pipeline.

    :param config: The config dictionary
    :type config: dict
    :param pipeline: The pipeline dictionary
    :type pipeline: dict
    :param run: The run dictionary
    :type run: dict
    :return: None
    :rtype: None
    """
    logging.info(json.dumps({
        "event_type": "post_analysis_started",
        "sequencing_run_id": run['sequencing_run_id'],
        "pipeline": pipeline,
        "run": run,
    }))
    sequencing_run_id = run['sequencing_run_id']
    analysis_run_output_dir = os.path.join(config['analysis_output_dir'], sequencing_run_id)

    return None


def post_analysis(config, pipeline, run, analysis_mode):
    """
    Perform post-analysis tasks for a pipeline.

    :param config: The config dictionary
    :type config: dict
    :param pipeline: The pipeline dictionary
    :type pipeline: dict
    :param run: The run dictionary
    :type run: dict
    :param analysis_mode: The analysis mode
    :type analysis_mode: str
    :return: None
    """
    pipeline_name = pipeline['name']
    pipeline_short_name = pipeline_name.split('/')[1]
    pipeline_version = pipeline['version']
    delete_pipeline_work_dir = pipeline.get('delete_work_dir', True)
    sequencing_run_id = run['sequencing_run_id']
    base_analysis_work_dir = config['analysis_work_dir']

    # The work_dir includes a timestamp, so we need to glob to find the most recent one
    work_dir_glob = os.path.join(base_analysis_work_dir, 'work-' + sequencing_run_id + '_' + pipeline_short_name + '_' + '*')
    work_dirs = glob.glob(work_dir_glob)
    if len(work_dirs) > 0:
        work_dir = work_dirs[-1]
    else:
        work_dir = None

    if work_dir and delete_pipeline_work_dir:
        try:
            shutil.rmtree(work_dir, ignore_errors=True)
            logging.info(json.dumps({
                "event_type": "analysis_work_dir_deleted",
                "sequencing_run_id": sequencing_run_id,
                "analysis_work_dir_path": work_dir
            }))
        except OSError as e:
            logging.error(json.dumps({
                "event_type": "delete_analysis_work_dir_failed",
                "sequencing_run_id": analysis_run_id,
                "analysis_work_dir_path": analysis_work_dir
            }))
    else:
        if not work_dir or not os.path.exists(work_dir):
            logging.warning(json.dumps({
                "event_type": "analysis_work_dir_not_found",
                "sequencing_run_id": sequencing_run_id,
                "analysis_work_dir_glob": work_dir_glob
            }))
        elif not delete_pipeline_work_dir:
            logging.info(json.dumps({
                "event_type": "skipped_deletion_of_analysis_work_dir",
                "sequencing_run_id": sequencing_run_id,
                "analysis_work_dir_path": work_dir
            }))

    if pipeline_name == 'BCCDC-PHL/basic-sequence-qc':
        return post_analysis_basic_sequence_qc(config, pipeline, run, analysis_mode)
    elif pipeline_name == 'BCCDC-PHL/basic-nanopore-qc':
        return post_analysis_basic_nanopore_qc(config, pipeline, run, analysis_mode)
    elif pipeline_name == 'BCCDC-PHL/taxon-abundance':
        return post_analysis_taxon_abundance(config, pipeline, run, analysis_mode)
    elif pipeline_name == 'BCCDC-PHL/routine-assembly':
        return post_analysis_routine_assembly(config, pipeline, run, analysis_mode)
    elif pipeline_name == 'BCCDC-PHL/mlst-nf':
        return post_analysis_mlst_nf(config, pipeline, run, analysis_mode)
    elif pipeline_name == 'BCCDC-PHL/plasmid-screen':
        return post_analysis_plasmid_screen(config, pipeline, run, analysis_mode)
    else:
        logging.warning(json.dumps({
            "event_type": "post_analysis_not_implemented",
            "sequencing_run_id": sequencing_run_id,
            "pipeline_name": pipeline_name
        }))
        return None
