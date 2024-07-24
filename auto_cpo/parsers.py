import csv
import json
import logging

from pathlib import Path


def parse_basic_qc_stats(basic_qc_stats_csv_path: Path) -> dict:
    """
    Parse the basic sequence qc stats file

    :param basic_qc_stats_file: The input basic sequence qc stats file (csv)
    :type basic_qc_stats_file: str
    :return: The parsed basic sequence qc stats
    :rtype: dict
    """
    basic_qc_stats_by_library = {}
    int_fields = [
        'total_bases_before_filtering',
        'total_bases_after_filtering',
    ]
    float_fields = [
        'q30_rate_before_filtering',
        'q30_rate_after_filtering',
    ]
    with open(basic_qc_stats_csv_path, 'r') as f:
        reader = csv.DictReader(f, dialect='unix')
        for row in reader:
            library_id = row['sample_id']
            for field in int_fields:
                try:
                    row[field] = int(row[field])
                except ValueError as e:
                    row[field] = None
            for field in float_fields:
                try:
                    row[field] = float(row[field])
                except ValueError as e:
                    row[field] = None
            estimated_depth_coverage_per_megabase = None
            try:
                filtered_bases_input = row['total_bases_after_filtering']
                estimated_depth_coverage_per_megabase = round(filtered_bases_input / 1_000_000, 2)
            except ZeroDivisionError as e:
                pass
            except ValueError as e:
                pass
            
            q30_percent_before_filtering = None
            q30_percent_after_filtering = None

            if row['q30_rate_before_filtering'] is not None:
                q30_percent_before_filtering = round(row['q30_rate_before_filtering'] * 100, 2)
            if row['q30_rate_after_filtering'] is not None:
                q30_percent_after_filtering = round(row['q30_rate_after_filtering'] * 100, 2)
            
            basic_qc_stats_by_library[library_id] = {
                "library_id": library_id,
                "total_bases_input": row['total_bases_before_filtering'],
                "filtered_bases_input": row['total_bases_after_filtering'],
                "pre_alignment_estimated_depth_coverage_per_mb": estimated_depth_coverage_per_megabase,
                "q30_percent_before_filtering": q30_percent_before_filtering,
                "q30_percent_after_filtering": q30_percent_after_filtering,
            }

            logging.debug(json.dumps({
                "event_type": "estimated_coverage_by_library",
                "library_id": library_id,
                "estimated_coverage_per_mb": estimated_depth_coverage_per_megabase
            }))

    return basic_qc_stats_by_library


def parse_basic_nanopore_qc_stats(basic_nanopore_qc_stats_csv_path: Path) -> dict:
    """
    Parse the basic nanopore sequence qc stats file

    :param basic_nanopore_qc_stats_csv: The input basic nanopore sequence qc stats file (csv)
    :type basic_nanopore_qc_stats_csv: Path
    :return: The parsed basic nanopore sequence qc stats
    :rtype: dict
    """
    basic_qc_stats_by_library = {}
    int_fields = [
        'reads',
        'bases',
        'n50',
        'longest',
        'shortest',
        'mean_length',
        'median_length',
    ]
    float_fields = [
        'mean_quality',
        'median_quality',
    ]
    with open(basic_nanopore_qc_stats_csv_path, 'r') as f:
        reader = csv.DictReader(f, dialect='unix')
        for row in reader:
            library_id = row['sample_id']
            for field in int_fields:
                try:
                    row[field] = int(row[field])
                except ValueError as e:
                    row[field] = None
            for field in float_fields:
                try:
                    row[field] = float(row[field])
                except ValueError as e:
                    row[field] = None

            basic_qc_stats_by_library[library_id] = {
                "library_id": library_id,
                "reads": row['reads'],
                "bases": row['bases'],
                "n50": row['n50'],
                "longest": row['longest'],
                "shortest": row['shortest'],
                "mean_length": row['mean_length'],
                "median_length": row['median_length'],
                "mean_quality": row['mean_quality'],
                "median_quality": row['median_quality'],
            }

    return basic_qc_stats_by_library
