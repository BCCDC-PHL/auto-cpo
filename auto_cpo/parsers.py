import csv
import json
import logging

from pathlib import Path


def parse_generic_csv(csv_path: Path, delimiter=',', fieldname_translation={}, int_fields=[], float_fields=[]):
    """
    Parse a csv file to a list of dicts. Optionally cast int and float fields to those types.
    Optionally translate fieldnames via a lookup table. Translation is done after casting, so
    fields listed in the `int_fields` and `float_fields` args should reflect the fieldnames
    in the original file (pre-translation, if applicable).

    `fieldname_translation` is a dict from original fieldname to translated fieldname.
    """
    parsed_rows = []
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        for row in reader:
            parsed_row = {}
            for field, value in row.items():
                if field in int_fields:
                    try:
                        parsed_row[field] = int(value)
                    except ValueError as e:
                        parsed_row[field] = None
                elif field in float_fields:
                    try:
                        parsed_row[field] = float(value)
                    except ValueError as e:
                        parsed_row[field] = None
                else:
                    parsed_row[field] = value

            for original_fieldname, translated_fieldname in fieldname_translation.items():
                if original_fieldname in parsed_row:
                    value = parsed_row.pop(original_fieldname)
                    parsed_row[translated_fieldname] = value

            parsed_rows.append(parsed_row)

    return parsed_rows


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


def parse_auto_cpo_sample_qc_summary(summary_path: Path):
    """
    """
    parsed_summary_by_library = {}
    int_fields = [
        'total_bases_input',
        'filtered_bases_input',
    ]
    float_fields = [
        'pre_alignment_estimated_depth_coverage_per_mb',
    ]
    with open(summary_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            library_id = row['library_id']
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
            parsed_summary_by_library[library_id] = row

    return parsed_summary_by_library


def parse_bracken_abundances_top_5(abundances_path: Path):
    """
    """
    int_fields = [
        'abundance_1_num_assigned_reads',
        'abundance_2_num_assigned_reads',
        'abundance_3_num_assigned_reads',
        'abundance_4_num_assigned_reads',
        'abundance_4_num_assigned_reads',
        'abundance_5_num_assigned_reads',
        'unclassified_num_assigned_reads',
    ]
    float_fields = [
        'abundance_1_fraction_total_reads',
        'abundance_2_fraction_total_reads',
        'abundance_3_fraction_total_reads',
        'abundance_4_fraction_total_reads',
        'abundance_5_fraction_total_reads',
        'unclassified_fraction_total_reads',
        
    ]
    top_5_abundances_by_library = {}
    with open(abundances_path, 'r') as f:
        reader = csv.DictReader(f)
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

            top_5_abundances_by_library[library_id] = row

    return top_5_abundances_by_library


def parse_mlst_sequence_type(mlst_sequence_type_path):
    """
    """
    mlst_sequence_type = parse_generic_csv(mlst_sequence_type_path)

    return mlst_sequence_type


def parse_abricate(abricate_report_path):
    """
    """
    abricate_records = []
    int_fields = [
        'start',
        'end',
    ]
    float_fields = [
        'percent_coverage',
        'percent_identity',
    ]
    with open(abricate_report_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            abricate_record = {}
            for k, v in row.items():
                fieldname = k.lower().replace('#', '').replace('%', 'percent_')
                abricate_record[fieldname] = v
            for field in int_fields:
                try:
                    abricate_record[field] = int(abricate_record[field])
                except ValueError as e:
                    abricate_record[field] = None
            for field in float_fields:
                try:
                    abricate_record[field] = float(abricate_record[field])
                except ValueError as e:
                    abricate_record[field] = None
            abricate_records.append(abricate_record)

    return abricate_records


def parse_resistance_gene_report(resistance_gene_report_path: Path):
    """
    """
    int_fields = [
        'resistance_gene_contig_size',
        'resistance_gene_contig_position_start',
        'resistance_gene_contig_position_end',
        'num_contigs_in_plasmid_reconstruction',
        'plasmid_reconstruction_size',
        'depth_coverage_threshold',
        'num_snps_vs_ref_plasmid'
    ]
    float_fields = [
        'percent_resistance_gene_coverage',
        'percent_resistance_gene_identity',
        'mash_neighbor_distance',
        'percent_ref_plasmid_coverage_above_depth_threshold',
    ]

    parsed_resistance_gene_report = parse_generic_csv(
        resistance_gene_report_path,
        delimiter='\t',
        int_fields=int_fields,
        float_fields=float_fields,
    )

    # We have an issue in the plasmid-screen pipeline that 
    # causes duplicate records in the resistance gene report.
    # This occurs when multiple carbapenemase genes are found
    # On the same plasmid.
    
    deduplicated_resistance_gene_report = []
    for record in parsed_resistance_gene_report:
        if record not in deduplicated_resistance_gene_report:
            deduplicated_resistance_gene_report.append(record)

    return deduplicated_resistance_gene_report


def parse_quast(quast_report_path: Path):
    """
    """
    int_fields = [
        'total_length',
        'num_contigs',
        'largest_contig',
        'assembly_N50',
        'assembly_N75',
        'assembly_L50',
        'assembly_L75',
        'num_contigs_gt_0_bp',
        'num_contigs_gt_1000_bp',
        'num_contigs_gt_5000_bp',
        'num_contigs_gt_10000_bp',
        'num_contigs_gt_25000_bp',
        'num_contigs_gt_50000_bp',
        'total_length_gt_0_bp',
        'total_length_gt_1000_bp',
        'total_length_gt_5000_bp',
        'total_length_gt_10000_bp',
        'total_length_gt_25000_bp',
        'total_length_gt_50000_bp',
    ]
    float_fields = [
        'num_N_per_100_kb'
    ]
    quast_report = parse_generic_csv(quast_report_path, int_fields=int_fields, float_fields=float_fields)

    return quast_report
