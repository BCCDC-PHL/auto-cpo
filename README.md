# auto-cpo
Automated analysis of carbapenemase-producing organism (CPO) sequence data.

# Installation

# Usage
Start the tool as follows:

```bash
auto-cpo --config config.json
```

See the Configuration section of this document for details on preparing a configuration file.

More detailed logs can be produced by controlling the log level using the `--log-level` flag:

```bash
auto-cpo --config config.json --log-level debug
```

# Configuration
This tool takes a single config file, in JSON format, with the following structure:

```json
{
    "fastq_by_run_dir": "/path/to/fastq_symlinks_by_run",
    "analysis_output_dir": "/path/to/analysis_by_run",
    "analysis_work_dir": "/path/to/auto-cpo-work",
    "conda_cache_dir": "/path/to/.conda/envs",
    "notification_email_addresses": [
	"someone@example.org"
    ],
    "send_notification_emails": false,
    "scan_interval_seconds": 3600,
    "analyze_runs_in_reverse_order": true,
    "qc_filters": {
	"samples": {
	    "minimum_pre_alignment_estimated_depth_coverage_per_mb": 250,
	    "warning_pre_alignment_estimated_depth_coverage_per_mb": 300,
	    "minimum_q30_percent_before_filtering": 75,
	    "warning_q30_percent_before_filtering": 80
	}
    },
    "pipelines_by_analysis_mode": {
	"short": [
	    {
		"name": "BCCDC-PHL/basic-sequence-qc",
		"version": "v0.2.1",
		"dependencies": null,
		"parameters": {
		    "fastq_input": null,
		    "prefix": null,
		    "outdir": null
		}
	    },
	    {
		"name": "BCCDC-PHL/taxon-abundance",
		"version": "v0.1.7",
		"dependencies": [
		    {
			"pipeline_name": "BCCDC-PHL/basic-sequence-qc",
			"pipeline_version": "v0.2.1"
		    }
		],
		"parameters": {
		    "samplesheet_input": null,
		    "kraken_db": "/path/to/kraken2_db",
		    "bracken_db": "/path/to/bracken_db",
		    "read_length": null,
		    "outdir": null
		}
	    },
	    {
		"name": "BCCDC-PHL/routine-assembly",
		"version": "v0.4.5",
		"dependencies": [
		    {
			"pipeline_name": "BCCDC-PHL/basic-sequence-qc",
			"pipeline_version": "v0.2.1"
		    }
		],
		"parameters": {
		    "samplesheet_input": null,
		    "bakta": null,
		    "outdir": null
		}
	    },
	    {
		"name": "BCCDC-PHL/mlst-nf",
		"version": "v0.1.4",
		"dependencies": [
		    {
			"pipeline_name": "BCCDC-PHL/routine-assembly",
			"pipeline_version": "v0.4.5"
		    }
		],
		"parameters": {
		    "assembly_input": null,
		    "outdir": null
		}
	    },
	    {
		"name": "BCCDC-PHL/plasmid-screen",
		"version": "v0.2.3",
		"dependencies": [
		    {
			"pipeline_name": "BCCDC-PHL/routine-assembly",
			"pipeline_version": "v0.4.5"
		    }
		],
		"parameters": {
		    "samplesheet_input": null,
		    "pre_assembled": null,
		    "mob_db": "/path/to/mob-suite_db",
		    "outdir": null
		}
	    }
	],
	"hybrid": [
	    {
		"name": "BCCDC-PHL/basic-nanopore-qc",
		"version": "v0.1.0",
		"dependencies": null,
		"parameters": {
		    "fastq_input": null,
		    "prefix": null,
		    "outdir": null
		}
	    },
	    {
		"name": "BCCDC-PHL/plasmid-assembly",
		"version": "v0.1.1",
		"dependencies": [
		    {
			"pipeline_name": "BCCDC-PHL/basic-nanopore-qc",
			"pipeline_version": "v0.1.0"
		    }
		],
		"parameters": {
		    "samplesheet_input": null,
		    "hybrid": null,
		    "chromosome_length": 250000,
		    "db": "/path/to/plassembler_db",
		    "outdir": null
		}
	    },
	    {
		"name": "BCCDC-PHL/mlst-nf",
		"version": "v0.1.4",
		"dependencies": [
		    {
			"pipeline_name": "BCCDC-PHL/plasmid-assembly",
			"pipeline_version": "v0.1.1"
		    }
		],
		"parameters": {
		    "assembly_input": null,
		    "outdir": null
		}
	    },
	    {
		"name": "BCCDC-PHL/plasmid-screen",
		"version": "v0.2.3",
		"dependencies": [
		    {
			"pipeline_name": "BCCDC-PHL/plasmid-assembly",
			"pipeline_version": "v0.1.1"
		    }
		],
		"parameters": {
                    "samplesheet_input": null,
		    "pre_assembled": null,
		    "mob_db": "/path/to/mob-suite_db",
		    "outdir": null
		}
	    }
	]
    }
}
```

# Logging
This tool outputs [structured logs](https://www.honeycomb.io/blog/structured-logging-and-your-team/) in [JSON Lines](https://jsonlines.org/) format:

Every log line should include the fields:

- `timestamp`
- `level`
- `module`
- `function_name`
- `line_num`
- `message`

...and the contents of the `message` key will be a JSON object that includes at `event_type`. The remaining keys inside the `message` will vary by event type.

```json
{"timestamp": "2022-09-22T11:32:52.287", "level": "INFO", "module": "core", "function_name": "scan", "line_num": 56, "message": {"event_type": "scan_start"}}
```
