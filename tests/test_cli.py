from __future__ import annotations

from pathlib import Path

from flexgeo2.cli.main import build_config, build_parser


def test_build_config_maps_cli_flags() -> None:
    parser = build_parser()
    args = parser.parse_args(
        [
            "ensemble.pdb",
            "--output-dir",
            "out",
            "--chain",
            "A",
            "--chain",
            "B",
            "--n-jobs",
            "2",
            "--max-models-in-plot",
            "5",
            "--hide-model-traces",
            "--reference-pdb",
            "reference.pdb",
            "--reference-pdb-model",
            "3",
            "--cluster-residues",
            "--cluster-min-size",
            "7",
            "--cluster-min-samples",
            "2",
            "--cluster-residue-range",
            "10-12",
            "--output-verbose",
        ]
    )

    config = build_config(args)

    assert config.pdb_file == Path("ensemble.pdb")
    assert config.output.output_dir == Path("out")
    assert config.output.verbose is True
    assert config.output.write_files is True
    assert config.chains == ["A", "B"]
    assert config.n_jobs == 2
    assert config.max_models_in_plot == 5
    assert config.hide_model_traces is True
    assert config.reference is not None
    assert config.reference.pdb_file == Path("reference.pdb")
    assert config.reference.pdb_model_id == "3"
    assert config.clustering.cluster_residues is True
    assert config.clustering.cluster_residue_ranges == ["10-12"]
    assert config.clustering.min_cluster_size == 7
    assert config.clustering.min_samples == 2
