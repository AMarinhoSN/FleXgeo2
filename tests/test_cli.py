from __future__ import annotations

from pathlib import Path

import pytest

from flexgeo2.cli.main import build_config, build_parser, main, parse_args
from flexgeo2.geometry import GeometryService


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
            "--dmax-outlier-fraction",
            "0.05",
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
            "--distance-matrices",
            "--overwrite",
            "--plot-residues",
            "10,A:11",
            "--plot-residues",
            "12-13",
            "--plot-format",
            "svg",
        ]
    )

    config = build_config(args)

    assert config.pdb_file == Path("ensemble.pdb")
    assert config.output.output_dir == Path("out")
    assert config.output.distance_matrices is True
    assert config.output.plot_residues == ["10,A:11", "12-13"]
    assert config.output.plot_format == "svg"
    assert config.output.overwrite is True
    assert config.chains == ["A", "B"]
    assert config.n_jobs == 2
    assert config.max_models_in_plot == 5
    assert config.hide_model_traces is True
    assert config.dmax_outlier_fraction == 0.05
    assert config.reference is not None
    assert config.reference.pdb_file == Path("reference.pdb")
    assert config.reference.pdb_model_id == "3"
    assert config.clustering.cluster_residues is True
    assert config.clustering.cluster_residue_ranges == ["10-12"]
    assert config.clustering.min_cluster_size == 7
    assert config.clustering.min_samples == 2


def test_clustering_sizes_default_to_the_library_defaults() -> None:
    # None lets ClusteringService scale the minimum cluster size with the ensemble.
    config = build_config(build_parser().parse_args(["ensemble.pdb", "--cluster-residues"]))

    assert config.clustering.min_cluster_size is None
    assert config.clustering.min_samples is None


@pytest.mark.parametrize("value", ["-0.1", "1", "nan"])
def test_parser_rejects_invalid_dmax_outlier_fraction(value: str) -> None:
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(["ensemble.pdb", "--dmax-outlier-fraction", value])


def test_parser_rejects_reference_pdb_model_without_reference_pdb(capsys) -> None:
    with pytest.raises(SystemExit) as excinfo:
        parse_args(build_parser(), ["ensemble.pdb", "--reference-pdb-model", "2"])

    assert excinfo.value.code == 2
    assert "--reference-pdb-model requires --reference-pdb" in capsys.readouterr().err


def test_parser_rejects_distance_matrices_without_reference(capsys) -> None:
    with pytest.raises(SystemExit) as excinfo:
        parse_args(build_parser(), ["ensemble.pdb", "--distance-matrices"])

    assert excinfo.value.code == 2
    assert "--distance-matrices requires --reference-model or --reference-pdb" in (
        capsys.readouterr().err
    )


def test_removed_output_verbose_points_to_its_replacements(capsys) -> None:
    with pytest.raises(SystemExit) as excinfo:
        parse_args(build_parser(), ["ensemble.pdb", "--output-verbose"])

    assert excinfo.value.code == 2
    err = capsys.readouterr().err
    assert "--output-verbose was removed" in err
    assert "--distance-matrices" in err
    assert "clusters/clusters.png" in err
    assert "--output-verbose" not in build_parser().format_help()


@pytest.mark.parametrize(
    "arguments",
    [
        ["--plot-format", "jpg"],
        ["--cluster-min-size", "1"],
        ["--cluster-min-samples", "0"],
        ["--max-models-in-plot", "-1"],
        ["--n-jobs", "0"],
        ["--n-jobs", "-2"],
    ],
)
def test_parser_rejects_invalid_integer_options(arguments: list[str]) -> None:
    with pytest.raises(SystemExit):
        build_parser().parse_args(["ensemble.pdb", *arguments])


def test_main_reports_missing_input_without_traceback(tmp_path: Path, capsys) -> None:
    exit_code = main([str(tmp_path / "missing.pdb"), "--output-dir", str(tmp_path / "out")])

    captured = capsys.readouterr()
    assert exit_code == 1
    assert captured.err.startswith("flexgeo2: error: Input PDB file not found")
    assert "Traceback" not in captured.err


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        (["--reference-model", "0"], "Reference model '0' was not found"),
        (["--chain", "Z"], "Chain(s) not found in the input: Z"),
        (["--cluster-residue-range", "8-12"], "8-12 on chain 'A' is incomplete"),
    ],
)
def test_main_validates_before_running_melodia(
    monkeypatch: pytest.MonkeyPatch,
    mini_ensemble_pdb: Path,
    tmp_path: Path,
    capsys,
    arguments: list[str],
    message: str,
) -> None:
    def fail_compute_geometry(*args, **kwargs):
        raise AssertionError("Melodia should not run when the input is invalid.")

    monkeypatch.setattr(GeometryService, "compute_geometry", fail_compute_geometry)

    exit_code = main([str(mini_ensemble_pdb), "--output-dir", str(tmp_path / "out"), *arguments])

    assert exit_code == 1
    assert message in capsys.readouterr().err
    assert not (tmp_path / "out").exists()


def test_main_refuses_to_reuse_output_dir_without_overwrite(
    monkeypatch: pytest.MonkeyPatch, mini_ensemble_pdb: Path, tmp_path: Path, capsys
) -> None:
    output_dir = tmp_path / "out"
    assert main([str(mini_ensemble_pdb), "--output-dir", str(output_dir)]) == 0
    first_run = (output_dir / "run.json").read_text()
    capsys.readouterr()

    def fail_compute_geometry(*args, **kwargs):
        raise AssertionError("Melodia should not run when the output folder is in use.")

    with monkeypatch.context() as patch:
        patch.setattr(GeometryService, "compute_geometry", fail_compute_geometry)
        exit_code = main([str(mini_ensemble_pdb), "--output-dir", str(output_dir)])

    captured = capsys.readouterr()
    assert exit_code == 1
    assert captured.err == (
        f"flexgeo2: error: output folder {output_dir.resolve()} is not empty. Use "
        "--overwrite to replace the outputs of an earlier run, or choose another "
        "--output-dir.\n"
    )
    assert (output_dir / "run.json").read_text() == first_run

    assert main([str(mini_ensemble_pdb), "--output-dir", str(output_dir), "--overwrite"]) == 0
    assert (output_dir / "run.json").read_text() != first_run
