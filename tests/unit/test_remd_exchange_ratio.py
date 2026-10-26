"""
Unit tests for REMD exchange ratio analysis module.

Tests cover exchange ratio calculation and file output functionality using
synthetic built-in REMD output directories (exchange.csv), so no OpenMM 8.6
or openmmtools is required.
"""

import pytest

from ctgomartini.analysis.remd_exchange_ratio import (
    ReportExchangeRatio,
)


def _make_remd_dir(tmp_path, pairs_rows=None):
    """Create a synthetic built-in REMD output directory with exchange.csv."""
    directory = tmp_path / "output_remd"
    directory.mkdir()
    if pairs_rows is None:
        pairs_rows = [
            "0,0-1,0,0",
            "0,0-2,0,0",
            "0,1-2,0,0",
            "1,0-1,1,2",
            "1,0-2,0,2",
            "1,1-2,2,2",
            "2,0-1,3,4",
            "2,0-2,1,4",
            "2,1-2,4,4",
        ]
    lines = ["Iteration,Pair,Accepted,Proposed"] + pairs_rows
    (directory / "exchange.csv").write_text("\n".join(lines) + "\n")
    return directory


class TestReportExchangeRatio:
    """Tests for ReportExchangeRatio function."""

    def test_report_exchange_ratio(self, tmp_path):
        """Test calculating exchange ratios from a synthetic REMD directory."""
        directory = _make_remd_dir(tmp_path)

        results = ReportExchangeRatio(directory)

        # Should return dictionary
        assert isinstance(results, dict)
        assert len(results) == 3

        # Each key should be "i-j" format
        for key in results:
            parts = key.split("-")
            assert len(parts) == 2
            assert parts[0].isdigit()
            assert parts[1].isdigit()

        # exchange.csv is cumulative: ratios come from the last row per pair
        assert results["0-1"] == pytest.approx(3 / 4)
        assert results["0-2"] == pytest.approx(1 / 4)
        assert results["1-2"] == pytest.approx(1.0)

    def test_report_exchange_ratio_with_output(self, tmp_path):
        """Test calculating exchange ratios and saving to file."""
        directory = _make_remd_dir(tmp_path)
        output_file = tmp_path / "exchange_ratios.dat"

        results = ReportExchangeRatio(directory, output_txt=output_file)

        # Check file was created
        assert output_file.exists()

        # Check file content
        with open(output_file) as f:
            lines = f.readlines()

        # Should have header lines
        assert lines[0].startswith("# REMD Exchange Ratios")
        assert lines[1].startswith("# Replica Pair")

        # Data lines should match results
        data_lines = [ln for ln in lines[2:] if ln.strip()]
        assert len(data_lines) == len(results)

        # Check each data line format
        for line in data_lines:
            parts = line.strip().split("\t")
            assert len(parts) == 2

            key = parts[0]
            ratio = float(parts[1])

            assert key in results
            assert pytest.approx(ratio, rel=0.001) == results[key]

    def test_exchange_ratio_values(self, tmp_path):
        """Test that exchange ratios are within reasonable range."""
        directory = _make_remd_dir(tmp_path)

        results = ReportExchangeRatio(directory)

        for key, ratio in results.items():
            assert 0 <= ratio < 10, f"Exchange ratio {key}={ratio} out of bounds"

    def test_exchange_ratio_consistency(self, tmp_path):
        """Test that results are consistent across multiple calls."""
        directory = _make_remd_dir(tmp_path)

        results1 = ReportExchangeRatio(directory)
        results2 = ReportExchangeRatio(directory)

        assert results1.keys() == results2.keys()

        for key in results1:
            assert pytest.approx(results1[key], rel=0.001) == results2[key]

    def test_rejects_netcdf_input(self, tmp_path):
        """Test that a non-directory (.nc) input raises ValueError."""
        nc_file = tmp_path / "output.nc"
        nc_file.write_text("not a directory")

        with pytest.raises(ValueError, match="no longer supported"):
            ReportExchangeRatio(nc_file)

    def test_missing_exchange_csv(self, tmp_path):
        """Test that a directory without exchange.csv raises FileNotFoundError."""
        directory = tmp_path / "output_remd"
        directory.mkdir()

        with pytest.raises(FileNotFoundError, match="exchange.csv"):
            ReportExchangeRatio(directory)
