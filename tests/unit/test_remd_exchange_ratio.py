"""
Unit tests for REMD exchange ratio analysis module.

Tests cover exchange ratio calculation and file output functionality.
"""

import numpy as np
import pytest
from pathlib import Path

from ctgomartini.analysis.remd_exchange_ratio import (
    ReportExchangeRatio,
)


class TestReportExchangeRatio:
    """Tests for ReportExchangeRatio function."""

    def test_report_exchange_ratio_from_fixture(self, remd_output_nc):
        """Test calculating exchange ratios from actual REMD fixture."""
        results = ReportExchangeRatio(remd_output_nc)
        
        # Should return dictionary
        assert isinstance(results, dict)
        assert len(results) > 0
        
        # Each key should be "i-j" format
        for key in results.keys():
            parts = key.split("-")
            assert len(parts) == 2
            assert parts[0].isdigit()
            assert parts[1].isdigit()
        
        # Exchange ratios should be non-negative
        for ratio in results.values():
            assert isinstance(ratio, float)
            assert ratio >= 0

    def test_report_exchange_ratio_with_output(self, remd_output_nc, tmp_path):
        """Test calculating exchange ratios and saving to file."""
        output_file = tmp_path / "exchange_ratios.dat"
        
        results = ReportExchangeRatio(remd_output_nc, output_txt=output_file)
        
        # Check file was created
        assert output_file.exists()
        
        # Check file content
        with open(output_file) as f:
            lines = f.readlines()
        
        # Should have header lines
        assert lines[0].startswith("# REMD Exchange Ratios")
        assert lines[1].startswith("# Replica Pair")
        
        # Data lines should match results
        data_lines = [l for l in lines[2:] if l.strip()]
        assert len(data_lines) == len(results)
        
        # Check each data line format
        for line in data_lines:
            parts = line.strip().split("\t")
            assert len(parts) == 2
            
            key = parts[0]
            ratio = float(parts[1])
            
            assert key in results
            assert pytest.approx(ratio, rel=0.001) == results[key]

    def test_exchange_ratio_values(self, remd_output_nc):
        """Test that exchange ratios are within reasonable range."""
        results = ReportExchangeRatio(remd_output_nc)
        
        for key, ratio in results.items():
            # Exchange ratios can be > 1.0 due to neighbor swap algorithm
            # but should be within reasonable bounds
            assert 0 <= ratio < 10, f"Exchange ratio {key}={ratio} out of bounds"

    def test_exchange_ratio_consistency(self, remd_output_nc):
        """Test that results are consistent across multiple calls."""
        results1 = ReportExchangeRatio(remd_output_nc)
        results2 = ReportExchangeRatio(remd_output_nc)
        
        assert results1.keys() == results2.keys()
        
        for key in results1:
            assert pytest.approx(results1[key], rel=0.001) == results2[key]
