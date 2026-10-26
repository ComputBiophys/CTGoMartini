"""Unit tests for analysis of OpenMM >= 8.6 built-in REMD output directories.

Uses synthetic output directories, so these tests run without OpenMM 8.6 and
without openmmtools (proving the directory loading path is decoupled).
"""

import json

import numpy as np
import pytest

from ctgomartini.analysis.remd_exchange_ratio import ReportExchangeRatio
from ctgomartini.analysis.remd_mbar import BAR_NM3_TO_KJMOL, FESAnalyzer, kB

TEMPERATURES = [300.0, 310.0]
N_FRAMES = 3


def _make_openmm_dir(tmp_path, extra_molecule: bool = False):
    """Create a synthetic 2-state, 2-replica built-in REMD output directory."""
    outdir = tmp_path / "output_remd"
    outdir.mkdir()
    metadata = {
        "backend": "openmm",
        "openmm_version": "8.6.1",
        "method": "exp",
        "molname": "gbp",
        "replica_count": 2,
        "states": [
            {
                "temperature": 300.0,
                "ctg_gbp_beta": 0.001,
                "ctg_gbp_C1": -50.0,
                "ctg_gbp_C2": 0.0,
            },
            {
                "temperature": 310.0,
                "ctg_gbp_beta": 0.002,
                "ctg_gbp_C1": -50.0,
                "ctg_gbp_C2": 0.0,
            },
        ],
        "dt": 0.002,
        "exc_freq": 5,
        "stepsPerIteration": 5,
        "reportInterval": 1,
        "n_basins": 2,
    }
    (outdir / "metadata.json").write_text(json.dumps(metadata))

    # energy.csv: one row per frame, [replica][state] flattened reduced energies
    energy = 100.0 + np.arange(N_FRAMES * 4, dtype=float).reshape(N_FRAMES, 4)
    np.savetxt(outdir / "energy.csv", energy, delimiter=",")

    # log.csv: replica 0 swaps to state 1 at iteration 3
    (outdir / "log.csv").write_text(
        "Iteration,Step,Replica_0_State,Replica_1_State\n"
        "1,5,0,1\n"
        "2,10,0,1\n"
        "3,15,1,0\n"
    )

    # cv_energies.csv: E1/E2 in kJ/mol per replica per iteration
    lines = ["Iteration,Replica,Molecule,E1,E2"]
    for it in range(1, N_FRAMES + 1):
        for replica in range(2):
            e1 = 1000.0 + 100 * replica + it
            e2 = 2000.0 + 100 * replica + it
            lines.append(f"{it},{replica},gbp,{e1},{e2}")
            if extra_molecule:
                lines.append(f"{it},{replica},gbp_1,5.0,7.0")
    (outdir / "cv_energies.csv").write_text("\n".join(lines) + "\n")

    # CV file: n_states + 1 columns, one row per frame
    cv_file = tmp_path / "cv.dat"
    np.savetxt(cv_file, np.arange(N_FRAMES * 3, dtype=float).reshape(N_FRAMES, 3))
    return outdir, cv_file


def _make_multi_state_dir(tmp_path, temperatures, n_frames, n_basins=2, npt=None):
    """Create a synthetic built-in REMD output directory with n_states replicas.

    Basin energies are deterministic: E1 = 1000 + 100*replica + iteration,
    E2 = 2000 + 100*replica + iteration (and E3 = 3000 + ... when n_basins=3).
    log.csv rotates the state assignment, state[frame][replica] =
    (replica + frame) % n_states, so replicas visit different temperatures.
    The CV file holds 100*replica + frame per column, making the sample order
    directly assertable.

    With npt given (dict with p_type/pressure_bar and optionally
    surface_tension_bar_nm), the metadata gets an "npt" block and
    volume.csv (V = 500 + 10*replica + frame) is written; membrane coupling
    additionally writes area.csv (A = 100 + 5*replica + 0.5*frame).
    """
    n_states = len(temperatures)
    outdir = tmp_path / "output_remd"
    outdir.mkdir()
    metadata = {
        "backend": "openmm",
        "openmm_version": "8.6.1",
        "method": "exp",
        "molname": "gbp",
        "replica_count": n_states,
        "states": [
            {
                "temperature": float(t),
                "ctg_gbp_beta": 0.001 * (i + 1),
                "ctg_gbp_C1": -50.0,
                "ctg_gbp_C2": 0.0,
            }
            for i, t in enumerate(temperatures)
        ],
        "dt": 0.002,
        "exc_freq": 5,
        "stepsPerIteration": 5,
        "reportInterval": 1,
        "n_basins": n_basins,
    }
    if npt is not None:
        metadata["npt"] = {
            "enabled": True,
            "p_type": npt["p_type"],
            "pressure_bar": npt["pressure_bar"],
            "surface_tension_bar_nm": npt.get("surface_tension_bar_nm"),
        }
    (outdir / "metadata.json").write_text(json.dumps(metadata))

    energy = 100.0 + np.arange(n_frames * n_states**2, dtype=float).reshape(
        n_frames, n_states**2
    )
    np.savetxt(outdir / "energy.csv", energy, delimiter=",")

    if npt is not None:
        volume = np.array(
            [
                [500.0 + 10 * replica + frame for replica in range(n_states)]
                for frame in range(n_frames)
            ]
        )
        np.savetxt(outdir / "volume.csv", volume, delimiter=",")
        if npt["p_type"] == "membrane":
            area = np.array(
                [
                    [100.0 + 5 * replica + 0.5 * frame for replica in range(n_states)]
                    for frame in range(n_frames)
                ]
            )
            np.savetxt(outdir / "area.csv", area, delimiter=",")

    log_lines = [
        "Iteration,Step," + ",".join(f"Replica_{i}_State" for i in range(n_states))
    ]
    for f in range(n_frames):
        states_str = ",".join(str((r + f) % n_states) for r in range(n_states))
        log_lines.append(f"{f + 1},{5 * (f + 1)},{states_str}")
    (outdir / "log.csv").write_text("\n".join(log_lines) + "\n")

    basin_header = ",".join(f"E{s}" for s in range(1, n_basins + 1))
    lines = [f"Iteration,Replica,Molecule,{basin_header}"]
    for it in range(1, n_frames + 1):
        for replica in range(n_states):
            values = ",".join(
                str(1000.0 * s + 100 * replica + it) for s in range(1, n_basins + 1)
            )
            lines.append(f"{it},{replica},gbp,{values}")
    (outdir / "cv_energies.csv").write_text("\n".join(lines) + "\n")

    cv = np.zeros((n_frames, n_states + 1))
    for f in range(n_frames):
        cv[f, 0] = f
        for r in range(n_states):
            cv[f, r + 1] = 100 * r + f
    cv_file = tmp_path / "cv.dat"
    np.savetxt(cv_file, cv)
    return outdir, cv_file


def _bypass_initialize_fes(analyzer):
    """Select all frames without pymbar (mirrors initialize_fes bookkeeping)."""
    n_frames = analyzer.replica_energies.shape[2]
    analyzer.subsampled_indices = list(range(n_frames))
    analyzer.concatenated_cv = np.concatenate(analyzer.cv_values_replica)
    analyzer.subsampled_unsampled_energies = analyzer.unsampled_energies
    if analyzer._basin_energies_kj is not None:
        analyzer._subsampled_basin_energies_kj = analyzer._basin_energies_kj
    if analyzer._work_energies_kj is not None:
        analyzer._subsampled_work_kj = analyzer._work_energies_kj


def _expected_basin_energies(n_replicas, n_frames):
    """Replica-major raw E1/E2 arrays matching the synthetic builders."""
    e1 = np.array(
        [1000.0 + 100 * r + f + 1 for r in range(n_replicas) for f in range(n_frames)]
    )
    e2 = np.array(
        [2000.0 + 100 * r + f + 1 for r in range(n_replicas) for f in range(n_frames)]
    )
    return e1, e2


class TestFESAnalyzerOpenMMDirectory:
    def test_shapes_and_temperatures(self, tmp_path):
        outdir, cv_file = _make_openmm_dir(tmp_path)
        analyzer = FESAnalyzer(outdir, cv_file)

        assert analyzer.n_states == 2
        assert analyzer.interval == 1
        assert analyzer.temperatures_k.tolist() == TEMPERATURES
        assert analyzer.replica_energies.shape == (2, 2, N_FRAMES)
        assert analyzer.unsampled_energies.shape == (2, 2, N_FRAMES)
        assert analyzer._basin_energies_kj.shape == (2, 2, N_FRAMES)
        assert analyzer.cv_values_replica.shape == (2, N_FRAMES)

    def test_reduced_energies_transposed(self, tmp_path):
        outdir, cv_file = _make_openmm_dir(tmp_path)
        analyzer = FESAnalyzer(outdir, cv_file)

        energy = np.loadtxt(outdir / "energy.csv", delimiter=",")
        for frame in range(N_FRAMES):
            for replica in range(2):
                for state in range(2):
                    assert analyzer.replica_energies[state, replica, frame] == pytest.approx(
                        energy[frame][replica * 2 + state]
                    )

    def test_unsampled_reduced_by_frame_temperature(self, tmp_path):
        outdir, cv_file = _make_openmm_dir(tmp_path)
        analyzer = FESAnalyzer(outdir, cv_file)

        # Layout is [replica][basin][frame]. Frames 1-2: replica 0 in state 0
        # (300 K), replica 1 in state 1 (310 K); frame 3: swapped. The
        # off-diagonal assertions pin the layout (E1/E2 differ per replica).
        for frame, (it, t0) in enumerate(zip((1, 2, 3), (300.0, 300.0, 310.0), strict=True)):
            t1 = 610.0 - t0  # replica 1 occupies the other state
            e1_r0, e2_r0 = 1000.0 + it, 2000.0 + it
            e1_r1, e2_r1 = 1100.0 + it, 2100.0 + it
            assert analyzer.unsampled_energies[0, 0, frame] == pytest.approx(
                e1_r0 / (kB * t0)
            )
            assert analyzer.unsampled_energies[0, 1, frame] == pytest.approx(
                e2_r0 / (kB * t0)
            )
            assert analyzer.unsampled_energies[1, 0, frame] == pytest.approx(
                e1_r1 / (kB * t1)
            )
            assert analyzer.unsampled_energies[1, 1, frame] == pytest.approx(
                e2_r1 / (kB * t1)
            )
            assert analyzer._basin_energies_kj[0, 0, frame] == pytest.approx(e1_r0)
            assert analyzer._basin_energies_kj[1, 1, frame] == pytest.approx(e2_r1)

    def test_user_interval_downsamples(self, tmp_path):
        outdir, cv_file = _make_openmm_dir(tmp_path)
        analyzer = FESAnalyzer(outdir, cv_file, interval=2)

        assert analyzer.interval == 2
        assert analyzer.replica_energies.shape == (2, 2, 2)
        energy = np.loadtxt(outdir / "energy.csv", delimiter=",")
        assert analyzer.replica_energies[1, 0, 1] == pytest.approx(energy[2][1])
        # The raw basin energies are subsampled with the same indices:
        # replica 0, basin 1, subsampled frame 1 is iteration 3.
        assert analyzer._basin_energies_kj.shape == (2, 2, 2)
        assert analyzer._basin_energies_kj[0, 0, 1] == pytest.approx(1003.0)

    def test_multiple_molecules_summed(self, tmp_path):
        outdir, cv_file = _make_openmm_dir(tmp_path, extra_molecule=True)
        analyzer = FESAnalyzer(outdir, cv_file)

        # E1 of replica 0 at iteration 1: (1000 + 1) + 5
        assert analyzer.unsampled_energies[0, 0, 0] == pytest.approx(
            (1001.0 + 5.0) / (kB * 300.0)
        )

    def test_missing_files_raise(self, tmp_path):
        outdir = tmp_path / "broken"
        outdir.mkdir()
        cv_file = tmp_path / "cv.dat"
        np.savetxt(cv_file, np.zeros((3, 3)))
        with pytest.raises(RuntimeError, match="Error loading simulation data"):
            FESAnalyzer(outdir, cv_file)

    def test_single_frame_cv_file(self, tmp_path):
        outdir, cv_file = _make_openmm_dir(tmp_path)
        np.savetxt(cv_file, np.array([[0.0, 1.0, 2.0]]))
        analyzer = FESAnalyzer(outdir, cv_file)
        assert analyzer.cv_values_replica.shape == (2, 1)
        assert analyzer.cv_values_replica[:, 0].tolist() == [1.0, 2.0]

    def test_more_than_two_basins_rejected(self, tmp_path):
        outdir, cv_file = _make_multi_state_dir(
            tmp_path, [300.0, 310.0], 3, n_basins=3
        )
        with pytest.raises(RuntimeError, match="two-basin"):
            FESAnalyzer(outdir, cv_file)


class TestOpenMMDirectoryMixingNumerics:
    def test_unsampled_layout_replica_basin_frame(self, tmp_path):
        temperatures = [300.0, 310.0, 320.0]
        n_frames = 4
        outdir, cv_file = _make_multi_state_dir(tmp_path, temperatures, n_frames)
        analyzer = FESAnalyzer(outdir, cv_file)

        assert analyzer.unsampled_energies.shape == (3, 2, n_frames)
        assert analyzer._basin_energies_kj.shape == (3, 2, n_frames)
        log = np.loadtxt(outdir / "log.csv", delimiter=",", skiprows=1)
        for frame in range(n_frames):
            it = frame + 1
            for replica in range(3):
                state = int(log[frame][2 + replica])
                kt = kB * temperatures[state]
                e1 = 1000.0 + 100 * replica + it
                e2 = 2000.0 + 100 * replica + it
                # [replica][basin][frame]: transposed layouts are caught
                # because E1/E2 differ between replicas.
                assert analyzer.unsampled_energies[replica, 0, frame] == pytest.approx(
                    e1 / kt
                )
                assert analyzer.unsampled_energies[replica, 1, frame] == pytest.approx(
                    e2 / kt
                )
                assert analyzer._basin_energies_kj[replica, 0, frame] == pytest.approx(
                    e1
                )
                assert analyzer._basin_energies_kj[replica, 1, frame] == pytest.approx(
                    e2
                )

    def test_exp_mixing_mixed_temperatures_matches_analytic(self, tmp_path):
        n_frames = 4
        outdir, cv_file = _make_multi_state_dir(
            tmp_path, [300.0, 310.0, 320.0], n_frames
        )
        analyzer = FESAnalyzer(outdir, cv_file)
        _bypass_initialize_fes(analyzer)

        # Sample order is replica-major and matches concatenated_cv.
        assert analyzer.concatenated_cv.tolist() == [
            float(100 * r + f) for r in range(3) for f in range(n_frames)
        ]

        beta, c1, c2, t_analysis = 0.0015, -20.0, 5.0, 305.0
        e1, e2 = _expected_basin_energies(3, n_frames)
        expected = (
            -np.log(np.exp(-beta * (e1 + c1)) + np.exp(-beta * (e2 + c2))) / beta
        ) / (kB * t_analysis)
        u_n = analyzer._exp_mixing(beta, c1, c2, t_analysis)
        assert u_n.shape == (3 * n_frames,)
        np.testing.assert_allclose(u_n, expected, rtol=1e-12)

    def test_ham_mixing_mixed_temperatures_matches_analytic(self, tmp_path):
        n_frames = 4
        outdir, cv_file = _make_multi_state_dir(
            tmp_path, [300.0, 310.0, 320.0], n_frames
        )
        analyzer = FESAnalyzer(outdir, cv_file)
        _bypass_initialize_fes(analyzer)

        delta, c1, c2, t_analysis = 700.0, -20.0, 5.0, 305.0
        e1, e2 = _expected_basin_energies(3, n_frames)
        dv = c2 - c1
        expected = (
            (e1 + e2 + dv) / 2 - np.sqrt(((e1 - e2 - dv) / 2) ** 2 + delta**2)
        ) / (kB * t_analysis)
        u_n = analyzer._ham_mixing(delta, c1, c2, t_analysis)
        assert u_n.shape == (3 * n_frames,)
        np.testing.assert_allclose(u_n, expected, rtol=1e-12)

    def test_mixing_same_temperature_consistent(self, tmp_path):
        # H-REMD: equal state temperatures reproduce the analytic mixing
        # energies, and the reduced channel un-scaled by kT matches the raw
        # kJ/mol channel exactly.
        n_frames = 4
        outdir, cv_file = _make_multi_state_dir(
            tmp_path, [310.0, 310.0, 310.0], n_frames
        )
        analyzer = FESAnalyzer(outdir, cv_file)
        _bypass_initialize_fes(analyzer)

        beta, c1, c2, t_analysis = 0.0015, -20.0, 5.0, 310.0
        e1, e2 = _expected_basin_energies(3, n_frames)
        expected = (
            -np.log(np.exp(-beta * (e1 + c1)) + np.exp(-beta * (e2 + c2))) / beta
        ) / (kB * t_analysis)
        np.testing.assert_allclose(
            analyzer._exp_mixing(beta, c1, c2, t_analysis), expected, rtol=1e-12
        )

        e1_reduced = np.concatenate(analyzer.subsampled_unsampled_energies[:, 0, :])
        np.testing.assert_allclose(e1_reduced * (kB * t_analysis), e1, rtol=1e-12)


def _expected_work_kj(n_replicas, n_frames, pressure, gamma=None):
    """[replica][frame] NPT work (p*V - gamma*A) in kJ/mol for the synthetic dirs."""
    volume = np.array(
        [
            [500.0 + 10 * replica + frame for replica in range(n_replicas)]
            for frame in range(n_frames)
        ]
    )
    work = pressure * volume
    if gamma is not None:
        area = np.array(
            [
                [100.0 + 5 * replica + 0.5 * frame for replica in range(n_replicas)]
                for frame in range(n_frames)
            ]
        )
        work = work - gamma * area
    return (work * BAR_NM3_TO_KJMOL).T


class TestOpenMMDirectoryNPT:
    TEMPERATURES = [300.0, 310.0, 320.0]
    NPT_MEMBRANE = {
        "p_type": "membrane",
        "pressure_bar": 1.5,
        "surface_tension_bar_nm": 8.0,
    }
    NPT_ISOTROPIC = {"p_type": "isotropic", "pressure_bar": 2.0}

    def test_membrane_work_added_to_replica_energies(self, tmp_path):
        n_frames = 4
        outdir, cv_file = _make_multi_state_dir(
            tmp_path, self.TEMPERATURES, n_frames, npt=self.NPT_MEMBRANE
        )
        analyzer = FESAnalyzer(outdir, cv_file)

        work = _expected_work_kj(3, n_frames, 1.5, gamma=8.0)
        np.testing.assert_allclose(analyzer._work_energies_kj, work, rtol=1e-12)

        # u_s = (U_s + W) / (kB*T_s): the same per-replica work enters every
        # state, reduced by that state's temperature.
        energy = np.loadtxt(outdir / "energy.csv", delimiter=",")
        bare = energy.reshape(n_frames, 3, 3).transpose(2, 1, 0)
        expected = bare + work[None, :, :] / (
            kB * np.asarray(self.TEMPERATURES)[:, None, None]
        )
        np.testing.assert_allclose(analyzer.replica_energies, expected, rtol=1e-12)

    def test_isotropic_work_ignores_area(self, tmp_path):
        n_frames = 4
        outdir, cv_file = _make_multi_state_dir(
            tmp_path, self.TEMPERATURES, n_frames, npt=self.NPT_ISOTROPIC
        )
        assert not (outdir / "area.csv").exists()
        analyzer = FESAnalyzer(outdir, cv_file)
        work = _expected_work_kj(3, n_frames, 2.0)
        np.testing.assert_allclose(analyzer._work_energies_kj, work, rtol=1e-12)

    def test_interval_downsamples_work(self, tmp_path):
        n_frames = 4
        outdir, cv_file = _make_multi_state_dir(
            tmp_path, self.TEMPERATURES, n_frames, npt=self.NPT_ISOTROPIC
        )
        analyzer = FESAnalyzer(outdir, cv_file, interval=2)
        work = _expected_work_kj(3, n_frames, 2.0)[:, ::2]
        np.testing.assert_allclose(analyzer._work_energies_kj, work, rtol=1e-12)
        assert analyzer.replica_energies.shape == (3, 3, 2)

    def test_mixing_u_n_includes_work_term(self, tmp_path):
        n_frames = 4
        outdir, cv_file = _make_multi_state_dir(
            tmp_path, self.TEMPERATURES, n_frames, npt=self.NPT_MEMBRANE
        )
        analyzer = FESAnalyzer(outdir, cv_file)
        _bypass_initialize_fes(analyzer)

        captured = {}

        def fake_analyze(u_n, cv_n, temperature, **kwargs):
            captured["u_n"] = u_n
            return {}

        analyzer.analyze = fake_analyze
        beta, c1, c2, t_analysis = 0.0015, -20.0, 5.0, 305.0
        analyzer.analyze_one_mixing_parameters(
            {"beta": beta, "C1": c1, "C2": c2}, t_analysis, method="EXP"
        )

        e1, e2 = _expected_basin_energies(3, n_frames)
        mixing = (
            -np.log(np.exp(-beta * (e1 + c1)) + np.exp(-beta * (e2 + c2))) / beta
        ) / (kB * t_analysis)
        work = _expected_work_kj(3, n_frames, 1.5, gamma=8.0)
        # Replica-major concatenation, reduced at the analysis temperature.
        expected = mixing + np.concatenate(work) / (kB * t_analysis)
        np.testing.assert_allclose(captured["u_n"], expected, rtol=1e-12)

    def test_analyze_onestate_carries_work_exactly_once(self, tmp_path):
        n_frames = 4
        outdir, cv_file = _make_multi_state_dir(
            tmp_path, self.TEMPERATURES, n_frames, npt=self.NPT_ISOTROPIC
        )
        analyzer = FESAnalyzer(outdir, cv_file)
        _bypass_initialize_fes(analyzer)
        # u_kn is normally set by initialize_fes (needs pymbar); kln_to_kn
        # flattens [replica][frame] replica-major, so reshape does the same.
        analyzer.u_kn = analyzer.replica_energies.reshape(analyzer.n_states, -1)

        captured = {}

        def fake_analyze(u_n, cv_n, temperature, **kwargs):
            captured["u_n"] = u_n
            return {}

        analyzer.analyze = fake_analyze
        analyzer.analyze_onestate(selected_state=1)

        energy = np.loadtxt(outdir / "energy.csv", delimiter=",")
        bare_kn = energy.reshape(n_frames, 3, 3).transpose(2, 1, 0).reshape(3, -1)
        work_kn = (
            _expected_work_kj(3, n_frames, 2.0)[None, :, :]
            / (kB * np.asarray(self.TEMPERATURES)[:, None, None])
        ).reshape(3, -1)
        # The work term is present (exactly once) at the state's temperature.
        np.testing.assert_allclose(
            captured["u_n"], bare_kn[1] + work_kn[1], rtol=1e-12
        )

    def test_volume_frame_mismatch_rejected(self, tmp_path):
        n_frames = 4
        outdir, cv_file = _make_multi_state_dir(
            tmp_path, self.TEMPERATURES, n_frames, npt=self.NPT_ISOTROPIC
        )
        volume = np.loadtxt(outdir / "volume.csv", delimiter=",")
        np.savetxt(outdir / "volume.csv", volume[:-1], delimiter=",")
        with pytest.raises(RuntimeError, match="Frame count mismatch"):
            FESAnalyzer(outdir, cv_file)

    def test_nvt_directory_has_no_work(self, tmp_path):
        outdir, cv_file = _make_multi_state_dir(tmp_path, [300.0, 310.0], 3)
        analyzer = FESAnalyzer(outdir, cv_file)
        assert analyzer._work_energies_kj is None


class TestExchangeRatioDirectory:
    def test_ratios_from_cumulative_counters(self, tmp_path):
        outdir = tmp_path / "output_remd"
        outdir.mkdir()
        (outdir / "exchange.csv").write_text(
            "Iteration,Pair,Accepted,Proposed\n"
            "1,0-1,0,4\n"
            "1,0-2,1,4\n"
            "2,0-1,3,8\n"
            "2,0-2,3,8\n"
        )
        results = ReportExchangeRatio(outdir)
        assert results == {"0-1": pytest.approx(3 / 8), "0-2": pytest.approx(3 / 8)}

    def test_output_txt_format(self, tmp_path):
        outdir = tmp_path / "output_remd"
        outdir.mkdir()
        (outdir / "exchange.csv").write_text(
            "Iteration,Pair,Accepted,Proposed\n"
            "2,0-1,7,8\n"
        )
        output_txt = tmp_path / "ratios.dat"
        results = ReportExchangeRatio(outdir, output_txt=output_txt)
        assert results["0-1"] == pytest.approx(0.875)
        assert output_txt.read_text().splitlines() == [
            "# REMD Exchange Ratios",
            "# Replica Pair\tExchange Ratio",
            "0-1\t0.875000",
        ]

    def test_missing_exchange_csv(self, tmp_path):
        outdir = tmp_path / "output_remd"
        outdir.mkdir()
        with pytest.raises(FileNotFoundError, match="exchange.csv"):
            ReportExchangeRatio(outdir)
