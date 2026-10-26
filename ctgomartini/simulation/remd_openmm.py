"""Built-in OpenMM REMD backend (requires OpenMM >= 8.6).

Provides :class:`BuiltinREMDRunner`, a replica exchange runner built on
``openmm.app.replicaexchangesampler.ReplicaExchangeSampler``. The whole system
is simulated in a single :class:`openmm.Context`; replicas differ in their
temperatures and in the values of the MBP global parameters
(``ctg_<mol>_beta``/``ctg_<mol>_delta`` and ``ctg_<mol>_C1``/``ctg_<mol>_C2``),
which the sampler switches through ``Context.setParameter``.

The OpenMM 8.6 sampler/reporter classes are imported lazily inside functions so
this module stays importable under OpenMM < 8.6.

With ``remd_devices`` set in the input file, the runner switches to the
parallel multi-Context driver in :mod:`ctgomartini.simulation.remd_parallel`
(one worker Context per replica, integrated concurrently in threads); output
files and formats are identical to the serial driver.

Output layout (``remd_output`` is a directory for this backend):

- ``log.csv``, ``energy.csv``, ``checkpoint_<i>.xml`` and
  optional ``replica_<i>.xtc`` are written by the built-in
  ``ReplicaExchangeReporter``. With ``pcouple = yes`` the reporter also
  writes ``volume.csv`` (per-replica box volumes in nm^3, one row per
  report interval, no header); ``p_type = membrane`` additionally
  produces ``area.csv`` (per-replica xy areas in nm^2, same layout,
  written by :class:`AreaReporter`).
- ``cv_energies.csv``: long table ``Iteration,Replica,Molecule,E1,E2[,E3...]``
  with the per-basin energies (kJ/mol) of every multiple-basin CustomCVForce
  of the target molecule (``replica_molname``), collected from each force's
  collective variables (``state<s>_force<j>``).
- ``exchange.csv``: long table ``Iteration,Pair,Accepted,Proposed`` with
  cumulative per-pair exchange counters.
- ``metadata.json``: states, method, molecule name, timestep and reporting
  intervals needed by the analysis tools, plus an ``npt`` block recording
  the pressure coupling settings.

NPT (``pcouple = yes``) is supported with the correct exchange criterion.
The reduced potential of state s is u_s = beta_s * (U + p*V - gamma*A)
(the gamma*A term is subtracted, matching openmmtools'
``states.py:_compute_reduced_potential``), so for mixed-temperature REMD the
acceptance exponent gains the correction term
``(1/kT_si - 1/kT_sj) * (W_i - W_j)`` with ``W = p*V - gamma*A`` in kJ/mol
(:func:`exchange_exponent`). At equal state temperatures the correction is
identically zero and the criterion reduces to the built-in one. The MBAR
analysis (:mod:`ctgomartini.analysis.remd_mbar`) adds the same work term back
from ``volume.csv``/``area.csv``.
"""

from __future__ import annotations

import csv
import datetime
import json
import math
import os
import random
import re
from pathlib import Path
from typing import Any

import openmm as mm
import openmm.unit as u
from openmm.app import Simulation

from ctgomartini.simulation.base import (
    SimulationRunner,
    load_structure,
    report_time,
)
from ctgomartini.topology.interactions.mixing import _sanitize_molname

_BASIN_CV_PATTERN = re.compile(r"^state(\d+)_force(\d+)$")


# ---------------------------------------------------------------------------
# Pure helpers (no OpenMM 8.6 dependency)
# ---------------------------------------------------------------------------


def build_replica_states(
    molname: str,
    method: str,
    temperatures: list[float],
    couplings: list[float],
    c1_list: list[float],
    c2_list: list[float],
) -> list[dict[str, Any]]:
    """Build the per-state dicts for the OpenMM 8.6 ReplicaExchangeSampler.

    Each state sets the simulation temperature and the MBP global parameters of
    the target molecule in its per-molecule namespace
    (``ctg_<mol>_beta`` for EXP, ``ctg_<mol>_delta`` for HAM, plus
    ``ctg_<mol>_C1``/``ctg_<mol>_C2``).

    Args:
        molname: Target molecule name (sanitized for parameter names).
        method: Mixing method, 'exp' or 'ham'.
        temperatures: Per-state temperatures in Kelvin.
        couplings: Per-state coupling constants (beta for EXP, delta for HAM).
        c1_list: Per-state C1 basin energy offsets.
        c2_list: Per-state C2 basin energy offsets.

    Returns:
        List of state dicts accepted by ReplicaExchangeSampler.

    Raises:
        ValueError: If the method is unknown or the list lengths differ.
    """
    method = method.lower()
    if method not in ("exp", "ham"):
        raise ValueError(
            f"Unsupported replica_method '{method}'. Must be 'exp' or 'ham'."
        )
    n_states = len(temperatures)
    for name, values in (("coupling", couplings), ("C1", c1_list), ("C2", c2_list)):
        if len(values) != n_states:
            raise ValueError(
                f"replica {name} has {len(values)} values but there are "
                f"{n_states} replica temperatures."
            )
    namespace = f"ctg_{_sanitize_molname(molname)}"
    coupling_key = f"{namespace}_beta" if method == "exp" else f"{namespace}_delta"
    states: list[dict[str, Any]] = []
    for i in range(n_states):
        states.append(
            {
                "temperature": float(temperatures[i]),
                coupling_key: float(couplings[i]),
                f"{namespace}_C1": float(c1_list[i]),
                f"{namespace}_C2": float(c2_list[i]),
            }
        )
    return states


def find_mbp_forces(system: mm.System, molname: str) -> list[mm.CustomCVForce]:
    """Find the multiple-basin CustomCVForces of a molecule in a System.

    Args:
        system: The OpenMM System to inspect.
        molname: Molecule name whose parameter namespace identifies the forces.

    Returns:
        List of matching CustomCVForce objects (one per molecule copy).
    """
    namespace = f"ctg_{_sanitize_molname(molname)}"
    markers = (f"{namespace}_beta", f"{namespace}_delta")
    forces: list[mm.CustomCVForce] = []
    for force in system.getForces():
        if not isinstance(force, mm.CustomCVForce):
            continue
        names = {
            force.getGlobalParameterName(k) for k in range(force.getNumGlobalParameters())
        }
        if any(marker in names for marker in markers):
            forces.append(force)
    return forces


def basin_cv_groups(force: mm.CustomCVForce) -> list[list[int]]:
    """Group a MBP force's collective variable indices by basin state.

    CVs are named ``state<s>_force<j>``; the energy of basin s is the sum of
    its CVs. This is independent of the beta/delta/C1/C2 global parameters.

    Args:
        force: The multiple-basin CustomCVForce.

    Returns:
        List of CV index lists, one per basin, ordered by state index.

    Raises:
        ValueError: If no basin CVs are found.
    """
    groups: dict[int, list[int]] = {}
    for k in range(force.getNumCollectiveVariables()):
        match = _BASIN_CV_PATTERN.match(force.getCollectiveVariableName(k))
        if match:
            groups.setdefault(int(match.group(1)), []).append(k)
    if not groups:
        raise ValueError(
            "No basin collective variables (state<s>_force<j>) found in the force."
        )
    return [groups[s] for s in sorted(groups)]


def validate_state_parameters(
    context: mm.Context,
    states: list[dict[str, Any]],
    molname: str,
    method: str,
) -> None:
    """Check that every state parameter exists as a Context global parameter.

    Args:
        context: The OpenMM Context of the single combined System.
        states: State dicts produced by :func:`build_replica_states`.
        molname: Target molecule name, used for error reporting.
        method: Mixing method, used for error reporting.

    Raises:
        ValueError: If a state parameter is missing from the Context.
    """
    available = set(context.getParameters())
    missing = sorted(
        {key for state in states for key in state if key != "temperature"} - available
    )
    if missing:
        raise ValueError(
            f"REMD state parameter(s) {missing} not found in the Context. Check that "
            f"replica_molname '{molname}' and replica_method '{method}' match the "
            f"multiple-basin definition in the topology: EXP expects "
            f"'ctg_<mol>_beta', HAM expects 'ctg_<mol>_delta', and both expect "
            f"'ctg_<mol>_C1'/'ctg_<mol>_C2' (non-alphanumeric characters in the "
            f"molecule name become '_')."
        )


def read_exchange_statistics(csv_path: str | Path) -> dict[tuple[int, int], list[int]]:
    """Read cumulative exchange counters from an exchange.csv file.

    The file is cumulative, so the last row of each pair holds the totals.

    Args:
        csv_path: Path to exchange.csv.

    Returns:
        Mapping of pair ``(i, j)`` with ``i < j`` to ``[accepted, proposed]``.
    """
    stats: dict[tuple[int, int], list[int]] = {}
    with open(csv_path) as f:
        for row in csv.reader(f):
            if not row or row[0] == "Iteration":
                continue
            i, j = row[1].split("-")
            stats[(int(i), int(j))] = [int(row[2]), int(row[3])]
    return stats


# 1 bar * nm^3 per particle, in kJ/mol (N_A * 1e5 Pa * 1e-27 m^3 / 1000).
BAR_NM3_TO_KJMOL = 0.0602214076


def npt_work_kjmol(
    volume_nm3: float,
    area_nm2: float,
    pressure_bar: float,
    surface_tension_bar_nm: float = 0.0,
) -> float:
    """NPT work term ``W = p*V - gamma*A`` of one replica conformation.

    The gamma*A term is subtracted, matching the openmmtools reduced
    potential ``u = beta * (U + p*V - gamma*A)``.

    Args:
        volume_nm3: Periodic box volume in nm^3.
        area_nm2: Box cross-sectional xy area in nm^2 (a.x * b.y).
        pressure_bar: Reference pressure in bar.
        surface_tension_bar_nm: Membrane surface tension in bar*nm
            (0 for isotropic pressure coupling).

    Returns:
        Work term in kJ/mol.
    """
    return (
        pressure_bar * volume_nm3 - surface_tension_bar_nm * area_nm2
    ) * BAR_NM3_TO_KJMOL


def _as_kjmol(value: Any) -> float:
    """Strip openmm units, assuming the value is in or convertible to kJ/mol."""
    if u.is_quantity(value):
        return float(value.value_in_unit(u.kilojoules_per_mole))
    return float(value)


def exchange_exponent(
    e_i_si: Any,
    e_i_sj: Any,
    e_j_si: Any,
    e_j_sj: Any,
    kt_si: Any,
    kt_sj: Any,
    work_i: float = 0.0,
    work_j: float = 0.0,
) -> float:
    """Exchange acceptance exponent for one replica pair (mixed temperatures).

    With u_s = beta_s * (U + W_s), W = p*V - gamma*A being equal across
    states, the exponent of min(1, exp(exponent)) is

        (U_i(si) - U_j(si))/kT_si + (U_j(sj) - U_i(sj))/kT_sj
        + (1/kT_si - 1/kT_sj) * (W_i - W_j)

    where replica r occupies state s_r. All energies may be plain floats
    (kJ/mol) or openmm Quantities; the return value is a dimensionless float.

    Args:
        e_i_si: Energy of replica i evaluated in state si.
        e_i_sj: Energy of replica i evaluated in state sj.
        e_j_si: Energy of replica j evaluated in state si.
        e_j_sj: Energy of replica j evaluated in state sj.
        kt_si: kT of state si (kJ/mol or Quantity).
        kt_sj: kT of state sj (kJ/mol or Quantity).
        work_i: NPT work term of replica i's conformation (kJ/mol).
        work_j: NPT work term of replica j's conformation (kJ/mol).

    Returns:
        Dimensionless acceptance exponent.
    """
    kt_si = _as_kjmol(kt_si)
    kt_sj = _as_kjmol(kt_sj)
    exponent = (_as_kjmol(e_i_si) - _as_kjmol(e_j_si)) / kt_si + (
        _as_kjmol(e_j_sj) - _as_kjmol(e_i_sj)
    ) / kt_sj
    delta_work = work_i - work_j
    if delta_work != 0.0:
        exponent += delta_work * (1.0 / kt_si - 1.0 / kt_sj)
    return exponent


def box_volume_nm3(state: mm.State) -> float:
    """Periodic box volume of a State, in nm^3."""
    return float(state.getPeriodicBoxVolume().value_in_unit(u.nanometer**3))


def box_area_xy_nm2(state: mm.State) -> float:
    """Box cross-sectional xy area of a State, in nm^2 (a.x * b.y)."""
    a, b, _ = state.getPeriodicBoxVectors()
    return float(
        a[0].value_in_unit(u.nanometer) * b[1].value_in_unit(u.nanometer)
    )


# ---------------------------------------------------------------------------
# Reporters
# ---------------------------------------------------------------------------


class BasinEnergyReporter:
    """Record per-basin energies of every MBP force for each replica.

    For each replica the reporter loads its conformation into the Context and
    evaluates the collective variables of every multiple-basin CustomCVForce,
    whose per-state sums are the basin energies E1, E2, ... in kJ/mol. The
    sampler always resets the Context before further use, so temporarily
    changing it here is safe.
    """

    def __init__(
        self,
        file_path: str | Path,
        report_interval: int,
        labeled_forces: list[tuple[str, mm.CustomCVForce]],
        resume: bool = False,
    ) -> None:
        """Create the reporter.

        Args:
            file_path: CSV output path (cv_energies.csv).
            report_interval: Reporting interval in iterations.
            labeled_forces: ``(label, force)`` pairs, one per MBP force.
            resume: Append without rewriting the header.
        """
        self.file_path = str(file_path)
        self.report_interval = int(report_interval)
        self.entries = [
            (label, force, basin_cv_groups(force)) for label, force in labeled_forces
        ]
        self.n_basins_max = max(len(groups) for _, _, groups in self.entries)
        self._file = open(self.file_path, "a" if resume else "w")  # noqa: SIM115
        if not resume:
            header = ["Iteration", "Replica", "Molecule"] + [
                f"E{s}" for s in range(1, self.n_basins_max + 1)
            ]
            print(",".join(header), file=self._file)
            self._file.flush()

    def __call__(self, sampler: Any) -> None:
        """Write one row per replica and MBP force for this iteration."""
        if sampler.currentIteration % self.report_interval != 0:
            return
        context = sampler.simulation.context
        for replica_index, conformation in enumerate(sampler.replicaConformation):
            context.setState(conformation)
            for label, force, groups in self.entries:
                cv_values = force.getCollectiveVariableValues(context)
                basin_energies = [
                    sum(cv_values[k] for k in indices) for indices in groups
                ]
                basin_energies += [""] * (self.n_basins_max - len(basin_energies))
                fields = [sampler.currentIteration, replica_index, label, *basin_energies]
                print(",".join(str(x) for x in fields), file=self._file)
        self._file.flush()

    def close(self) -> None:
        """Close the output file."""
        if self._file is not None:
            self._file.close()
            self._file = None

    def __del__(self) -> None:
        self.close()


class ExchangeStatisticsReporter:
    """Record cumulative per-pair exchange counters to exchange.csv.

    The built-in sampler calls reporters before that iteration's exchange
    attempts, so ``__call__`` writes the counters as of the *previous*
    iteration (labeled ``currentIteration - 1``). The runner must call
    :meth:`finalize` after ``simulate()`` returns to write the last
    iteration's totals. Every row therefore reads "counters after completing
    iteration X", and the last row per pair always holds the exact totals
    (used to restore the counters when resuming).
    """

    def __init__(
        self, file_path: str | Path, report_interval: int, resume: bool = False
    ) -> None:
        """Create the reporter.

        Args:
            file_path: CSV output path (exchange.csv).
            report_interval: Reporting interval in iterations.
            resume: Append without rewriting the header.
        """
        self.file_path = str(file_path)
        self.report_interval = int(report_interval)
        self._file = open(self.file_path, "a" if resume else "w")  # noqa: SIM115
        if not resume:
            print("Iteration,Pair,Accepted,Proposed", file=self._file)
            self._file.flush()

    def _write_rows(self, sampler: Any, iteration_label: int) -> None:
        n_states = len(sampler.states)
        for i in range(n_states):
            for j in range(i + 1, n_states):
                fields = [
                    iteration_label,
                    f"{i}-{j}",
                    sampler.n_accepted[i][j],
                    sampler.n_proposed[i][j],
                ]
                print(",".join(str(x) for x in fields), file=self._file)
        self._file.flush()

    def __call__(self, sampler: Any) -> None:
        """Write the counters as of the previous iteration."""
        if sampler.currentIteration % self.report_interval != 0:
            return
        self._write_rows(sampler, sampler.currentIteration - 1)

    def finalize(self, sampler: Any) -> None:
        """Write the totals after the final completed iteration."""
        self._write_rows(sampler, sampler.currentIteration)

    def close(self) -> None:
        """Close the output file."""
        if self._file is not None:
            self._file.close()
            self._file = None

    def __del__(self) -> None:
        self.close()


class AreaReporter:
    """Record per-replica box xy areas to area.csv (membrane NPT only).

    One row per report interval holding one area per replica in nm^2, no
    header — the same layout as the built-in reporter's volume.csv, so the
    two files stay row-aligned. Attached only for ``p_type = membrane``.
    """

    def __init__(
        self, file_path: str | Path, report_interval: int, resume: bool = False
    ) -> None:
        """Create the reporter.

        Args:
            file_path: CSV output path (area.csv).
            report_interval: Reporting interval in iterations.
            resume: Append to the existing file.
        """
        self.file_path = str(file_path)
        self.report_interval = int(report_interval)
        self._file = open(self.file_path, "a" if resume else "w")  # noqa: SIM115

    def __call__(self, sampler: Any) -> None:
        """Write one row of per-replica xy areas for this iteration."""
        if sampler.currentIteration % self.report_interval != 0:
            return
        areas = [box_area_xy_nm2(conf) for conf in sampler.replicaConformation]
        print(",".join(str(area) for area in areas), file=self._file)
        self._file.flush()

    def close(self) -> None:
        """Close the output file."""
        if self._file is not None:
            self._file.close()
            self._file = None

    def __del__(self) -> None:
        self.close()


# ---------------------------------------------------------------------------
# Sampler (OpenMM >= 8.6, imported lazily)
# ---------------------------------------------------------------------------

_COUNTING_SAMPLER_CLASS = None


def _require_openmm86() -> None:
    """Raise a descriptive ImportError if the OpenMM 8.6 sampler is missing."""
    try:
        import openmm.app.replicaexchangereporter  # noqa: F401
        import openmm.app.replicaexchangesampler  # noqa: F401
    except ImportError as e:
        raise ImportError(
            "The 'openmm' REMD backend requires OpenMM >= 8.6 "
            f"(found {mm.__version__}), which provides "
            "openmm.app.replicaexchangesampler. Upgrade OpenMM or use "
            "remd_backend = openmmtools."
        ) from e


def _get_counting_sampler_class() -> type:
    """Return a ReplicaExchangeSampler subclass with per-pair exchange counters.

    The exchange logic is identical to the built-in random-pair scheme; only
    the ``n_proposed``/``n_accepted`` bookkeeping (unordered pairs i<j), a
    restore hook for resuming, and the NPT work-term correction
    (:func:`exchange_exponent`) are added.
    """
    global _COUNTING_SAMPLER_CLASS
    if _COUNTING_SAMPLER_CLASS is not None:
        return _COUNTING_SAMPLER_CLASS

    _require_openmm86()
    from openmm.app.replicaexchangesampler import ReplicaExchangeSampler

    class CountingReplicaExchangeSampler(ReplicaExchangeSampler):
        """ReplicaExchangeSampler that records exchange statistics.

        Args:
            pressure_bar: Reference pressure (bar) for the NPT exchange
                correction; None disables it (NVT).
            surface_tension_bar_nm: Membrane surface tension (bar*nm) for
                the NPT exchange correction; 0 for isotropic coupling.
        """

        def __init__(
            self,
            *args: Any,
            pressure_bar: float | None = None,
            surface_tension_bar_nm: float = 0.0,
            **kwargs: Any,
        ) -> None:
            super().__init__(*args, **kwargs)
            n = len(self.states)
            self.n_proposed = [[0] * n for _ in range(n)]
            self.n_accepted = [[0] * n for _ in range(n)]
            self.pressure_bar = pressure_bar
            self.surface_tension_bar_nm = float(surface_tension_bar_nm)
            self._kT_list: list[Any] | None = None
            if "temperature" in self.states[0]:
                temperatures = []
                for state in self.states:
                    t = state["temperature"]
                    temperatures.append(t if u.is_quantity(t) else t * u.kelvin)
                self._kT_list = [u.MOLAR_GAS_CONSTANT_R * t for t in temperatures]

        def set_exchange_statistics(
            self, stats: dict[tuple[int, int], list[int]]
        ) -> None:
            """Restore cumulative counters read from a previous exchange.csv."""
            for (i, j), (accepted, proposed) in stats.items():
                self.n_accepted[i][j] = accepted
                self.n_proposed[i][j] = proposed

        def _npt_works_kjmol(self) -> list[float] | None:
            """Per-replica NPT work terms (p*V - gamma*A) in kJ/mol.

            Computed from the current replica conformations (which exchange
            attempts only reassign, never modify), or None for NVT runs.
            """
            if self.pressure_bar is None:
                return None
            return [
                npt_work_kjmol(
                    box_volume_nm3(conformation),
                    box_area_xy_nm2(conformation),
                    self.pressure_bar,
                    self.surface_tension_bar_nm,
                )
                for conformation in self.replicaConformation
            ]

        def exchangeReplicas(self) -> None:  # noqa: N802
            """Attempt exchanges between random replica pairs (built-in logic).

            Under NPT with mixed state temperatures the acceptance exponent
            includes the (1/kT_si - 1/kT_sj)*(W_i - W_j) work correction
            (:func:`exchange_exponent`); same-temperature states and NVT runs
            are exactly the built-in criterion.
            """
            n = len(self.states)
            if n < 2:
                return
            # The correction is identically zero at equal temperatures, so it
            # is only needed on the per-state-temperature branch.
            works = self._npt_works_kjmol() if self._kT_list is not None else None
            for _ in range(self.exchangesPerIteration):
                # Select a random pair of replicas.
                i = j = 0
                while i == j:
                    i = random.randrange(n)
                    j = random.randrange(n)
                si = self.replicaStateIndex[i]
                sj = self.replicaStateIndex[j]

                # Compute the probability of exchanging their states.
                if self._kT_list is None:
                    kt = (
                        u.MOLAR_GAS_CONSTANT_R
                        * self.simulation.integrator.getTemperature()
                    )
                    exponent = (
                        self.replicaStateEnergy[i][si]
                        - self.replicaStateEnergy[i][sj]
                        + self.replicaStateEnergy[j][sj]
                        - self.replicaStateEnergy[j][si]
                    ) / kt
                else:
                    work_i = works[i] if works is not None else 0.0
                    work_j = works[j] if works is not None else 0.0
                    exponent = exchange_exponent(
                        self.replicaStateEnergy[i][si],
                        self.replicaStateEnergy[i][sj],
                        self.replicaStateEnergy[j][si],
                        self.replicaStateEnergy[j][sj],
                        self._kT_list[si],
                        self._kT_list[sj],
                        work_i,
                        work_j,
                    )
                probability = min(1.0, math.exp(exponent))
                a, b = (i, j) if i < j else (j, i)
                self.n_proposed[a][b] += 1
                if random.random() <= probability:
                    # Accept the exchange.
                    self.replicaStateIndex[i] = sj
                    self.replicaStateIndex[j] = si
                    self.n_accepted[a][b] += 1

    _COUNTING_SAMPLER_CLASS = CountingReplicaExchangeSampler
    return _COUNTING_SAMPLER_CLASS


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

# metadata.json "npt" block for NVT runs (and the assumed value when extending
# a run directory written before NPT support existed).
_NPT_DISABLED_METADATA: dict[str, Any] = {
    "enabled": False,
    "p_type": None,
    "pressure_bar": None,
    "surface_tension_bar_nm": None,
}


class BuiltinREMDRunner(SimulationRunner):
    """Replica exchange runner based on the OpenMM >= 8.6 built-in sampler.

    All replicas share a single System built from the topology (the MBP
    coupling/offset values default to the topology values); replicas differ
    through the temperatures and global-parameter values in the sampler states.

    Attributes:
        config: Parsed simulation configuration object.
        output_data: Output directory (``remd_output``).
        states: State dicts passed to the sampler.
        sampler: The counting ReplicaExchangeSampler (set after run/extend).
    """

    def __init__(
        self,
        inpfile: str,
        replica_params: dict[str, Any] | None = None,
        output_data: str = "output_remd",
    ) -> None:
        """Initialize the runner.

        Args:
            inpfile: Path to input parameter file.
            replica_params: Optional dict with keys 'molname', 'beta', 'C1',
                'C2' (same layout as REMDRunner). If None, the values are taken
                from the input file's replica_* parameters.
            output_data: Output directory.

        Raises:
            ImportError: If OpenMM < 8.6 (no built-in replica exchange sampler).
        """
        _require_openmm86()
        super().__init__(inpfile)
        self.output_data = output_data
        if replica_params is not None:
            molname = replica_params["molname"]
            couplings = replica_params["beta"]
            c1_list = replica_params["C1"]
            c2_list = replica_params["C2"]
        else:
            molname = self.config.replica_molname
            couplings = self.config.replica_coupling
            c1_list = self.config.replica_c1
            c2_list = self.config.replica_c2
        self.molname: str = molname
        self.states: list[dict[str, Any]] = build_replica_states(
            molname,
            self.config.replica_method,
            self.config.replica_temp,
            couplings,
            c1_list,
            c2_list,
        )
        self.sampler: Any | None = None
        self.simulation: Simulation | None = None
        self._mbp_forces: list[mm.CustomCVForce] = []
        self._exchange_reporter: ExchangeStatisticsReporter | None = None
        self._reporters_closed = False
        if self.config.remd_unsampled_topfiles:
            print(
                "Warning: remd_unsampled_topfiles is ignored by the 'openmm' REMD "
                "backend; basin energies are recorded directly in cv_energies.csv."
            )

    def detect_existing_output(self) -> bool:
        """Detect if a previous run's output directory exists.

        Returns:
            True if the output directory exists and contains log.csv.
        """
        return os.path.isdir(self.output_data) and os.path.isfile(
            os.path.join(self.output_data, "log.csv")
        )

    def _npt_metadata(self) -> dict[str, Any]:
        """NPT settings block for metadata.json (single source for resume checks).

        Returns the disabled sentinel for NVT runs so old run directories
        without an "npt" key compare equal.
        """
        if self.config.pcouple != "yes":
            return dict(_NPT_DISABLED_METADATA)
        return {
            "enabled": True,
            "p_type": self.config.p_type,
            "pressure_bar": self.config.p_ref,
            "surface_tension_bar_nm": (
                self.config.p_tens if self.config.p_type == "membrane" else None
            ),
        }

    def _check_resume_metadata(self) -> None:
        """Validate that the existing run's metadata matches this input file.

        The states (temperatures and every parameter value), method, molname,
        replica count and report interval must be identical; nstep may differ.

        Raises:
            FileNotFoundError: If metadata.json is missing.
            ValueError: If any state definition differs, listing the
                differences.
        """
        metadata_path = os.path.join(self.output_data, "metadata.json")
        if not os.path.isfile(metadata_path):
            raise FileNotFoundError(
                f"Cannot extend REMD: {metadata_path} not found. The output "
                f"directory was not written by the 'openmm' backend (or is "
                f"incomplete); use run() to start a new simulation."
            )
        with open(metadata_path) as f:
            metadata = json.load(f)

        differences = []
        for key, current in (
            ("method", self.config.replica_method),
            ("molname", self.molname),
            ("replica_count", len(self.states)),
            ("reportInterval", self.config.remd_checkpoint_interval),
        ):
            if metadata.get(key) != current:
                differences.append(
                    f"{key}: existing {metadata.get(key)!r} vs input {current!r}"
                )
        existing_states = metadata.get("states", [])
        if len(existing_states) != len(self.states):
            differences.append(
                f"state count: existing {len(existing_states)} "
                f"vs input {len(self.states)}"
            )
        else:
            for i, (old, new) in enumerate(zip(existing_states, self.states, strict=True)):
                for key in sorted(set(old) | set(new)):
                    if old.get(key) != new.get(key):
                        differences.append(
                            f"state {i} parameter '{key}': existing "
                            f"{old.get(key)!r} vs input {new.get(key)!r}"
                        )
        existing_npt = metadata.get("npt", _NPT_DISABLED_METADATA)
        current_npt = self._npt_metadata()
        if existing_npt != current_npt:
            differences.append(f"npt: existing {existing_npt!r} vs input {current_npt!r}")
        if differences:
            raise ValueError(
                f"REMD state definition differs from the existing run in "
                f"{self.output_data}:\n  " + "\n  ".join(differences)
            )

    def _exchange_attempts(self) -> int:
        """Number of exchange iterations corresponding to nstep."""
        return int(math.floor(self.config.nstep / self.config.exc_freq))

    def _setup_simulation(self) -> Simulation:
        """Build the single System/Simulation shared by all replicas.

        Returns:
            Configured OpenMM Simulation object.
        """
        start_time = datetime.datetime.now()

        system, topology, _ = self.generate_topology()

        integrator = mm.LangevinIntegrator(
            self.config.temp * u.kelvin,
            self.config.fric_coeff / u.picosecond,
            self.config.dt * u.picosecond,
        )
        if self.config.const_tol:
            integrator.setConstraintTolerance(self.config.const_tol)
            print("Set Constraint Tolerance:", integrator.getConstraintTolerance())

        simulation = Simulation(
            topology, system, integrator, self.platform, self.platform_properties
        )

        conf, _ = load_structure(self.config.input)
        positions = conf.getPositions()
        assert len(positions) == topology.getNumAtoms(), (
            f"Error: Number of atoms in {self.config.input} is not the same as "
            f"that from {self.config.topol}!"
        )

        self._velocities_from_checkpoint = False
        if self.config.ichk:
            print(f"\n[Checkpoint] {self.config.ichk}", flush=True)
            with open(self.config.ichk) as f:
                state = mm.XmlSerializer.deserialize(f.read())
            simulation.context.setPositions(state.getPositions())
            simulation.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
            try:
                velocities = state.getVelocities()
            except mm.OpenMMException:
                velocities = None
            if velocities is not None:
                simulation.context.setVelocities(velocities)
                self._velocities_from_checkpoint = True
        else:
            simulation.context.setPositions(positions)

        report_time("REMD Setup", start_time)
        return simulation

    def _create_sampler(self, simulation: Simulation) -> Any:
        """Create the ReplicaExchangeSampler for the given Simulation.

        With remd_devices set, returns the parallel multi-Context driver (one
        worker Context per replica); otherwise the serial counting sampler.
        Under NPT both get the work-term exchange correction parameters.
        """
        validate_state_parameters(
            simulation.context, self.states, self.molname, self.config.replica_method
        )
        self._mbp_forces = find_mbp_forces(simulation.system, self.molname)
        if not self._mbp_forces:
            raise ValueError(
                f"No multiple-basin CustomCVForce found for molecule "
                f"'{self.molname}' in the system. Check replica_molname."
            )
        npt = self._npt_metadata()
        npt_kwargs: dict[str, Any] = {}
        if npt["enabled"]:
            npt_kwargs = {
                "pressure_bar": npt["pressure_bar"],
                "surface_tension_bar_nm": npt["surface_tension_bar_nm"] or 0.0,
            }
        if self.config.remd_devices:
            from ctgomartini.simulation.remd_parallel import (
                _get_parallel_sampler_class,
                assign_devices,
            )

            devices = assign_devices(len(self.states), self.config.remd_devices)
            base_seed = self.config.gen_seed
            if base_seed is None:
                base_seed = random.SystemRandom().randrange(2**31)
            print(
                f"  Parallel driver: {len(self.states)} worker Contexts, "
                f"replica -> device {devices}, worker integrator seeds "
                f"{base_seed}+i"
            )
            return _get_parallel_sampler_class()(
                self.states,
                simulation,
                self.config.exc_freq,
                reinitializeVelocities=self.config.remd_reassign_velocities,
                devices=devices,
                platform=self.platform,
                platform_properties=self.platform_properties,
                base_seed=base_seed,
                **npt_kwargs,
            )
        sampler_class = _get_counting_sampler_class()
        return sampler_class(
            self.states,
            simulation,
            self.config.exc_freq,
            reinitializeVelocities=self.config.remd_reassign_velocities,
            **npt_kwargs,
        )

    def _initialize_velocities(self, sampler: Any) -> None:
        """Assign per-replica velocities at each replica's initial state temperature.

        Uses gen_temp/gen_seed when gen_vel = yes, otherwise config.temp. A
        per-replica offset keeps seeded initializations distinct.
        """
        if self._velocities_from_checkpoint:
            return
        context = sampler.simulation.context
        if self.config.gen_vel == "yes":
            base_temperature = self.config.gen_temp
            seed = self.config.gen_seed
        else:
            base_temperature = self.config.temp
            seed = None
        for i, state_index in enumerate(sampler.replicaStateIndex):
            temperature = self.states[state_index].get("temperature", base_temperature)
            if seed is None:
                context.setVelocitiesToTemperature(temperature * u.kelvin)
            else:
                # Per-replica offset keeps seeded initializations distinct.
                context.setVelocitiesToTemperature(temperature * u.kelvin, seed + i)
            sampler.replicaConformation[i] = context.getState(
                getPositions=True,
                getVelocities=True,
                getParameters=True,
                getIntegratorParameters=True,
            )

    def _create_reporters(self, sampler: Any, resume: bool) -> None:
        """Attach the built-in and CTGoMartini reporters to the sampler."""
        from openmm.app.replicaexchangereporter import ReplicaExchangeReporter

        self._reporters_closed = False
        interval = self.config.remd_checkpoint_interval
        builtin_reporter = ReplicaExchangeReporter(
            self.output_data,
            interval,
            sampler,
            trajectoryPerReplica=(self.config.remd_xtc_output == "yes"),
            energy=True,
            volume=(self.config.pcouple == "yes"),
            checkpoints=True,
            resume=resume,
        )
        sampler.reporters.append(builtin_reporter)

        if len(self._mbp_forces) == 1:
            labels = [self.molname]
        else:
            labels = [f"{self.molname}_{k}" for k in range(len(self._mbp_forces))]
        labeled_forces = list(zip(labels, self._mbp_forces, strict=True))
        sampler.reporters.append(
            BasinEnergyReporter(
                os.path.join(self.output_data, "cv_energies.csv"),
                interval,
                labeled_forces,
                resume=resume,
            )
        )
        self._exchange_reporter = ExchangeStatisticsReporter(
            os.path.join(self.output_data, "exchange.csv"),
            interval,
            resume=resume,
        )
        sampler.reporters.append(self._exchange_reporter)

        if self.config.pcouple == "yes" and self.config.p_type == "membrane":
            # volume.csv is written by the built-in reporter; area.csv keeps
            # the same row-aligned layout for the membrane work term.
            sampler.reporters.append(
                AreaReporter(
                    os.path.join(self.output_data, "area.csv"),
                    interval,
                    resume=resume,
                )
            )

    def _finalize_reporters(self) -> None:
        """Flush the final exchange totals, then close all reporters."""
        if self._exchange_reporter is not None and self.sampler is not None:
            self._exchange_reporter.finalize(self.sampler)
        self._close_reporters()

    def _close_reporters(self) -> None:
        """Close every sampler reporter that exposes close() (idempotent)."""
        if self.sampler is None or self._reporters_closed:
            return
        for reporter in self.sampler.reporters:
            close = getattr(reporter, "close", None)
            if close is not None:
                close()
        self._reporters_closed = True

    def _write_metadata(self, report_interval: int) -> None:
        """Write metadata.json describing states and run parameters."""
        coupling_key = next(
            key for key in self.states[0] if key not in ("temperature",)
            and key.endswith(("_beta", "_delta"))
        )
        namespace = coupling_key.rsplit("_", 1)[0]
        metadata = {
            "backend": "openmm",
            "openmm_version": mm.__version__,
            "method": self.config.replica_method,
            "molname": self.molname,
            "replica_count": len(self.states),
            "parameter_names": {
                "coupling": coupling_key,
                "C1": f"{namespace}_C1",
                "C2": f"{namespace}_C2",
            },
            "states": self.states,
            "dt": self.config.dt,
            "exc_freq": self.config.exc_freq,
            "stepsPerIteration": self.config.exc_freq,
            "reportInterval": report_interval,
            "n_basins": max(
                len(basin_cv_groups(force)) for force in self._mbp_forces
            ),
            "devices": self.config.remd_devices,
            "parallel": bool(self.config.remd_devices),
            "npt": self._npt_metadata(),
        }
        worker_seeds = getattr(self.sampler, "worker_seeds", None)
        if worker_seeds is not None:
            metadata["worker_seeds"] = list(worker_seeds)
        metadata_path = os.path.join(self.output_data, "metadata.json")
        with open(metadata_path, "w") as f:
            json.dump(metadata, f, indent=2)

    def run(self) -> None:
        """Execute a new REMD simulation from scratch."""
        start_time = datetime.datetime.now()
        exchange_attempts = self._exchange_attempts()
        if exchange_attempts < 1:
            raise ValueError(
                f"nstep ({self.config.nstep}) must be >= exc_freq "
                f"({self.config.exc_freq}) for REMD: at least one exchange "
                f"iteration is required."
            )

        output_dir = Path(self.output_data)
        if output_dir.exists() and any(output_dir.iterdir()):
            raise ValueError(
                f"REMD output directory {self.output_data} exists and is not empty. "
                f"Use --append (extend) to continue it, or remove it manually."
            )

        self.simulation = self._setup_simulation()
        self.sampler = self._create_sampler(self.simulation)
        self._initialize_velocities(self.sampler)
        # The built-in reporter creates the (empty) directory; metadata and the
        # custom reporters are written afterwards.
        self._create_reporters(self.sampler, resume=False)
        self._write_metadata(self.config.remd_checkpoint_interval)

        total_steps = exchange_attempts * self.config.exc_freq
        print("\n[Production REMD] OpenMM built-in backend", flush=True)
        print(f"  Replicas: {len(self.states)}", flush=True)
        print(f"  Iterations: {exchange_attempts} ({total_steps:,} steps)", flush=True)
        print(f"  Output directory: {self.output_data}", flush=True)

        try:
            self.sampler.simulate(exchange_attempts)
            self._finalize_reporters()
        finally:
            self._close_reporters()
        report_time("REMD", start_time)

    def extend(self, n_iterations: int | None = None) -> None:
        """Continue an existing REMD simulation.

        Args:
            n_iterations: Number of additional iterations to run. If None,
                runs until the target implied by nstep is reached.
        """
        if not self.detect_existing_output():
            raise FileNotFoundError(
                f"Cannot extend REMD: output directory {self.output_data} with "
                f"log.csv does not exist. Use run() to start a new simulation."
            )
        self._check_resume_metadata()
        start_time = datetime.datetime.now()
        exchange_attempts = self._exchange_attempts()

        print(f"\n[Extend REMD] {self.output_data}", flush=True)
        self.simulation = self._setup_simulation()
        self.sampler = self._create_sampler(self.simulation)
        # resume=True restores replica conformations, state assignments and
        # currentIteration from checkpoint_*.xml and log.csv.
        self._create_reporters(self.sampler, resume=True)
        self.sampler.set_exchange_statistics(
            read_exchange_statistics(os.path.join(self.output_data, "exchange.csv"))
        )

        try:
            if n_iterations is not None:
                run_iterations = n_iterations
            else:
                run_iterations = exchange_attempts - self.sampler.currentIteration
                if run_iterations <= 0:
                    print(
                        f"  Simulation already complete: "
                        f"{self.sampler.currentIteration}/{exchange_attempts} iterations",
                        flush=True,
                    )
                    return
            print(
                f"  Progress: iteration {self.sampler.currentIteration} -> "
                f"{self.sampler.currentIteration + run_iterations}",
                flush=True,
            )
            self.sampler.simulate(run_iterations)
            self._finalize_reporters()
            report_time("REMD Extend", start_time)
        finally:
            self._close_reporters()
