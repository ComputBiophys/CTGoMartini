"""Multi-Context parallel driver for the built-in REMD backend (OpenMM >= 8.6).

Selected by setting ``remd_devices`` in the input file (used only by
``remd_backend = openmm``). One worker Context is created per replica; each
worker owns its Context, LangevinIntegrator and ``MultistateSampler``, so all
replicas of an iteration are integrated concurrently in a thread pool (OpenMM
releases the GIL while integrating). The main Context is never stepped by the
parallel driver and stays available to the reporters, so output files and
formats are identical to the serial driver.

The per-iteration work mirrors ``ReplicaExchangeSampler.simulateReplica()``
plus ``MultistateSampler.computeAllEnergies()`` exactly (velocity
reinitialization on state change, sqrt(kT2/k1) rescaling on temperature
change, ``_previousReplicaStateIndex`` bookkeeping); the iteration structure
(reporters after ``currentIteration`` is incremented, then
``exchangeReplicas()``) is unchanged. Worker integrators get independent
random seeds ``base_seed + i`` with ``base_seed = gen_seed`` (or a random
draw, recorded in ``metadata.json`` as ``worker_seeds``). The NPT
work-term exchange correction (``pressure_bar``/``surface_tension_bar_nm``)
is inherited from the counting sampler through ``**kwargs``. Since the
workers step their own Contexts, the main Context's step count is synced to
the last replica's conformation count before reporters run, matching the
serial driver's ``log.csv`` Step semantics.

The OpenMM 8.6 sampler classes are imported lazily inside functions so this
module stays importable under OpenMM < 8.6.
"""

from __future__ import annotations

import math
from concurrent.futures import ThreadPoolExecutor
from typing import Any

import openmm as mm
import openmm.unit as u


def assign_devices(n_replicas: int, devices: list[int]) -> list[int]:
    """Assign each replica a device, round-robin over the given device list.

    Args:
        n_replicas: Number of replicas (one worker Context each).
        devices: Device indices from ``remd_devices``; must be non-empty.

    Returns:
        List of length n_replicas: replica i uses devices[i % len(devices)].

    Raises:
        ValueError: If devices is empty.
    """
    if not devices:
        raise ValueError("remd_devices must not be empty for the parallel driver.")
    return [devices[i % len(devices)] for i in range(n_replicas)]


def worker_platform_properties(
    platform: mm.Platform, base_properties: dict[str, str], device: int
) -> dict[str, str]:
    """Platform properties for a worker Context bound to the given device.

    ``DeviceIndex`` is set (overriding any base value) only for the GPU
    platforms CUDA and OpenCL; CPU and Reference ignore device assignment.

    Args:
        platform: The OpenMM Platform shared by all workers.
        base_properties: Base platform properties from the runner.
        device: Device index assigned to this worker.

    Returns:
        A new properties dict for the worker Context.
    """
    properties = dict(base_properties)
    if platform.getName() in ("CUDA", "OpenCL"):
        properties["DeviceIndex"] = str(device)
    return properties


_PARALLEL_SAMPLER_CLASS = None


def _get_parallel_sampler_class() -> type:
    """Return the parallel ReplicaExchangeSampler subclass (lazy, OpenMM >= 8.6).

    The class subclasses the counting sampler from remd_openmm (per-pair
    exchange counters) and overrides ``simulate()`` to integrate every replica
    on its own worker Context in a shared thread pool.
    """
    global _PARALLEL_SAMPLER_CLASS
    if _PARALLEL_SAMPLER_CLASS is not None:
        return _PARALLEL_SAMPLER_CLASS

    from ctgomartini.simulation.remd_openmm import (
        _get_counting_sampler_class,
        _require_openmm86,
    )

    _require_openmm86()
    from openmm.app.internal.multistatesampler import MultistateSampler

    class _Worker:
        """Per-replica worker: its own Context, integrator and MultistateSampler."""

        def __init__(
            self,
            context: mm.Context,
            integrator: mm.LangevinIntegrator,
            multistate_sampler: Any,
            device: int,
            seed: int,
        ) -> None:
            self.context = context
            self.integrator = integrator
            self.multistate_sampler = multistate_sampler
            self.device = device
            self.seed = seed

    class ParallelReplicaExchangeSampler(_get_counting_sampler_class()):
        """ReplicaExchangeSampler driving one worker Context per replica.

        Attributes:
            worker_devices: Device index assigned to each replica (round-robin).
            worker_seeds: Integrator random seed of each worker (base_seed + i).
        """

        def __init__(
            self,
            *args: Any,
            devices: list[int],
            platform: mm.Platform,
            platform_properties: dict[str, str],
            base_seed: int,
            **kwargs: Any,
        ) -> None:
            super().__init__(*args, **kwargs)
            main_integrator = self.simulation.integrator
            self.worker_devices = list(devices)
            self.worker_seeds = [base_seed + i for i in range(len(self.states))]
            self._workers: list[_Worker] = []
            for i in range(len(self.states)):
                integrator = mm.LangevinIntegrator(
                    main_integrator.getTemperature(),
                    main_integrator.getFriction(),
                    main_integrator.getStepSize(),
                )
                integrator.setConstraintTolerance(
                    main_integrator.getConstraintTolerance()
                )
                integrator.setRandomNumberSeed(self.worker_seeds[i])
                properties = worker_platform_properties(
                    platform, platform_properties, devices[i]
                )
                if properties:
                    context = mm.Context(
                        self.simulation.system, integrator, platform, properties
                    )
                else:
                    context = mm.Context(self.simulation.system, integrator, platform)
                mss = MultistateSampler(self.states, context)
                self._workers.append(
                    _Worker(context, integrator, mss, devices[i], self.worker_seeds[i])
                )

        def _simulate_replica(self, index: int) -> Any:
            """Simulate one replica on its worker Context; return its energies.

            Mirrors ReplicaExchangeSampler.simulateReplica() followed by
            MultistateSampler.computeAllEnergies(), but runs entirely on the
            worker's own Context/integrator so threads never share mutable
            OpenMM objects. Called only from the thread pool inside simulate().
            """
            worker = self._workers[index]
            worker.context.setState(self.replicaConformation[index])
            if self.reinitializeVelocities:
                # Reinitialize velocities for replicas whose states have changed.
                previous = self._previousReplicaStateIndex[index]
                if previous != self.replicaStateIndex[index]:
                    state = self.states[self.replicaStateIndex[index]]
                    if "temperature" in state:
                        temperature = state["temperature"]
                    else:
                        temperature = worker.integrator.getTemperature()
                    worker.context.setVelocitiesToTemperature(temperature)
            elif self._kT is not None:
                # Rescale velocities for replicas whose temperatures have changed.
                kt1 = self._kT[self._previousReplicaStateIndex[index]]
                kt2 = self._kT[self.replicaStateIndex[index]]
                if kt1 != kt2:
                    scale = math.sqrt(kt2 / kt1)
                    worker.context.setVelocities(
                        self.replicaConformation[index].getVelocities(asNumpy=True)
                        * scale
                    )
            self._previousReplicaStateIndex[index] = self.replicaStateIndex[index]
            worker.multistate_sampler.applyState(self.replicaStateIndex[index])
            worker.integrator.step(self.stepsPerIteration)
            self.replicaConformation[index] = worker.context.getState(
                positions=True,
                velocities=True,
                parameters=True,
                integratorParameters=True,
            )
            return worker.multistate_sampler.computeAllEnergies()

        def simulate(self, iterations: int) -> None:
            """Run iterations, integrating all replicas in parallel threads.

            The iteration structure is identical to the serial sampler:
            per-replica simulation and energy evaluation (parallel), then
            reporters, then exchange attempts (both on the main thread).
            Worker exceptions propagate through their futures.
            """
            n = len(self.states)
            with ThreadPoolExecutor(max_workers=n) as pool:
                for _ in range(iterations):
                    self.replicaStateEnergy = None
                    futures = [pool.submit(self._simulate_replica, i) for i in range(n)]
                    energies = [] * u.kilojoules_per_mole
                    for future in futures:
                        energies.append(future.result())
                    self.replicaStateEnergy = energies
                    self.currentIteration += 1
                    # Reporters read the main Context's step count (log.csv's
                    # Step column). The serial driver's main Context ends each
                    # iteration holding the last replica's conformation count;
                    # mirror that here since the workers step their own
                    # Contexts. Positions/parameters are left untouched.
                    self.simulation.context.setStepCount(
                        self.replicaConformation[-1].getStepCount()
                    )
                    for reporter in self.reporters:
                        reporter(self)
                    self.exchangeReplicas()

    _PARALLEL_SAMPLER_CLASS = ParallelReplicaExchangeSampler
    return _PARALLEL_SAMPLER_CLASS
