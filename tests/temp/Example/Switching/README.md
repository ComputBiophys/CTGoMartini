# Switching Example

Switching potential topology generation for conformational transitions.

## Overview

Switching mode generates independent single-basin topologies for each state.
Users manually switch between states during simulation by changing the topology.

## Files

- `State1.pdb`, `State1.map` - Starting conformation and contacts
- `State2.pdb`, `State2.map` - Target conformation and contacts
- `martini_v3.0.0.itp` - Martini 3.0 force field
- `generate_topology.py` - Topology generation script

## Usage

```bash
python generate_topology.py
```

Generates:
- `Open.itp`, `Open_params.itp` - State 1 topology
- `Closed.itp`, `Closed_params.itp` - State 2 topology

## Simulation Strategy

1. Start simulation with Open.itp
2. Run equilibrium simulation
3. Switch to Closed.itp for target state simulation
4. Compare results between states
