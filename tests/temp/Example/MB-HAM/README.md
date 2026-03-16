# MB-HAM Example

This example demonstrates Multiple-Basin HAM topology generation.

## Method

**HAM (Harmonic Approximation)**: Uses harmonic potential for multiple-basin systems.

## Files

- `State1.pdb`, `State1.map` - Starting conformation and contacts
- `State2.pdb`, `State2.map` - Target conformation and contacts
- `martini_v3.0.0.itp` - Martini force field
- `generate_topology.py` - Topology generation script

## Usage

```bash
python generate_topology.py
```

This creates `Protein.itp` and `Protein_params.itp`.

## Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `cutoff_BBB_angles` | 15.0° | Angle tolerance |
| `cutoff_BBBB_dihedrals` | 30.0° | Dihedral tolerance |
| `cutoff_contacts` | 0.06 | Contact overlap tolerance |

## Method: HAM

Uses harmonic approximation for multiple-basin potential.
