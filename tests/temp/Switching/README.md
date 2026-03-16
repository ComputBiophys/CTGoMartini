# Switching Go-Martini Test

This directory contains a test case for the Switching Go-Martini mode using vermouth >= 0.15.0.

## Files

### Input Files
- `State1.pdb`, `State1.map` - First state structure and contact map
- `State2.pdb`, `State2.map` - Second state structure and contact map
- `martini_v3.0.0.itp` - Martini 3.0 force field

### Scripts
- `run_switching.py` - Main script to generate Switching Go-Martini topology
- `verify_contacts.py` - Verification script to check contacts conversion

### Output Files (generated)
- `Open/` - Directory containing Open state files
  - `Open.itp` - Final topology with [contacts] section
  - `go_atomtypes.itp`, `go_nbparams.itp` - Go model parameters
  - `Open_cg.pdb` - Coarse-grained structure
- `Closed/` - Directory containing Closed state files
  - `Closed.itp` - Final topology with [contacts] section
  - Similar Go model files as Open state

## Usage

### 1. Run Switching Go-Martini Generation

```bash
python run_switching.py
```

This will:
1. Convert `.map` files to `contact_map.out` format
2. Run `martinize2 -go` to generate Go contacts
3. Convert LJ `nonbond_params` to `[contacts]` bonds
4. Generate `Open.itp` and `Closed.itp`

### 2. Verify Results

```bash
python verify_contacts.py
```

This compares:
- `go_nbparams.itp` (LJ nonbond_params) 
- `Open.itp` / `Closed.itp` ([contacts] section)

Expected output:
```
【Open 状态】
  go_nbparams.itp 接触对: 479
  Open.itp contacts: 479
  ✅ 完全匹配 (479 contacts)

【Closed 状态】
  go_nbparams.itp 接触对: 516
  Closed.itp contacts: 516
  ✅ 完全匹配 (516 contacts)
```

## New Workflow (vermouth >= 0.15.0)

### Old Workflow (0.9.6)
```
martinize2 -dssp <cmd> -govs-include ...
create_go_virt_for_multimer -f contact.map
→ BB-part-def_VirtGoSites.itp, go-table_VirtGoSites.itp
```

### New Workflow (0.15.0)
```
convert_map_format input.map output.out
martinize2 -dssp -go output.out -go-low 0.3 -go-up 1.1 -go-eps 12.0
→ go_atomtypes.itp, go_nbparams.itp
GenSBPTop converts nonbond_params → [contacts]
```

## Key Differences

| Feature | Old (0.9.6) | New (0.15.0) |
|---------|-------------|--------------|
| DSSP | External command | MDTraj (built-in) |
| Side-chain fix | `-scfix` to enable | Enabled by default |
| Go contacts | `create_go_virt_for_multimer` | `martinize2 -go` |
| Output files | `go-table_VirtGoSites.itp` | `go_nbparams.itp` |
| Contact format | LJ nonbond_params | [contacts] bonds |
