"""Generator-level tests for the global ``contact_type`` switch.

``contact_type`` selects the Go contact potential written by the topology
generators (combiner.combine_contacts, single_basin.create_sb_topology,
multi_basin.create_mb_topology and the ctgomartinize CLI):

- type 1: 6-12 LJ with fixed 1.1 nm cutoff (legacy functype 1, no rcut
  field); an explicit non-default contact_rcut is ignored with a notice.
- type 2 (default): 6-12 LJ, functype 2, per-contact rcut (contact_rcut if
  given, else 1.1) - identical to the pre-switch CLI default output.
- type 3: 10-12 LJ, functype 3, per-contact rcut (same rule as type 2).
"""

import sys

import pytest

import ctgomartini.cli.ctgomartinize as ctgomartinize
from ctgomartini.topology.generator.combiner import _resolve_contact_format, combine_contacts
from ctgomartini.topology.generator.single_basin import (
    _apply_contact_format,
    _apply_contact_rcut,
)

# Two-state contacts: pair (1,5) is shared (|dSigma| <= cutoff), pairs
# (2,8)/(3,9) are state-specific and land in multi_contacts.
TWO_STATE_CONTACTS = [
    [
        ['1', '5', '1', '0.50', '15.0'],
        ['2', '8', '1', '0.70', '15.0'],
    ],
    [
        ['1', '5', '1', '0.54', '15.0'],
        ['3', '9', '1', '0.80', '15.0'],
    ],
]


# --------------------------------------------------------------------------- format resolution
def test_resolve_contact_format_mapping() -> None:
    assert _resolve_contact_format(1, None) == ('1', None)
    assert _resolve_contact_format(1, 1.1) == ('1', None)
    assert _resolve_contact_format(2, None) == ('2', '1.1')
    assert _resolve_contact_format(2, 1.5) == ('2', '1.5')
    assert _resolve_contact_format(3, None) == ('3', '1.1')
    assert _resolve_contact_format(3, 1.5) == ('3', '1.5')


def test_resolve_contact_format_type1_ignores_rcut_with_notice(capsys) -> None:
    assert _resolve_contact_format(1, 1.5) == ('1', None)
    assert 'contact_type=1' in capsys.readouterr().out


@pytest.mark.parametrize("bad", [0, 4, -1])
def test_resolve_contact_format_rejects_invalid_type(bad: int) -> None:
    with pytest.raises(ValueError, match="unsupported contact_type"):
        _resolve_contact_format(bad, None)


# --------------------------------------------------------------------------- combine_contacts
def test_combine_contacts_type_2_default_matches_current_output() -> None:
    """Default contact_type=2 with an explicit rcut is byte-identical to the
    pre-switch adjustable output; with rcut omitted it defaults to 1.1."""
    reference = combine_contacts(2, TWO_STATE_CONTACTS, cutoff=0.06, contact_rcut=1.5)
    explicit = combine_contacts(
        2, TWO_STATE_CONTACTS, cutoff=0.06, contact_rcut=1.5, contact_type=2
    )
    assert reference == explicit
    assert reference == (
        [['1', '5', '2', '0.52', '15.0', '1.5']],
        [
            ['2', '8', '2', '1', '2', '0.70', '15.0', '1.5'],
            ['3', '9', '2', '2', '2', '0.80', '15.0', '1.5'],
        ],
    )

    default_rcut = combine_contacts(2, TWO_STATE_CONTACTS, cutoff=0.06, contact_type=2)
    assert default_rcut == (
        [['1', '5', '2', '0.52', '15.0', '1.1']],
        [
            ['2', '8', '2', '1', '2', '0.70', '15.0', '1.1'],
            ['3', '9', '2', '2', '2', '0.80', '15.0', '1.1'],
        ],
    )


def test_combine_contacts_type_3_fields() -> None:
    """functype 3 keeps the functype-2 field layout (sigma/epsilon + rcut);
    the sigma/epsilon -> C10/C12 conversion happens in the handler."""
    contacts, multi_contacts = combine_contacts(
        2, TWO_STATE_CONTACTS, cutoff=0.06, contact_rcut=1.3, contact_type=3
    )
    assert contacts == [['1', '5', '3', '0.52', '15.0', '1.3']]
    assert multi_contacts == [
        ['2', '8', '2', '1', '3', '0.70', '15.0', '1.3'],
        ['3', '9', '2', '2', '3', '0.80', '15.0', '1.3'],
    ]


def test_combine_contacts_type_1_fixed_cutoff() -> None:
    contacts, multi_contacts = combine_contacts(
        2, TWO_STATE_CONTACTS, cutoff=0.06, contact_type=1
    )
    assert contacts == [['1', '5', '1', '0.52', '15.0']]
    assert multi_contacts == [
        ['2', '8', '2', '1', '1', '0.70', '15.0'],
        ['3', '9', '2', '2', '1', '0.80', '15.0'],
    ]


def test_combine_contacts_type_1_ignores_rcut_with_notice(capsys) -> None:
    contacts, _ = combine_contacts(
        2, TWO_STATE_CONTACTS, cutoff=0.06, contact_rcut=1.5, contact_type=1
    )
    assert contacts == [['1', '5', '1', '0.52', '15.0']]
    assert 'ignoring contact_rcut=1.5' in capsys.readouterr().out


def test_combine_contacts_rejects_invalid_type() -> None:
    with pytest.raises(ValueError, match="unsupported contact_type"):
        combine_contacts(2, TWO_STATE_CONTACTS, cutoff=0.06, contact_type=5)


# --------------------------------------------------------------------------- SB contact conversion
def test_apply_contact_format_type_3() -> None:
    contacts = [['1', '5', '1', '0.50', '15.0'], ['2', '8', '1', '0.70', '12.0']]
    assert _apply_contact_format(contacts, 3, 1.5) == [
        ['1', '5', '3', '0.50', '15.0', '1.5'],
        ['2', '8', '3', '0.70', '12.0', '1.5'],
    ]


def test_apply_contact_format_type_1_passthrough() -> None:
    contacts = [['1', '5', '1', '0.50', '15.0']]
    assert _apply_contact_format(contacts, 1) == contacts
    assert _apply_contact_format(contacts, 1, 1.5) == contacts  # rcut ignored


def test_apply_contact_format_default_matches_legacy_helper() -> None:
    """type 2 with/without explicit rcut reproduces _apply_contact_rcut at
    the same cutoff - the pre-switch default output is unchanged."""
    contacts = [['1', '5', '1', '0.50', '15.0'], ['2', '8', '1', '0.70', '12.0']]
    assert _apply_contact_format(contacts, 2, 1.5) == _apply_contact_rcut(contacts, 1.5)
    assert _apply_contact_format(contacts, 2) == _apply_contact_rcut(contacts, 1.1)


def test_apply_contact_format_rejects_invalid_type() -> None:
    with pytest.raises(ValueError, match="unsupported contact_type"):
        _apply_contact_format([['1', '5', '1', '0.50', '15.0']], 0)


# --------------------------------------------------------------------------- CLI plumbing
def test_ctgomartinize_forwards_contact_type_to_all_methods(monkeypatch) -> None:
    for method, fn_name in [
        ('sbp', 'SBGOMartinize'),
        ('exp', 'MBGOMartinize'),
        ('ham', 'MBGOMartinize'),
        ('switching', 'SwitchingGOMartinize'),
    ]:
        captured: dict = {}
        monkeypatch.setattr(
            ctgomartinize, fn_name, lambda *a, _c=captured, **k: _c.update(k)
        )
        ctgomartinize.CTGOMartinize(
            ['a.pdb'], ['auto'], ['A'], 'mbmol', {}, method=method, contact_type=3
        )
        assert captured['contact_type'] == 3, method


def test_cli_contact_type_argument(monkeypatch) -> None:
    captured: dict = {}
    monkeypatch.setattr(
        ctgomartinize, 'CTGOMartinize', lambda *a, **k: captured.update(k)
    )
    monkeypatch.setattr(sys, 'argv', [
        'ctgomartinize', '-s', 'a.pdb', '-mol', 'A', '-method', 'sbp',
        '-contact-type', '3',
    ])
    ctgomartinize.main()
    assert captured['contact_type'] == 3


def test_cli_contact_type_default_is_2(monkeypatch) -> None:
    captured: dict = {}
    monkeypatch.setattr(
        ctgomartinize, 'CTGOMartinize', lambda *a, **k: captured.update(k)
    )
    monkeypatch.setattr(
        sys, 'argv', ['ctgomartinize', '-s', 'a.pdb', '-mol', 'A', '-method', 'sbp']
    )
    ctgomartinize.main()
    assert captured['contact_type'] == 2


def test_cli_contact_type_rejects_invalid_value(monkeypatch) -> None:
    monkeypatch.setattr(sys, 'argv', [
        'ctgomartinize', '-s', 'a.pdb', '-mol', 'A', '-method', 'sbp',
        '-contact-type', '7',
    ])
    with pytest.raises(SystemExit):
        ctgomartinize.main()
