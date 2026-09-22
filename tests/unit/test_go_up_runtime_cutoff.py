import pytest

from ctgomartini.topology.generator.combiner import combine_contacts
from ctgomartini.topology.interactions.bonds import Contacts6_12
from ctgomartini.topology.interactions.multi_state import MultiAllBonds


def test_combine_contacts_writes_runtime_cutoff() -> None:
    state_contacts = [
        [
            ['1', '5', '1', '0.50', '15.0'],
            ['2', '8', '1', '0.70', '15.0'],
        ],
        [
            ['1', '5', '1', '0.54', '15.0'],
            ['3', '9', '1', '0.80', '15.0'],
        ],
    ]

    contacts, multi_contacts = combine_contacts(
        2,
        state_contacts,
        cutoff=0.06,
        contact_rcut=1.5,
    )

    assert contacts == [['1', '5', '2', '0.52', '15.0', '1.5']]
    assert multi_contacts == [
        ['2', '8', '2', '1', '2', '0.70', '15.0', '1.5'],
        ['3', '9', '2', '2', '2', '0.80', '15.0', '1.5'],
    ]


def test_combine_contacts_keeps_legacy_format_when_cutoff_is_omitted() -> None:
    state_contacts = [
        [['1', '5', '1', '0.50', '15.0']],
        [['1', '5', '1', '0.54', '15.0']],
    ]

    contacts, multi_contacts = combine_contacts(
        2,
        state_contacts,
        cutoff=0.06,
    )

    assert contacts == [['1', '5', '1', '0.52', '15.0']]
    assert multi_contacts == []


def test_adjustable_contact_force_reads_runtime_cutoff() -> None:
    interaction = Contacts6_12()
    interaction.add_interaction(['1', '2', '2', '0.55', '15.0', '1.5'])

    _, _, parameters = interaction.mm_force.getBondParameters(0)

    assert float(parameters[2]) == pytest.approx(1.5)


def test_multi_contact_force_reads_runtime_cutoff() -> None:
    interaction = MultiAllBonds()
    fields = ['1', '2', '2', '1', '2', '0.55', '15.0', '1.5']

    assert interaction.add_interaction('1', 'multi_contacts', fields)

    force = interaction.mm_force
    parameter_names = [
        force.getPerBondParameterName(index)
        for index in range(force.getNumPerBondParameters())
    ]
    _, parameters = force.getBondParameters(0)

    assert 'rcut' in parameter_names
    assert float(parameters[parameter_names.index('rcut')]) == pytest.approx(1.5)
