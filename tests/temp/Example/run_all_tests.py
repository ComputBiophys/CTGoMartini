"""Test all CTGoMartini modes."""
import os
import sys
import shutil

sys.path.insert(0, '/home/ys/CommonUse/Martini/CTGoMartini')

from ctgomartini.data.ctgomartinize import SBGOMartinize, SwitchingGOMartinize, MBGOMartinize


def clean_dir(path, keep_files=None):
    """Remove generated directories and files, keep inputs."""
    if keep_files is None:
        keep_files = ['.pdb', '.map', '.py', '.md']
    
    for item in os.listdir(path):
        item_path = os.path.join(path, item)
        if os.path.isdir(item_path):
            shutil.rmtree(item_path)
        elif os.path.isfile(item_path):
            ext = os.path.splitext(item)[1].lower()
            if ext not in keep_files:
                os.remove(item_path)


def test_sb_mode():
    """Test Single-Basin mode."""
    print("\n" + "="*70)
    print("TEST 1: Single-Basin (SB) Mode")
    print("="*70)
    
    os.chdir('/home/ys/CommonUse/Martini/CTGoMartini/tests/temp/Example/SB')
    clean_dir('.')
    
    SBGOMartinize(
        aa_strfile_list=['HAMP_Bicyclomycin.pdb'],
        aa_map_list=['HAMP_Bicyclomycin.map'],
        state_name_list=['Protein'],
        sbmol_name='Protein',
        method='SBP',
        dssp=None,
        ff='martini3001',
        other_params=''
    )
    
    assert os.path.exists('Protein.itp'), "Protein.itp not generated"
    print("✓ SB mode test passed - Protein.itp generated")


def test_switching_mode():
    """Test Switching mode - generates independent state topologies."""
    print("\n" + "="*70)
    print("TEST 2: Switching Mode")
    print("="*70)
    
    os.chdir('/home/ys/CommonUse/Martini/CTGoMartini/tests/temp/Example/Switching')
    clean_dir('.')
    
    SwitchingGOMartinize(
        aa_strfile_list=['State1.pdb', 'State2.pdb'],
        aa_map_list=['State1.map', 'State2.map'],
        state_name_list=['Open', 'Closed'],
        mbmol_name='Protein',
        dict_cutoffs={
            'cutoff_BBB_angles': 15.0,
            'cutoff_BBBB_dihedrals': 30.0,
            'cutoff_SBBS_dihedrals': 30.0,
            'cutoff_contacts': 0.06
        },
        method='switching',
        dssp=None,
        ff='martini3001',
        other_params=''
    )
    
    # Switching mode generates separate state topologies
    assert os.path.exists('Open.itp'), "Open.itp not generated"
    assert os.path.exists('Closed.itp'), "Closed.itp not generated"
    print("✓ Switching mode test passed - Open.itp and Closed.itp generated")


def test_mb_exp_mode():
    """Test MB EXP mode."""
    print("\n" + "="*70)
    print("TEST 3: MB-EXP Mode")
    print("="*70)
    
    os.chdir('/home/ys/CommonUse/Martini/CTGoMartini/tests/temp/Example/MB-EXP')
    clean_dir('.')
    
    MBGOMartinize(
        aa_strfile_list=['State1.pdb', 'State2.pdb'],
        aa_map_list=['State1.map', 'State2.map'],
        state_name_list=['State1', 'State2'],
        mbmol_name='Protein',
        dict_cutoffs={
            'cutoff_BBB_angles': 15.0,
            'cutoff_BBBB_dihedrals': 30.0,
            'cutoff_SBBS_dihedrals': 30.0,
            'cutoff_contacts': 0.06
        },
        method='exp',
        dssp=None,
        ff='martini3001',
        other_params=''
    )
    
    assert os.path.exists('Protein.itp'), "Protein.itp not generated"
    
    with open('Protein.itp') as f:
        content = f.read()
        assert 'exp' in content, "EXP method not found"
    
    print("✓ MB-EXP mode test passed - Protein.itp with 'exp' method")


def test_mb_ham_mode():
    """Test MB HAM mode."""
    print("\n" + "="*70)
    print("TEST 4: MB-HAM Mode")
    print("="*70)
    
    os.chdir('/home/ys/CommonUse/Martini/CTGoMartini/tests/temp/Example/MB-HAM')
    clean_dir('.')
    
    MBGOMartinize(
        aa_strfile_list=['State1.pdb', 'State2.pdb'],
        aa_map_list=['State1.map', 'State2.map'],
        state_name_list=['Open', 'Closed'],
        mbmol_name='Protein',
        dict_cutoffs={
            'cutoff_BBB_angles': 15.0,
            'cutoff_BBBB_dihedrals': 30.0,
            'cutoff_SBBS_dihedrals': 30.0,
            'cutoff_contacts': 0.06
        },
        method='ham',
        dssp=None,
        ff='martini3001',
        other_params=''
    )
    
    assert os.path.exists('Protein.itp'), "Protein.itp not generated"
    
    with open('Protein.itp') as f:
        content = f.read()
        assert 'ham' in content, "HAM method not found"
    
    print("✓ MB-HAM mode test passed - Protein.itp with 'ham' method")


def main():
    """Run all tests."""
    print("="*70)
    print("CTGoMartini Test Suite - 4 Modes (SB/Switching/MB-EXP/MB-HAM)")
    print("="*70)
    
    original_dir = os.getcwd()
    tests = [
        test_sb_mode,
        test_switching_mode,
        test_mb_exp_mode,
        test_mb_ham_mode
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            test()
            passed += 1
        except Exception as e:
            failed += 1
            print(f"\n✗ {test.__name__} failed: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "="*70)
    print(f"Results: {passed} passed, {failed} failed")
    if failed == 0:
        print("All tests passed successfully!")
    print("="*70)
    
    os.chdir(original_dir)
    return 0 if failed == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
