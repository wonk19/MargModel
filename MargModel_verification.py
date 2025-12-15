# -*- coding: utf-8 -*-
"""
MargModel_verification.py

Verification script to confirm MargModel.py produces identical results
to the original K03 implementation.

This compares the standalone MargModel with the original gf_paper_surv19_alex module.
"""

import sys
import os
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_dir, 'src'))

import numpy as np
import MargModel
import gf_paper_surv19_alex as gf
import gv_ocean_optics as gvoo

def verify_forward_model():
    """
    Compare MargModel.runForward() with gf.get_ccrr_rrs_rt()
    """
    print("=" * 70)
    print("Verification: MargModel vs Original Implementation")
    print("=" * 70)
    
    # Test parameters
    wavelengths = np.arange(400., 750.)
    chla = 2.0
    mineral = 0.5
    aCDOM = 0.02
    
    # Model parameters (from K03)
    eta = 0.001
    g0 = 0.001
    nu = 0.0
    
    fresnel_add = 0.02
    sol_zen = 30.0
    
    # rho_sky (zero for simplicity)
    rho_sky = np.zeros_like(wavelengths)
    
    print(f"\nTest Parameters:")
    print(f"  Wavelengths: {wavelengths[0]:.0f} - {wavelengths[-1]:.0f} nm")
    print(f"  Chla: {chla} mg/m3")
    print(f"  Mineral: {mineral} g/m3")
    print(f"  aCDOM(443): {aCDOM} 1/m")
    print(f"  eta: {eta}, g0: {g0}, nu: {nu}")
    
    # Run MargModel (without fluorescence for comparison with original)
    print(f"\nRunning MargModel.runForward()...", end=" ")
    Rrs_new = MargModel.runForward(
        wavelengths, chla, mineral, aCDOM,
        eta_param=eta, g0_param=g0, nu_param=nu,
        fresnel_add=fresnel_add, rho_sky=rho_sky, sol_zen=sol_zen,
        add_fluorescence_flag=False  # Disable fluorescence for fair comparison
    )
    print("Done")
    
    # Run original implementation
    print(f"Running original gf.get_ccrr_rrs_rt()...", end=" ")
    gvoo.eta = eta
    gvoo.g0 = g0
    gvoo.nu = nu
    
    Rrs_old = gf.get_ccrr_rrs_rt(
        wavelengths, chla, mineral, aCDOM,
        fresnel_add, wavelengths, rho_sky, sol_zen
    )
    print("Done")
    
    # Compare results
    print(f"\n{'='*70}")
    print("Comparison Results:")
    print(f"{'='*70}")
    
    diff = Rrs_new - Rrs_old
    rel_diff = np.abs(diff / Rrs_old) * 100
    
    print(f"\nAbsolute Difference (MargModel - Original):")
    print(f"  Mean: {np.mean(diff)*1000:.6f} x10^-3 sr^-1")
    print(f"  Max: {np.max(np.abs(diff))*1000:.6f} x10^-3 sr^-1")
    print(f"  RMS: {np.sqrt(np.mean(diff**2))*1000:.6f} x10^-3 sr^-1")
    
    print(f"\nRelative Difference:")
    print(f"  Mean: {np.mean(rel_diff):.6f} %")
    print(f"  Max: {np.max(rel_diff):.6f} %")
    
    # Check if differences are within tolerance
    tolerance = 1e-10  # Very small tolerance for numerical precision
    max_abs_diff = np.max(np.abs(diff))
    
    if max_abs_diff < tolerance:
        print(f"\n*** VERIFICATION PASSED ***")
        print(f"Maximum difference ({max_abs_diff:.2e}) is within tolerance ({tolerance:.2e})")
    else:
        print(f"\n*** WARNING ***")
        print(f"Maximum difference ({max_abs_diff:.2e}) exceeds tolerance ({tolerance:.2e})")
        print(f"This may indicate implementation differences.")
    
    # Show some example values
    print(f"\n{'Wavelength (nm)':<20s} {'MargModel':<15s} {'Original':<15s} {'Diff':<15s}")
    print("-" * 70)
    test_wavelengths = [443, 490, 555, 670]
    for wl in test_wavelengths:
        idx = int(wl - 400)
        if idx < len(Rrs_new):
            print(f"{wl:<20d} {Rrs_new[idx]*1000:<15.6f} {Rrs_old[idx]*1000:<15.6f} {diff[idx]*1000:<15.6e}")
    
    return max_abs_diff < tolerance


def verify_multiple_conditions():
    """
    Test multiple water quality conditions
    """
    print(f"\n\n{'='*70}")
    print("Extended Verification: Multiple Conditions")
    print(f"{'='*70}")
    
    wavelengths = np.arange(400., 750.)
    rho_sky = np.zeros_like(wavelengths)
    
    # Test cases
    test_cases = [
        {'name': 'Clear water', 'chla': 0.5, 'mineral': 0.1, 'aCDOM': 0.005},
        {'name': 'Coastal water', 'chla': 2.0, 'mineral': 0.5, 'aCDOM': 0.02},
        {'name': 'Turbid water', 'chla': 5.0, 'mineral': 2.0, 'aCDOM': 0.05},
        {'name': 'High Chla', 'chla': 10.0, 'mineral': 0.5, 'aCDOM': 0.02},
        {'name': 'High mineral', 'chla': 2.0, 'mineral': 5.0, 'aCDOM': 0.02},
    ]
    
    eta, g0, nu = 0.001, 0.001, 0.0
    gvoo.eta, gvoo.g0, gvoo.nu = eta, g0, nu
    
    all_passed = True
    
    print(f"\n{'Test Case':<20s} {'Max Diff':<15s} {'Status':<10s}")
    print("-" * 70)
    
    for tc in test_cases:
        # Run both models (without fluorescence for fair comparison)
        Rrs_new = MargModel.runForward(
            wavelengths, tc['chla'], tc['mineral'], tc['aCDOM'],
            eta_param=eta, g0_param=g0, nu_param=nu,
            add_fluorescence_flag=False
        )
        
        Rrs_old = gf.get_ccrr_rrs_rt(
            wavelengths, tc['chla'], tc['mineral'], tc['aCDOM'],
            0.02, wavelengths, rho_sky, 30.0
        )
        
        # Compare
        max_diff = np.max(np.abs(Rrs_new - Rrs_old))
        tolerance = 1e-10
        passed = max_diff < tolerance
        all_passed = all_passed and passed
        
        status = "PASS" if passed else "FAIL"
        print(f"{tc['name']:<20s} {max_diff:<15.2e} {status:<10s}")
    
    print("-" * 70)
    
    if all_passed:
        print(f"\n*** ALL TESTS PASSED ***")
        return True
    else:
        print(f"\n*** SOME TESTS FAILED ***")
        return False


def verify_parameter_sensitivity():
    """
    Test different model parameters (eta, g0, nu)
    """
    print(f"\n\n{'='*70}")
    print("Parameter Sensitivity Verification")
    print(f"{'='*70}")
    
    wavelengths = np.arange(400., 750.)
    rho_sky = np.zeros_like(wavelengths)
    
    chla, mineral, aCDOM = 2.0, 0.5, 0.02
    
    # Test different parameter sets
    param_sets = [
        {'eta': 0.0001, 'g0': 0.001, 'nu': 0.0},
        {'eta': 0.001, 'g0': 0.001, 'nu': 0.0},
        {'eta': 0.005, 'g0': 0.001, 'nu': 0.0},
        {'eta': 0.001, 'g0': 0.0001, 'nu': 0.0},
        {'eta': 0.001, 'g0': 0.01, 'nu': 0.0},
    ]
    
    all_passed = True
    
    print(f"\n{'Parameters':<30s} {'Max Diff':<15s} {'Status':<10s}")
    print("-" * 70)
    
    for params in param_sets:
        eta, g0, nu = params['eta'], params['g0'], params['nu']
        param_str = f"eta={eta:.4f}, g0={g0:.4f}"
        
        # Run both models (without fluorescence for fair comparison)
        Rrs_new = MargModel.runForward(
            wavelengths, chla, mineral, aCDOM,
            eta_param=eta, g0_param=g0, nu_param=nu,
            add_fluorescence_flag=False
        )
        
        gvoo.eta, gvoo.g0, gvoo.nu = eta, g0, nu
        Rrs_old = gf.get_ccrr_rrs_rt(
            wavelengths, chla, mineral, aCDOM,
            0.02, wavelengths, rho_sky, 30.0
        )
        
        # Compare
        max_diff = np.max(np.abs(Rrs_new - Rrs_old))
        tolerance = 1e-10
        passed = max_diff < tolerance
        all_passed = all_passed and passed
        
        status = "PASS" if passed else "FAIL"
        print(f"{param_str:<30s} {max_diff:<15.2e} {status:<10s}")
    
    print("-" * 70)
    
    if all_passed:
        print(f"\n*** ALL PARAMETER TESTS PASSED ***")
        return True
    else:
        print(f"\n*** SOME PARAMETER TESTS FAILED ***")
        return False


if __name__ == "__main__":
    """
    Run all verification tests
    """
    print("\n" + "=" * 70)
    print("MargModel.py Verification Suite")
    print("=" * 70)
    print("\nThis script verifies that MargModel.py produces identical results")
    print("to the original K03 implementation (gf_paper_surv19_alex.py).\n")
    
    # Run verification tests
    test1 = verify_forward_model()
    test2 = verify_multiple_conditions()
    test3 = verify_parameter_sensitivity()
    
    # Final summary
    print(f"\n\n{'='*70}")
    print("VERIFICATION SUMMARY")
    print(f"{'='*70}")
    
    tests = [
        ('Basic Forward Model', test1),
        ('Multiple Conditions', test2),
        ('Parameter Sensitivity', test3),
    ]
    
    all_passed = True
    for test_name, result in tests:
        status = "PASS" if result else "FAIL"
        print(f"  {test_name:<30s}: {status}")
        all_passed = all_passed and result
    
    print(f"{'='*70}")
    
    if all_passed:
        print("\n*** ALL VERIFICATION TESTS PASSED ***")
        print("\nMargModel.py is confirmed to produce identical results")
        print("to the original K03 implementation.")
    else:
        print("\n*** SOME VERIFICATION TESTS FAILED ***")
        print("\nPlease review the differences above.")
    
    print("\n")

