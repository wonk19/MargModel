# -*- coding: utf-8 -*-
"""
MargModel_example.py

Simple example script demonstrating how to use MargModel.py
to compute Rrs spectra from water quality parameters.

This script can be given to third parties as a starting point.
"""

import numpy as np
import matplotlib.pyplot as plt
import MargModel

def example_single_spectrum():
    """
    Example 1: Compute Rrs for a single set of water quality parameters
    """
    print("=" * 70)
    print("Example 1: Single Rrs Spectrum")
    print("=" * 70)
    
    # Define wavelength range (400-750 nm)
    wavelengths = np.arange(400., 750.)
    
    # Set water quality parameters
    chla = 2.0       # Chlorophyll-a (mg/m3)
    mineral = 0.5    # Suspended minerals (g/m3)
    aCDOM = 0.02     # CDOM absorption at 443nm (1/m)
    
    print(f"\nWater Quality Parameters:")
    print(f"  Chlorophyll-a: {chla} mg/m3")
    print(f"  Mineral: {mineral} g/m3")
    print(f"  CDOM(443nm): {aCDOM} 1/m")
    
    # Compute Rrs
    Rrs = MargModel.runForward(wavelengths, chla, mineral, aCDOM)
    
    # Print some statistics
    print(f"\nRrs Statistics:")
    print(f"  Max Rrs: {np.max(Rrs)*1000:.3f} x10^-3 sr^-1")
    print(f"  Rrs at 443nm: {Rrs[43]*1000:.3f} x10^-3 sr^-1")
    print(f"  Rrs at 555nm: {Rrs[155]*1000:.3f} x10^-3 sr^-1")
    print(f"  Rrs at 670nm: {Rrs[270]*1000:.3f} x10^-3 sr^-1")
    
    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(wavelengths, Rrs * 1000, 'b-', linewidth=2)
    plt.xlabel('Wavelength (nm)', fontsize=12)
    plt.ylabel('Rrs (x10$^{-3}$ sr$^{-1}$)', fontsize=12)
    plt.title('Remote Sensing Reflectance Spectrum', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('example1_single_spectrum.png', dpi=150)
    print(f"\nPlot saved as: example1_single_spectrum.png")
    plt.close()


def example_chlorophyll_sensitivity():
    """
    Example 2: Sensitivity to chlorophyll-a concentration
    """
    print("\n" + "=" * 70)
    print("Example 2: Chlorophyll-a Sensitivity Analysis")
    print("=" * 70)
    
    wavelengths = np.arange(400., 750.)
    
    # Fixed parameters
    mineral = 0.5
    aCDOM = 0.02
    
    # Varying chlorophyll
    chla_values = [0.5, 1.0, 2.0, 5.0, 10.0]
    
    print(f"\nFixed parameters:")
    print(f"  Mineral: {mineral} g/m3")
    print(f"  CDOM(443nm): {aCDOM} 1/m")
    print(f"\nVarying Chlorophyll-a: {chla_values} mg/m3")
    
    plt.figure(figsize=(10, 6))
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    for i, chla in enumerate(chla_values):
        Rrs = MargModel.runForward(wavelengths, chla, mineral, aCDOM)
        plt.plot(wavelengths, Rrs * 1000, linewidth=2, 
                color=colors[i], label=f'Chla = {chla} mg/m³')
        
        # Print peak wavelength and value
        peak_idx = np.argmax(Rrs)
        peak_wl = wavelengths[peak_idx]
        peak_rrs = Rrs[peak_idx] * 1000
        print(f"  Chla = {chla:5.1f} mg/m3: Peak at {peak_wl:.0f} nm, Rrs = {peak_rrs:.3f} x10^-3 sr^-1")
    
    plt.xlabel('Wavelength (nm)', fontsize=12)
    plt.ylabel('Rrs (x10$^{-3}$ sr$^{-1}$)', fontsize=12)
    plt.title('Effect of Chlorophyll-a on Rrs Spectrum', fontsize=14, fontweight='bold')
    plt.legend(fontsize=10, loc='upper right')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('example2_chlorophyll_sensitivity.png', dpi=150)
    print(f"\nPlot saved as: example2_chlorophyll_sensitivity.png")
    plt.close()


def example_mineral_sensitivity():
    """
    Example 3: Sensitivity to mineral concentration
    """
    print("\n" + "=" * 70)
    print("Example 3: Mineral Concentration Sensitivity Analysis")
    print("=" * 70)
    
    wavelengths = np.arange(400., 750.)
    
    # Fixed parameters
    chla = 2.0
    aCDOM = 0.02
    
    # Varying mineral
    mineral_values = [0.1, 0.5, 1.0, 2.0, 5.0]
    
    print(f"\nFixed parameters:")
    print(f"  Chlorophyll-a: {chla} mg/m3")
    print(f"  CDOM(443nm): {aCDOM} 1/m")
    print(f"\nVarying Mineral: {mineral_values} g/m3")
    
    plt.figure(figsize=(10, 6))
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    for i, mineral in enumerate(mineral_values):
        Rrs = MargModel.runForward(wavelengths, chla, mineral, aCDOM)
        plt.plot(wavelengths, Rrs * 1000, linewidth=2,
                color=colors[i], label=f'Mineral = {mineral} g/m³')
        
        # Print Rrs at 555 nm
        rrs_555 = Rrs[155] * 1000
        print(f"  Mineral = {mineral:5.1f} g/m3: Rrs(555nm) = {rrs_555:.3f} x10^-3 sr^-1")
    
    plt.xlabel('Wavelength (nm)', fontsize=12)
    plt.ylabel('Rrs (x10$^{-3}$ sr$^{-1}$)', fontsize=12)
    plt.title('Effect of Mineral Concentration on Rrs Spectrum', fontsize=14, fontweight='bold')
    plt.legend(fontsize=10, loc='upper right')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('example3_mineral_sensitivity.png', dpi=150)
    print(f"\nPlot saved as: example3_mineral_sensitivity.png")
    plt.close()


def example_batch_processing():
    """
    Example 4: Batch processing for multiple stations
    """
    print("\n" + "=" * 70)
    print("Example 4: Batch Processing Multiple Stations")
    print("=" * 70)
    
    wavelengths = np.arange(400., 750.)
    
    # Water quality data for multiple stations
    stations = {
        'Clear_Water': {'chla': 0.5, 'mineral': 0.1, 'aCDOM': 0.005},
        'Coastal_A': {'chla': 2.0, 'mineral': 0.5, 'aCDOM': 0.020},
        'Coastal_B': {'chla': 3.5, 'mineral': 1.2, 'aCDOM': 0.030},
        'Turbid_Water': {'chla': 5.0, 'mineral': 2.0, 'aCDOM': 0.050},
    }
    
    print(f"\nProcessing {len(stations)} stations...")
    print(f"\nStation Data:")
    for name, params in stations.items():
        print(f"  {name:15s}: Chla={params['chla']:.1f}, Mineral={params['mineral']:.1f}, CDOM={params['aCDOM']:.3f}")
    
    # Process all stations
    results = {}
    plt.figure(figsize=(12, 6))
    
    for i, (station_name, params) in enumerate(stations.items()):
        # Compute Rrs
        Rrs = MargModel.runForward(
            wavelengths,
            params['chla'],
            params['mineral'],
            params['aCDOM']
        )
        results[station_name] = Rrs
        
        # Plot
        plt.plot(wavelengths, Rrs * 1000, linewidth=2, label=station_name)
    
    plt.xlabel('Wavelength (nm)', fontsize=12)
    plt.ylabel('Rrs (x10$^{-3}$ sr$^{-1}$)', fontsize=12)
    plt.title('Rrs Spectra for Multiple Stations', fontsize=14, fontweight='bold')
    plt.legend(fontsize=10, loc='upper right')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('example4_batch_processing.png', dpi=150)
    print(f"\nPlot saved as: example4_batch_processing.png")
    plt.close()
    
    # Print summary statistics
    print(f"\nSummary Statistics:")
    print(f"{'Station':<15s} {'Rrs(443nm)':<12s} {'Rrs(555nm)':<12s} {'Rrs(670nm)':<12s} {'Max Rrs':<12s}")
    print("-" * 70)
    for station_name, Rrs in results.items():
        rrs_443 = Rrs[43] * 1000
        rrs_555 = Rrs[155] * 1000
        rrs_670 = Rrs[270] * 1000
        max_rrs = np.max(Rrs) * 1000
        print(f"{station_name:<15s} {rrs_443:<12.3f} {rrs_555:<12.3f} {rrs_670:<12.3f} {max_rrs:<12.3f}")


def example_custom_parameters():
    """
    Example 5: Using custom model parameters (eta, g0, nu)
    """
    print("\n" + "=" * 70)
    print("Example 5: Custom Model Parameters")
    print("=" * 70)
    
    wavelengths = np.arange(400., 750.)
    
    # Water quality
    chla = 2.0
    mineral = 0.5
    aCDOM = 0.02
    
    # Different parameter sets
    param_sets = [
        {'eta': 0.0001, 'g0': 0.001, 'nu': 0.0, 'name': 'Set 1 (eta=0.0001)'},
        {'eta': 0.001, 'g0': 0.001, 'nu': 0.0, 'name': 'Set 2 (eta=0.001)'},
        {'eta': 0.005, 'g0': 0.001, 'nu': 0.0, 'name': 'Set 3 (eta=0.005)'},
    ]
    
    print(f"\nWater Quality: Chla={chla}, Mineral={mineral}, CDOM={aCDOM}")
    print(f"\nTesting different parameter sets:")
    
    plt.figure(figsize=(10, 6))
    
    for params in param_sets:
        Rrs = MargModel.runForward(
            wavelengths, chla, mineral, aCDOM,
            eta_param=params['eta'],
            g0_param=params['g0'],
            nu_param=params['nu']
        )
        
        plt.plot(wavelengths, Rrs * 1000, linewidth=2, label=params['name'])
        
        max_rrs = np.max(Rrs) * 1000
        print(f"  {params['name']}: Max Rrs = {max_rrs:.3f} x10^-3 sr^-1")
    
    plt.xlabel('Wavelength (nm)', fontsize=12)
    plt.ylabel('Rrs (x10$^{-3}$ sr$^{-1}$)', fontsize=12)
    plt.title('Effect of Model Parameters on Rrs', fontsize=14, fontweight='bold')
    plt.legend(fontsize=10, loc='upper right')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('example5_custom_parameters.png', dpi=150)
    print(f"\nPlot saved as: example5_custom_parameters.png")
    plt.close()


if __name__ == "__main__":
    """
    Run all examples
    """
    print("\n" + "=" * 70)
    print("MargModel.py - Usage Examples")
    print("=" * 70)
    print("\nThis script demonstrates various ways to use MargModel.py")
    print("to compute Rrs spectra from water quality parameters.\n")
    
    # Run all examples
    example_single_spectrum()
    example_chlorophyll_sensitivity()
    example_mineral_sensitivity()
    example_batch_processing()
    example_custom_parameters()
    
    print("\n" + "=" * 70)
    print("All examples completed successfully!")
    print("=" * 70)
    print("\nGenerated files:")
    print("  - example1_single_spectrum.png")
    print("  - example2_chlorophyll_sensitivity.png")
    print("  - example3_mineral_sensitivity.png")
    print("  - example4_batch_processing.png")
    print("  - example5_custom_parameters.png")
    print("\n")

