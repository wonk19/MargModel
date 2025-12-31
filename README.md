# MargModel.py - Standalone CCRR-RT Forward Model

## Overview

`MargModel.py` is a standalone Python module that implements the CCRR-RT (Case 2 Regional CoastColour Radiative Transfer) forward model. This model computes remote sensing reflectance (Rrs) spectra from water quality parameters including chlorophyll-a, suspended minerals, and colored dissolved organic matter (CDOM).

This module was extracted from the K03 optimization pipeline and can be distributed independently to third parties who need to generate Rrs spectra from water quality parameters.

## Requirements

- Python 3.6+
- numpy
- scipy

Install dependencies:
```bash
pip install numpy scipy
```

## File Structure

```
MargModel.py              # Main module
independent/              # Required data files
├── AB_box_Margalef.pkl   # Pigment absorption parameters
└── Water_a_b.txt         # Water optical properties
```

**Important:** The `independent/` folder and its contents must be in the same directory as `MargModel.py`.

## Quick Start

### Basic Usage

```python
import numpy as np
import MargModel

# Define wavelength range (400-750 nm recommended)
wavelengths = np.arange(400., 750.)

# Set water quality parameters
chla = 2.0       # Chlorophyll-a (mg/m3)
mineral = 0.5    # Suspended minerals (g/m3)
aCDOM = 0.02     # CDOM absorption at 443nm (1/m)

# Run forward model (includes fluorescence by default)
Rrs = MargModel.runForward(wavelengths, chla, mineral, aCDOM)

# To exclude fluorescence:
Rrs_no_fl = MargModel.runForward(wavelengths, chla, mineral, aCDOM, 
                                  add_fluorescence_flag=False)

# Rrs is now a numpy array with Rrs values (sr^-1) for each wavelength
```

### Plotting Results

```python
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
plt.plot(wavelengths, Rrs * 1000)  # Convert to 10^-3 sr^-1
plt.xlabel('Wavelength (nm)')
plt.ylabel('Rrs (×10⁻³ sr⁻¹)')
plt.title('Remote Sensing Reflectance Spectrum')
plt.grid(True, alpha=0.3)
plt.show()
```

## API Reference

### Main Function: `runForward()`

```python
Rrs = MargModel.runForward(wavelengths, chla, mineral, aCDOM, 
                           eta_param=0.7, g0_param=0.0002, nu_param=0.0,
                           fresnel_add=0.02, rho_sky=None, sol_zen=30.0)
```

**Parameters:**

- **wavelengths** : array-like
  - Wavelengths for which to compute Rrs (nm)
  - Recommended range: 400-750 nm
  - Example: `np.arange(400., 750.)`

- **chla** : float
  - Chlorophyll-a concentration (mg/m³)
  - Typical range: 0.01 - 100 mg/m³
  - Controls pigment absorption and backscattering

- **mineral** : float
  - Mineral (suspended sediment) concentration (g/m³)
  - Typical range: 0.1 - 10 g/m³
  - Controls mineral absorption and backscattering

- **aCDOM** : float
  - CDOM absorption coefficient at 443nm (1/m)
  - Typical range: 0.001 - 0.5 (1/m)
  - Controls yellow substance absorption

- **eta_param** : float, optional (default: 0.7)
  - Pigment backscattering exponent
  - Typical range: 0.0001 - 0.05
  - Optimized value from K06 analysis: 0.7

- **g0_param** : float, optional (default: 0.0002)
  - Pigment backscattering coefficient
  - Typical range: 0.0001 - 0.01
  - Optimized value from K06 analysis: 0.0002

- **nu_param** : float, optional (default: 0.0)
  - Pigment backscattering spectral slope
  - Typical range: 0.0 - 0.1
  - Optimized value from K06 analysis: 0.0

- **fresnel_add** : float, optional (default: 0.02)
  - Fresnel reflection correction term
  - Accounts for surface reflection effects

- **rho_sky** : array-like or None, optional (default: None)
  - Sky reflectance spectrum
  - If None, assumed to be zero (no sky contribution)
  - Should have same length as wavelengths

- **sol_zen** : float, optional (default: 30.0)
  - Solar zenith angle in degrees
  - Note: Currently fixed at 30° internally in the BRDF calculation

- **add_fluorescence_flag** : bool, optional (default: True)
  - Whether to add chlorophyll fluorescence to Rrs
  - Fluorescence parameterization based on K06 field analysis
  - Recommended to keep True for realistic simulation

- **fluorescence_method** : str, optional (default: 'semilog')
  - Method for fluorescence calculation
  - 'semilog': Peak wavelength and intensity scale with ln(Chla) (recommended)
  - 'linear': Peak wavelength and intensity scale linearly with Chla

**Returns:**

- **Rrs** : ndarray
  - Remote sensing reflectance (sr⁻¹) at each wavelength
  - Same length as input wavelengths array

### Utility Functions

#### Get Default Parameters

```python
params = MargModel.getDefaultParameters()
# Returns: {'eta': 0.7, 'g0': 0.0002, 'nu': 0.0}
```

Get default model parameters. These are the default values used if not specified in `runForward()`.

## Example Use Cases

### Example 1: Generate Rrs for Different Chlorophyll Levels

```python
import numpy as np
import matplotlib.pyplot as plt
import MargModel

wavelengths = np.arange(400., 750.)

# Fixed mineral and CDOM
mineral = 0.5
aCDOM = 0.02

# Varying chlorophyll
chla_values = [0.5, 1.0, 2.0, 5.0, 10.0]

plt.figure(figsize=(10, 6))
for chla in chla_values:
    Rrs = MargModel.runForward(wavelengths, chla, mineral, aCDOM)
    plt.plot(wavelengths, Rrs * 1000, linewidth=2, label=f'Chla = {chla} mg/m³')

plt.xlabel('Wavelength (nm)', fontsize=12)
plt.ylabel('Rrs (×10⁻³ sr⁻¹)', fontsize=12)
plt.title('Effect of Chlorophyll-a on Rrs Spectrum')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
```

### Example 2: Sensitivity Analysis for Mineral Concentration

```python
import numpy as np
import matplotlib.pyplot as plt
import MargModel

wavelengths = np.arange(400., 750.)

# Fixed chlorophyll and CDOM
chla = 2.0
aCDOM = 0.02

# Varying mineral concentration
mineral_values = [0.1, 0.5, 1.0, 2.0, 5.0]

plt.figure(figsize=(10, 6))
for mineral in mineral_values:
    Rrs = MargModel.runForward(wavelengths, chla, mineral, aCDOM)
    plt.plot(wavelengths, Rrs * 1000, linewidth=2, label=f'Mineral = {mineral} g/m³')

plt.xlabel('Wavelength (nm)', fontsize=12)
plt.ylabel('Rrs (×10⁻³ sr⁻¹)', fontsize=12)
plt.title('Effect of Mineral Concentration on Rrs Spectrum')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
```

### Example 3: Batch Processing Multiple Stations

```python
import numpy as np
import pandas as pd
import MargModel

wavelengths = np.arange(400., 750.)

# Water quality data for multiple stations
stations_data = {
    'Station_A': {'chla': 1.5, 'mineral': 0.3, 'aCDOM': 0.015},
    'Station_B': {'chla': 3.2, 'mineral': 0.8, 'aCDOM': 0.025},
    'Station_C': {'chla': 0.8, 'mineral': 0.2, 'aCDOM': 0.010},
}

# Process all stations
results = {}
for station_name, params in stations_data.items():
    Rrs = MargModel.runForward(
        wavelengths, 
        params['chla'], 
        params['mineral'], 
        params['aCDOM']
    )
    results[station_name] = Rrs
    
    # Calculate Rrs at specific bands
    idx_443 = np.argmin(np.abs(wavelengths - 443))
    idx_555 = np.argmin(np.abs(wavelengths - 555))
    
    print(f"{station_name}:")
    print(f"  Rrs(443nm) = {Rrs[idx_443]*1000:.3f} ×10⁻³ sr⁻¹")
    print(f"  Rrs(555nm) = {Rrs[idx_555]*1000:.3f} ×10⁻³ sr⁻¹")
```

### Example 4: Using Custom Model Parameters

```python
import numpy as np
import MargModel

wavelengths = np.arange(400., 750.)

# Set custom model parameters (from optimization results)
eta = 0.005
g0 = 0.001
nu = 0.0

# Water quality parameters
chla = 2.0
mineral = 0.5
aCDOM = 0.02

# Run with custom parameters
Rrs = MargModel.runForward(
    wavelengths, chla, mineral, aCDOM,
    eta_param=eta, g0_param=g0, nu_param=nu
)

# For multiple runs with same parameters, just specify them each time
Rrs1 = MargModel.runForward(wavelengths, 1.0, 0.3, 0.01, 
                            eta_param=eta, g0_param=g0, nu_param=nu)
Rrs2 = MargModel.runForward(wavelengths, 2.0, 0.5, 0.02, 
                            eta_param=eta, g0_param=g0, nu_param=nu)
Rrs3 = MargModel.runForward(wavelengths, 3.0, 0.8, 0.03, 
                            eta_param=eta, g0_param=g0, nu_param=nu)
```

## Model Details

### CCRR-RT Model Components

The CCRR-RT forward model computes Rrs based on inherent optical properties (IOPs):

1. **Absorption Components:**
   - Pigment absorption (from Margalef-type phytoplankton)
   - Mineral absorption (Babin et al., 2003)
   - CDOM absorption (exponential decay from 443nm)
   - Pure water absorption

2. **Backscattering Components:**
   - Pigment backscattering (Morel 2002 parameterization)
   - Mineral backscattering (based on KORUS data)
   - Pure water backscattering

3. **Bidirectional Reflectance:**
   - Morel's f/Q lookup tables (Morel 2002)
   - Accounts for solar zenith angle effects

### Key References

- Morel, A., & Gentili, B. (2002). A simple band ratio technique to quantify the colored dissolved and detrital organic material from ocean color remotely sensed data. Remote Sensing of Environment, 113(5), 998-1011.

- Babin, M., Stramski, D., Ferrari, G. M., Claustre, H., Bricaud, A., Obolensky, G., & Hoepffner, N. (2003). Variations in the light absorption coefficients of phytoplankton, nonalgal particles, and dissolved organic matter in coastal waters around Europe. Journal of Geophysical Research: Oceans, 108(C7).

## Testing

Run the built-in test suite:

```bash
python MargModel.py
```

This will:
1. Run three test cases (coastal, clear, turbid water)
2. Generate a comparison plot saved as `MargModel_test_output.png`
3. Display results in the console

Expected output:
```
Test Case 1: Typical coastal water
  Chla: 2.0 mg/m3, Mineral: 0.5 g/m3, aCDOM(443): 0.02 (1/m)
  Max Rrs: 7.762 x10^-3 sr^-1
  
Test Case 2: Clear water
  Chla: 0.5 mg/m3, Mineral: 0.1 g/m3, aCDOM(443): 0.005 (1/m)
  Max Rrs: 11.699 x10^-3 sr^-1
  
Test Case 3: Turbid water
  Chla: 5.0 mg/m3, Mineral: 2.0 g/m3, aCDOM(443): 0.05 (1/m)
  Max Rrs: 10.904 x10^-3 sr^-1
```

## Distribution

To share this model with third parties, provide:

1. `MargModel.py` - The main Python module
2. `independent/` folder with:
   - `AB_box_Margalef.pkl`
   - `Water_a_b.txt`
3. `MargModel_README.md` - This documentation file

The recipient can then:
1. Place all files in the same directory
2. Install numpy and scipy
3. Import and use `MargModel` in their Python scripts

## Performance Notes

The module uses several caching mechanisms to improve performance:
- Water optical properties are loaded once and cached
- Absorption parameters are loaded once and cached
- Morel interpolators are created once and reused
- Wavelength-dependent exponential terms are cached

For batch processing of many spectra, the first call may be slower due to initialization, but subsequent calls will be much faster.

## Troubleshooting

### Error: "FileNotFoundError: independent/AB_box_Margalef.pkl"

**Solution:** Ensure the `independent/` folder and its contents are in the same directory as `MargModel.py`.

### Error: "ImportError: No module named scipy"

**Solution:** Install scipy: `pip install scipy`

### Warning: "Wavelength outside valid range"

**Solution:** Use wavelengths in the range 400-750 nm for best results. The model will clip wavelengths outside the range of the lookup tables (412.5-660 nm for Morel parameters).

### Unexpected Rrs values (too high or too low)

**Possible causes:**
1. Check parameter units:
   - chla in mg/m³ (not μg/L)
   - mineral in g/m³
   - aCDOM in 1/m
2. Verify parameter ranges are realistic
3. Ensure wavelengths are in nanometers (not micrometers)

## Contact

For questions or issues related to this model, please contact the author or refer to the original K03 optimization pipeline documentation.

## Chlorophyll-a Correction

### Overview

When using the forward model to estimate chlorophyll-a from observed Rrs spectra (inverse problem), the raw estimates may show systematic bias. A 3rd-order polynomial correction in log-space can significantly improve accuracy.

### Correction Function

```python
chla_corrected = MargModel.apply_chla_correction(chla_est, poly_coeffs)
```

The correction is applied as:
```
log(chla_corrected) = a3*log(chla_est)³ + a2*log(chla_est)² + a1*log(chla_est) + a0
```

### Recommended Polynomial Coefficients

The polynomial coefficients are parameter-set specific. Below are the top 10 parameter sets based on two criteria:

#### Top 10 by R² (Best Chla Estimation Accuracy)

| Rank | eta   | g0      | nu   | R²     | RMSE_Rrs | poly_a3   | poly_a2   | poly_a1  | poly_a0   |
|------|-------|---------|------|--------|----------|-----------|-----------|----------|-----------|
| 1    | 0.500 | 0.00100 | 1.00 | 0.9390 | 0.175120 | 0.101620  | -0.322535 | 1.073302 | -0.221764 |
| 2    | 0.775 | 0.00020 | 1.00 | 0.9334 | 0.165825 | 0.085890  | -0.210978 | 0.866397 | -0.164719 |
| 3    | 0.950 | 0.00005 | 1.00 | 0.9332 | 0.145160 | 0.010728  | 0.137408  | 0.352919 | 0.034673  |
| 4    | 0.650 | 0.00040 | 1.00 | 0.9312 | 0.171452 | 0.083595  | -0.222705 | 0.921967 | -0.201947 |
| 5    | 0.625 | 0.00040 | 1.00 | 0.9310 | 0.174814 | 0.056867  | -0.099050 | 0.747213 | -0.139404 |
| 6    | 0.850 | 0.00010 | 1.00 | 0.9302 | 0.154771 | 0.012258  | 0.140254  | 0.334053 | 0.061684  |
| 7    | 0.600 | 0.00060 | 1.00 | 0.9286 | 0.180918 | 0.047271  | -0.084712 | 0.784730 | -0.150628 |
| 8    | 0.575 | 0.00060 | 1.00 | 0.9267 | 0.179130 | 0.036921  | -0.030353 | 0.684607 | -0.103030 |
| 9    | 0.675 | 0.00040 | 0.75 | 0.9252 | 0.181476 | 0.066129  | -0.144183 | 0.815283 | -0.164588 |
| 10   | 0.675 | 0.00040 | 1.00 | 0.9246 | 0.175862 | 0.075362  | -0.181135 | 0.868817 | -0.173960 |

#### Top 10 by RMSE(Rrs) (Best Spectral Fitting)

| Rank | eta   | g0      | nu   | R²     | RMSE_Rrs | poly_a3   | poly_a2  | poly_a1  | poly_a0  |
|------|-------|---------|------|--------|----------|-----------|----------|----------|----------|
| 1    | 0.625 | 0.00010 | 0.25 | 0.8317 | 0.138459 | -0.023557 | 0.239656 | 0.241357 | 0.074092 |
| 2    | 0.600 | 0.00010 | 0.00 | 0.8657 | 0.138800 | -0.016444 | 0.229737 | 0.223725 | 0.090080 |
| 3    | 0.750 | 0.00005 | 0.25 | 0.8697 | 0.138962 | 0.006367  | 0.121789 | 0.381765 | 0.024732 |
| 4    | 0.775 | 0.00005 | 0.75 | 0.8764 | 0.140155 | 0.006764  | 0.124196 | 0.383149 | 0.021258 |
| 5    | 0.600 | 0.00010 | 0.25 | 0.8681 | 0.140367 | -0.016027 | 0.230239 | 0.222851 | 0.089623 |
| 6    | 0.575 | 0.00010 | 0.00 | 0.8646 | 0.140465 | -0.017433 | 0.232398 | 0.225316 | 0.084033 |
| 7    | 0.850 | 0.00005 | 0.00 | 0.8680 | 0.140480 | 0.004128  | 0.115731 | 0.404868 | 0.012019 |
| 8    | 0.850 | 0.00005 | 0.75 | 0.8775 | 0.140969 | 0.005941  | 0.122976 | 0.396324 | 0.014580 |
| 9    | 0.625 | 0.00010 | 0.75 | 0.8723 | 0.141006 | -0.016906 | 0.233909 | 0.227529 | 0.085500 |
| 10   | 0.900 | 0.00005 | 0.50 | 0.8769 | 0.141029 | 0.006582  | 0.124217 | 0.387385 | 0.019099 |

### Usage Example

```python
import numpy as np
import MargModel

# Use the top R² parameter set (Rank 1)
eta = 0.500
g0 = 0.00100
nu = 1.00
poly_coeffs = [0.101620, -0.322535, 1.073302, -0.221764]

# Simulate forward model
wavelengths = np.arange(400., 750.)
chla_true = 5.0
mineral = 0.5
aCDOM = 0.02

Rrs = MargModel.runForward(wavelengths, chla_true, mineral, aCDOM,
                           eta=eta, g0=g0, nu=nu)

# In inverse problem, you would estimate chla from Rrs
# For demonstration, assume we got a raw estimate
chla_est = 6.5  # Raw estimate (biased)

# Apply correction
chla_corrected = MargModel.apply_chla_correction(chla_est, poly_coeffs)

print(f"True Chla: {chla_true:.2f} mg/m³")
print(f"Raw estimate: {chla_est:.2f} mg/m³")
print(f"Corrected estimate: {chla_corrected:.2f} mg/m³")
```

### Notes on Parameter Selection

- **For Chla estimation**: Use parameter sets from "Top 10 by R²" table
- **For spectral fitting**: Use parameter sets from "Top 10 by RMSE(Rrs)" table
- The correction coefficients are derived from field validation data (81 stations)
- Correction is most effective for Chla range: 0.5-50 mg/m³

## Version History

- **v1.1** (2024-12-31) - Added polynomial correction
  - Added `apply_chla_correction()` function
  - Included calibrated polynomial coefficients for top 10 parameter sets
  - Updated documentation with correction usage examples

- **v1.0** (2024) - Initial standalone release
  - Extracted from K03_optimize_CCRR_RT.py
  - Self-contained with independent data files
  - Comprehensive API and documentation

