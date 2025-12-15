# MargModel Distribution Package

## Contents

This package contains a standalone implementation of the CCRR-RT forward model for computing remote sensing reflectance (Rrs) from water quality parameters.

### Files to Distribute

```
MargModel.py                   # Main module
MargModel_README.md            # Complete documentation
MargModel_example.py           # Usage examples
independent/                   # Required data files (must be included)
├── AB_box_Margalef.pkl        # Pigment absorption parameters
└── Water_a_b.txt              # Water optical properties
```

### Optional Files (for reference/testing)

```
MargModel_verification.py      # Verification against original implementation
MargModel_DISTRIBUTION.md      # This file
```

## Quick Installation Guide for Recipients

1. **Receive the package** - Get all files listed above

2. **Ensure directory structure** - Keep files organized as:
   ```
   your_directory/
   ├── MargModel.py
   ├── MargModel_README.md
   ├── MargModel_example.py
   └── independent/
       ├── AB_box_Margalef.pkl
       └── Water_a_b.txt
   ```

3. **Install Python dependencies**:
   ```bash
   pip install numpy scipy matplotlib
   ```

4. **Test the installation**:
   ```bash
   python MargModel.py
   ```
   
   This should generate `MargModel_test_output.png` and display test results.

5. **Run examples**:
   ```bash
   python MargModel_example.py
   ```
   
   This generates 5 example plots demonstrating various uses.

## Basic Usage

```python
import numpy as np
import MargModel

# Define wavelengths
wavelengths = np.arange(400., 750.)

# Set water quality parameters
chla = 2.0       # Chlorophyll-a (mg/m3)
mineral = 0.5    # Suspended minerals (g/m3)
aCDOM = 0.02     # CDOM absorption at 443nm (1/m)

# Compute Rrs
Rrs = MargModel.runForward(wavelengths, chla, mineral, aCDOM)

# Rrs is now a numpy array with Rrs values (sr^-1)
```

See `MargModel_README.md` for complete documentation.

## Model Parameters

The forward model uses three tunable parameters that affect pigment backscattering:

- **eta** (default: 0.001) - Pigment backscattering exponent
- **g0** (default: 0.001) - Pigment backscattering coefficient  
- **nu** (default: 0.0) - Pigment backscattering spectral slope

These defaults are from optimization results. To use different values:

```python
Rrs = MargModel.runForward(
    wavelengths, chla, mineral, aCDOM,
    eta_param=0.005, g0_param=0.001, nu_param=0.0
)
```

## Input Parameter Ranges

| Parameter | Unit | Typical Range | Description |
|-----------|------|---------------|-------------|
| wavelengths | nm | 400-750 | Wavelengths for Rrs computation |
| chla | mg/m³ | 0.01-100 | Chlorophyll-a concentration |
| mineral | g/m³ | 0.1-10 | Suspended mineral concentration |
| aCDOM | 1/m | 0.001-0.5 | CDOM absorption at 443nm |
| eta_param | - | 0.0001-0.05 | Pigment backscattering exponent |
| g0_param | - | 0.0001-0.01 | Pigment backscattering coefficient |
| nu_param | - | 0.0-0.1 | Pigment backscattering spectral slope |

## Output

The `runForward()` function returns:
- **Rrs** : numpy array of remote sensing reflectance (sr⁻¹)
- Same length as input wavelengths
- Typical values: 0.001 - 0.02 sr⁻¹ (or 1-20 ×10⁻³ sr⁻¹)

## Verification

The module has been verified to produce **identical results** to the original K03 implementation:
- Maximum difference < 1×10⁻¹⁰ sr⁻¹ (numerical precision limit)
- Verified across multiple water quality conditions
- Verified across different parameter settings

Run `MargModel_verification.py` to confirm (requires original K03 codebase).

## Support

For questions or issues:
1. Check `MargModel_README.md` for detailed documentation
2. Review `MargModel_example.py` for usage examples
3. Contact the original author or refer to K03 documentation

## Citation

If you use this model in your research, please cite:

[Add appropriate citation information here]

## Model Description

This is the CCRR-RT (Case 2 Regional CoastColour Radiative Transfer) forward model, which computes remote sensing reflectance from:

1. **Water Quality Inputs:**
   - Chlorophyll-a concentration
   - Suspended mineral concentration
   - CDOM absorption

2. **Physical Parameterizations:**
   - Margalef-type pigment absorption (from in-situ measurements)
   - Mineral absorption and scattering (Babin et al., 2003)
   - CDOM exponential absorption
   - Pure water absorption and scattering
   - Morel's f/Q bidirectional effects (Morel, 2002)

3. **Output:**
   - Remote sensing reflectance (Rrs) spectrum

The model is optimized for coastal/Case-2 waters with moderate to high turbidity.

## License

[Add license information if applicable]

## Version

Version 1.0 (2024)

Extracted from K03_optimize_CCRR_RT.py pipeline.

