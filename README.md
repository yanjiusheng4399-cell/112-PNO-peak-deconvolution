# PrNiO2 Sigma RIXS Peak Deconvolution

MATLAB workflow for reproducing Sigma-polarized PrNiO2 RIXS peak
deconvolution and dispersion-style plots from the provided BNL data.

## Contents

- `make_PNO_sigma_RIXS_figures.m`: main MATLAB script.
- `PNO_STO_Sigma_BNL.csv`: PrNiO2/STO Sigma BNL input spectra.
- `PNO_LSAT_Sigma_BNL.csv`: PrNiO2/LSAT Sigma BNL input spectra.
- `PNO_sigma_epsilon_gamma.txt`: fitted epsilon_q and gamma_q table.
- `PNO_sigma_component_diagnostics.txt`: fitted component diagnostics and
  component hierarchy checks.
- `PNO_sigma_fig4a_point_compare.txt`: comparison against the Fig. 4a branch
  prior used for point selection.
- `Fig2_sigma_fits.png/.fig`: peak-deconvolution panels.
- `Fig3_sigma_magnetic_map.png/.fig`: nonmagnetic-subtracted magnetic map.
- `Fig4_sigma_epsilon_gamma_scatter.png/.fig`: epsilon_q and gamma_q scatter
  plot with error bars.

## Run

Open MATLAB in this folder and run:

```matlab
make_PNO_sigma_RIXS_figures
```

The script reads only the two Sigma BNL CSV files and regenerates the result
tables and three figure outputs.
