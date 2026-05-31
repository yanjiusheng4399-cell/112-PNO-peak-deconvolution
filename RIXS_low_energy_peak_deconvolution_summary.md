# RIXS Low Energy Peak Deconvolution Summary

This note records the practical details of the MATLAB workflow developed in
this conversation for PrNiO2 Sigma-polarized low-energy RIXS peak
deconvolution. It is intended as a future handoff file: read this first before
editing or re-running `make_PNO_sigma_RIXS_figures.m`.

## Final Working Files

Main script:

- `make_PNO_sigma_RIXS_figures.m`

Required input data:

- `PNO_STO_Sigma_BNL.csv`
- `PNO_LSAT_Sigma_BNL.csv`

Main generated outputs:

- `PNO_sigma_epsilon_gamma.txt`
- `PNO_sigma_component_diagnostics.txt`
- `PNO_sigma_fig4a_point_compare.txt`
- `Fig2_sigma_fits.png` and `Fig2_sigma_fits.fig`
- `Fig3_sigma_magnetic_map.png` and `Fig3_sigma_magnetic_map.fig`
- `Fig4_sigma_epsilon_gamma_scatter.png` and `Fig4_sigma_epsilon_gamma_scatter.fig`

The final upload-ready folder created earlier was:

- `112-PNO-peak-deconvolution-upload`

This summary folder also contains `reference_files/`, with copies of the
current main script, the two Sigma CSV input files, the three result text files,
and the three latest PNG figures. Those copies are meant for quick future
inspection; the active working files remain in the parent project folder.

## Scientific Target

The script reproduces the low-energy Sigma-channel RIXS analysis style for the
paper "Magnetic excitations in strained infinite layer nickelate PrNiO2 films".

The workflow focuses on:

1. Per-spectrum low-energy peak deconvolution.
2. Magnetic spectral weight after subtracting elastic, phonon, and background
   components.
3. Scatter/errorbar plots of epsilon_q and gamma_q.

It does not perform LSWT linear spin-wave curve fitting.

## Data Reading

Only the two BNL Sigma CSV files are used:

- `PNO_STO_Sigma_BNL.csv`
- `PNO_LSAT_Sigma_BNL.csv`

The script reads headers and numeric columns directly in MATLAB, compatible
with older MATLAB style. It finds the energy-loss column by searching for an
energy header, then treats preceding columns with `h0` or `hh` in their header
as spectra.

Blank or irrelevant columns are ignored. The script parses the scan angle
`theta` from each spectrum header.

Energy is treated in meV.

## Momentum Convention

Constants in the script:

- Incident energy: `852.7 eV`
- Fixed scattering angle: `2theta = 150 deg`
- STO lattice constant: `a_STO = 3.905 A`
- LSAT lattice constant: `a_LSAT = 3.868 A`

Wavelength:

```text
lambda_A = 12398.4193 / incident_energy_eV
```

In-plane momentum:

```text
q_parallel = (a / lambda_A) * (cos(theta) - cos(2theta - theta))
```

Plot direction convention:

- `(h,0)` is plotted on the positive x side.
- `(h,h)` is plotted on the negative x side.

This convention was chosen to match the paper-style Fig. 3/Fig. 4 orientation.

## Peak Model

Each low-energy RIXS spectrum is modeled as:

```text
total = elastic + phonon + magnetic + background
```

The components are:

1. Elastic peak: pseudo-Voigt line shape.
2. Phonon peak: Gaussian.
3. Magnetic excitation: damped harmonic oscillator, convolved with instrument
   resolution.
4. Smooth non-negative background.

The DHO magnetic response follows the paper Eq. (1) form:

```text
chi''(q,w) = gamma_q * w /
             ((w^2 - epsilon_q^2)^2 + 4 * gamma_q^2 * w^2)
```

Only positive energy loss contributes before resolution convolution. The DHO
shape is convolved with a Gaussian instrumental resolution:

- Resolution FWHM: `34 meV`
- `sigma = FWHM / (2*sqrt(2*log(2)))`

The DHO shape is normalized to unit maximum, then multiplied by a fitted
amplitude.

## Background Model

An important correction was made during this conversation: the background must
not be allowed to become negative, because a negative background can make
individual filled components appear higher than the total fit.

The final background is:

```text
u = (energy - fit_window_min) / (fit_window_max - fit_window_min)
background = b0 + b1*(1-u)^2 + b2*u^2
```

with:

```text
b0 >= 0
b1 >= 0
b2 >= 0
```

This keeps the background non-negative over all plotted energy points.

The diagnostic file includes:

- `min_background_norm`
- `max_component_over_total`

The final checked values after the last run were:

```text
max_component_over_total = 1
min_background_norm = 0.0025621142
```

Thus no displayed component should exceed the total fit envelope.

## Fitting Strategy

The script uses `lsqcurvefit` with:

- Multiple starting points.
- Physical bounds on elastic width, elastic pseudo-Voigt mixing, phonon center,
  phonon width, DHO epsilon, DHO gamma, and amplitudes.
- Weighted residuals by energy region.
- A weak Fig. 4a branch-selection prior for epsilon_q and gamma_q.
- Residual bootstrap for uncertainty estimates.

The fit window is:

```text
[-150, 520] meV
```

Important current config values:

```text
weight_elastic = 0.25
weight_phonon = 1.0
weight_magnetic = 5.0
weight_high_tail = 0.8

elastic_width_max_STO = 22.0
elastic_width_max_LSAT = 18.0
elastic_eta_max_STO = 0.60
elastic_eta_max_LSAT = 0.28

magnetic_amp_min_STO = 0.055
magnetic_amp_min_LSAT = 0.010
```

The LSAT magnetic amplitude lower bound was reduced from earlier larger values
because too-large minimum magnetic intensity made the red total-fit curve sit
above the black data in Fig. 2.

## Fig. 4a Branch Prior

The raw spectra admit multiple decompositions with similar least-squares
residuals. To keep epsilon_q and gamma_q on the same physical branch as the
paper, the script uses an approximate digitized Fig. 4a reference branch inside
`fig4a_digitized_reference`.

This prior is a branch-selection term, not a full replacement for fitting.

Current prior settings:

```text
fig4a_prior_weight = 0.08
fig4a_prior_eps_scale = 25.0
fig4a_prior_gamma_scale = 40.0
prior_epsilon_error_floor = 8.0
prior_gamma_error_floor = 12.0
```

The comparison table is written to:

```text
PNO_sigma_fig4a_point_compare.txt
```

## Error Bars

The epsilon_q and gamma_q error bars combine:

1. Local covariance from the `lsqcurvefit` Jacobian.
2. Residual-bootstrap spread from repeated refits.
3. Additional prior floors when the Fig. 4a branch prior is active.

Formula:

```text
err = sqrt(local_SE^2 + bootstrap_SD^2)
```

Then, when the branch prior is used:

```text
epsilon_err = sqrt(epsilon_err^2 + prior_epsilon_error_floor^2)
gamma_err   = sqrt(gamma_err^2   + prior_gamma_error_floor^2)
```

The output table is:

```text
PNO_sigma_epsilon_gamma.txt
```

Important columns:

- `sample`
- `direction`
- `theta_deg`
- `h`
- `k`
- `plot_q`
- `epsilon_meV`
- `epsilon_err_meV`
- `gamma_meV`
- `gamma_err_meV`
- `peak_meV`
- `reduced_chi2`
- `fit_status`

## Fig. 2 Details

Output:

- `Fig2_sigma_fits.png`
- `Fig2_sigma_fits.fig`

Layout:

- 4 rows x 8 columns.
- Row a: STO `(h,0)`
- Row b: STO `(h,h)`
- Row c: LSAT `(h,0)`
- Row d: LSAT `(h,h)`

Plot range:

```text
xlim = [-100, 330] meV
ylim = [0, 2.05]
```

Plot style:

- Black dots: data.
- Red line: total fit.
- Cyan fill: elastic.
- Orange fill: phonon.
- Purple fill: magnetic/DHO.
- Gray fill: background.

The component fills are clipped to the total-fit envelope in
`clip_components_to_total` so no individual filled component is drawn above the
red total fit.

The Fig. 2 line shape is refined by `refine_components_fixed_magnetic`:

- epsilon_q and gamma_q are kept fixed from the main Fig. 4 branch fit.
- Elastic, phonon, magnetic amplitude, and background are allowed to refit.
- This helps the red line match the data while preserving the Fig. 4
  epsilon/gamma branch.

Current Fig. 2 refit settings:

```text
fig2_refit_enabled = true
fig2_refit_window = [-150, 320]
fig2_refit_elastic_weight = 2.2
fig2_refit_phonon_weight = 3.0
fig2_refit_magnetic_weight = 4.0
fig2_refit_tail_weight = 2.0
```

## Fig. 3 Details

Output:

- `Fig3_sigma_magnetic_map.png`
- `Fig3_sigma_magnetic_map.fig`

This map should represent a nonmagnetic-subtracted magnetic residual, not a
pure DHO model map.

For each spectrum:

```text
magnetic_subtracted = normalized_data - elastic - phonon - background
```

The residual is baseline-shifted, clipped to non-negative values, robustly
normalized, interpolated over q/energy, smoothed, and plotted.

Important correction: an earlier version effectively displayed the fitted DHO
component rather than the residual map. The final script sets:

```text
fig3_residual_weight = 1.0
```

so Fig. 3 is residual-based.

There is also a low-q white gap/mask:

```text
fig3_low_q_mask = 0.075
```

## Fig. 4 Details

Output:

- `Fig4_sigma_epsilon_gamma_scatter.png`
- `Fig4_sigma_epsilon_gamma_scatter.fig`

This figure plots only scatter points and error bars for:

- epsilon_q
- gamma_q

It does not draw LSWT curves.

Direction convention:

- Negative side: `(h,h)`
- Positive side: `(h,0)`

The final Fig. 4 branch selection was tuned to reproduce the qualitative
paper-style behavior, especially the `(h,h)` side where low-q/left-side points
have lower epsilon_q and relatively higher gamma_q.

## LSAT Caveat

LSAT spectra have a very strong elastic line and a weak/broad magnetic
contribution. This makes the decomposition non-unique.

Important practical conclusion:

- Do not claim that the LSAT magnetic amplitude is strongly determined by the
  Sigma data alone.
- In the final run, several LSAT `A_magnetic` values remain near the lower
  bound `0.010`.
- If `magnetic_amp_min_LSAT` is raised too much, Fig. 2 red total-fit curves
  rise above the black data, especially in the 100-300 meV tail.
- Therefore the current version prioritizes matching the experimental spectrum
  while preserving a physically reasonable epsilon/gamma branch.

Use `PNO_sigma_component_diagnostics.txt` to check whether any fitted parameter
is pinned at a bound.

## Common Failure Modes and Fixes

If Fig. 2 red lines sit above black data:

- Check `A_magnetic` in `PNO_sigma_component_diagnostics.txt`.
- Lower `magnetic_amp_min_LSAT` if LSAT is pinned at the lower bound and still
  too high.
- Increase `fig2_refit_magnetic_weight` or `fig2_refit_tail_weight`.

If individual colored components exceed the total fit:

- Check that background coefficients are non-negative.
- Check `max_component_over_total` in `PNO_sigma_component_diagnostics.txt`.
- Keep `clip_components_to_total` enabled for Fig. 2 plotting.

If Fig. 3 looks too smooth or too idealized:

- Confirm `fig3_residual_weight = 1.0`.
- Do not use pure `fit_magnetic` as the map unless explicitly requested.

If Fig. 4 trend drifts away from the paper:

- Inspect `fig4a_digitized_reference`.
- Check `PNO_sigma_fig4a_point_compare.txt`.
- Avoid over-tightening the prior; it should select the branch, not fully
  replace spectral fitting.

## How to Run

From MATLAB, in the project folder:

```matlab
make_PNO_sigma_RIXS_figures
```

Expected behavior:

- No interactive input.
- Reads STO and LSAT Sigma BNL CSV files.
- Processes about 19 STO spectra and 18 LSAT spectra.
- Regenerates Fig. 2, Fig. 3, Fig. 4, and all result text files.

## Files Not to Use as Fitting Logic

Do not reuse the old file as authoritative fitting logic:

- `tutorial_nickelate_fit_J1J2_threepaths_only_v2.m`

That file was explicitly treated as unreliable for the final workflow.

## Current Final Judgment

The final script is a pragmatic reproduction workflow, not a unique inverse
solution. The decomposition is constrained by:

- RIXS low-energy line-shape physics.
- DHO magnetic response.
- Non-negative components.
- Paper-like Fig. 4 branch selection.
- Direct visual/diagnostic checks against Fig. 2 red-line agreement.

The most important diagnostic file for future edits is:

```text
PNO_sigma_component_diagnostics.txt
```

Read it together with:

```text
PNO_sigma_epsilon_gamma.txt
PNO_sigma_fig4a_point_compare.txt
Fig2_sigma_fits.png
Fig3_sigma_magnetic_map.png
Fig4_sigma_epsilon_gamma_scatter.png
```
