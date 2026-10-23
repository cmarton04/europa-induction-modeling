# europa-induction-modeling

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

Modeling Europa's time-dependent magnetic induction response and subsurface ocean characterization — NASA ICONS Internship 2026

---

# About This Project

Jupiter's magnetosphere acts as a natural laboratory for probing the interior of Europa.
The tilt of Jupiter's magnetic field relative to its rotation axis causes the magnetic equator
to sweep up and down across Europa's orbital path, creating a time-varying magnetic wave.
This wave drives electric currents in Europa's subsurface ocean, producing a measurable
magnetic induction response: a signal that carries information about the ocean's
depth, thickness, and conductivity.

The goal of this project is to build a physics-based model of that induction signal and
use it to better characterize Europa's subsurface ocean.

## Method

The model implements the analytical induction framework of Zimmer et al. (2000): an
elliptically polarized primary field from Jupiter, written in complex phasor notation,
drives a secondary (induced) field in Europa's ocean. The induced field behaves like the
response of a perfect conductor, scaled down by an amplitude factor and delayed by a phase lag.


The model is flown along the real trajectories of the Galileo E4, E14, E19, and E26 flybys
and compared with the measured magnetometer data (residuals = model − measured). A
parameter-sensitivity test checks how the predicted field changes under three extreme
interior scenarios: weak induction, thick ice with a deep ocean, and a minimal ocean shell.

---

# Key Results

- **The model reproduces Europa's induction signature in real Galileo data.** The
  predicted field components follow the measured trends along the flybys, including the
  structure near closest approach. To line the model up with the data, its clock is
  phase-aligned to Jupiter's rotational phase at the time of each flyby.

- **The data favor strong induction over weak induction.** In the parameter-sensitivity
  test, the weak-induction scenario misses the closest-approach structure in the field,
  while the stronger-response scenarios capture it.

- **The results are consistent with a conducting subsurface ocean.** This agrees with
  Zimmer et al. (2000), who found that matching the Galileo data requires a salty,
  conductive layer close to the surface. The amplitude and phase lag together act as a
  probe of ocean depth and salinity.

- **The model diverges most in B_z.** B_x and B_y are dominated by the induction signal
  and match more cleanly. B_z is more sensitive to physics this model does not include:
  interaction with the ambient plasma (Alfvén-wing currents and local plasma
  perturbations), uncertainty in the fitted Jovian background field, and small-scale,
  short-timescale fluctuations too abrupt to come from a smooth internal field.

## Limitations and Next Steps

This is a first-principles analytical induction model, not a full plasma simulation.
Planned extensions include the ~85-hour periodic component in B_z and coupling the
induction model to an MHD simulation to capture plasma contributions.

---

# Getting Started

**Requirements:** MATLAB R2026a (MATLAB Online also works).

**To run:**
1. Clone this repository.
2. In MATLAB, set the repository folder as the current folder.
3. Open `matlab/main proj/europa_flyby_model.mlx` (the main model).
4. Select flyby id or use custom file
5. Click **Run**.

The visualization scripts in `secondary field codes/` and `total field codes/` are run the
same way. Pre-rendered outputs are in `figures/main proj/`. Plain `.m` versions of the
main project code are in `matlab/main proj/.m files/`.

---

# Repository Structure

```
europa-induction-modeling/
├── figures/                                      ← Output plots and visualizations
│   ├── main proj/
│   │   ├── flyby figures/                        ← Model outputs for Galileo flybys E4, E14, E19, E26
│   │   │   ├── E4/ E14/ E19/ E26/                ← One folder per flyby, each containing:
│   │   │   │     <flyby>_field_graphs.png        Stackplot of measured field with model prediction overlaid
│   │   │   │     <flyby>_flyby_geometry.mp4      Flyby trajectory through Europa's magnetic environment
│   │   │   │     <flyby>_residuals.png           Residuals (model − measured)
│   │   │   │     <flyby>_parameter_test.png      Model under 3 extreme parameter sets, measured in grey
│   │   │   │     <flyby>_parameter_error.png     Error under the 3 extreme parameter sets
│   │   │   └── parameter_testing_legend.png      The three extreme scenarios (legend)
│   │   ├── old figures/                          ← Figures from earlier stages before final output
│   │   ├── secondary field figures/              ← Europa's secondary (induced) field
│   │   │   ├── bsec_heatmaps.mp4                 Equatorial-plane heatmaps over time
│   │   │   ├── bsec_heatmaps_fieldlines.mp4      Same, with field lines
│   │   │   ├── bsec_t0_2d.png                    Snapshot of secondary field lines
│   │   │   ├── bsec_t_complx.mp4                 Over time (complex B_sec equation)
│   │   │   └── bsec_t_real.mp4                   Over time (real B_sec equation, see notes)
│   │   └── total field figures/                  ← Total field = Jupiter's primary + secondary
│   │       ├── btotal_t_colormap.mp4             Streamlines colored by field strength
│   │       ├── btotal_withvectors.mp4            Plus induced moment & B_prim vectors
│   │       ├── by_vs_bx_param.png                Parametric curve of B_y vs B_x
│   │       └── delta_b_colormap.mp4              Streamlines colored by induced deviation
│   └── practice codes/                           ← Figures from preparation exercises
├── matlab/                                       ← MATLAB scripts and Live Scripts (.mlx)
│   ├── main proj/
│   │   ├── europa_flyby_model.mlx                ← MAIN MODEL: flies spacecraft along Galileo trajectories
│   │   ├── .m files/                             ← .m versions of the main project code
│   │   ├── flyby data/                           ← Galileo trajectory/magnetometer TAB files (model input)
│   │   ├── secondary field codes/                ← Heatmaps, field lines, B_sec validity checks
│   │   ├── total field codes/                    ← Colormaps, vectors, B_y vs B_x ellipse
│   │   └── europaimage.jpg                       ← Europa surface texture for visualizations
│   └── practice codes/                           ← Preparation exercise code
├── notes/                                        ← Reference materials
│   ├── Analytical equations for bsec & bprim.pdf   Derivation of the real primary & secondary equations
│   └── Gyromotion exercise (analytical).pdf        Practice exercise derivation
├── .gitignore
├── CITATION.cff
├── LICENSE
└── README.md
```

---

# Data

Galileo trajectory and magnetometer data are from the NASA Planetary Data System (PDS).
Coordinates are in the EPhiO frame.

# Reference

Zimmer, C., Khurana, K. K., & Kivelson, M. G. (2000). Subsurface oceans on Europa and
Callisto: Constraints from Galileo magnetometer observations. *Icarus*, 147, 329–347.

# Citation

If you use this code, please cite it using `CITATION.cff`, or:

> Marton, C. (2026). *Europa Induction Modeling* (v1.0.0) [Software].
> Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX

# License

MIT. See `LICENSE`.

---

# Author

**Chayla Marton**

Astrophysics student | Florida Institute of Technology
NASA ICONS Internship, Summer 2026
GitHub: [@cmarton04](https://github.com/cmarton04)

# Acknowledgments

Developed during the NASA Europa ICONS internship under Dr. Carol Paty
(University of Oregon).
