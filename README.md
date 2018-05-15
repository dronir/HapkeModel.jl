
# Hapke model

Based on Hapke (2002). Optionally including the shadow-hiding opposition effect term.

## Usage

Example:

```
using HapkeModel
model = ScatteringModel(HenyeyGreenstein(0.5), 0.1)
r = BRDF(model, 0.0, 0.0, 0.0)
```

To define a scattering model without the Shadow-Hiding Opposition Effect (SHOE) term, call `model = ScatteringModel(P, w)`, where P is a `PhaseFunction` and `w` is the single-scattering albedo.

To use the SHOE term, define the model with:
* `model = ScatteringModel(P, w, hs)`, where P is a `PhaseFunction`, `w` is the single-scattering albedo, and `hs` is the shadowing parameter (eq. 28-30 in Hapke, 2000), or
* `model = ScatteringModel(P, w, E, a, phi)`, where P is a `PhaseFunction`, `w` is the single-scattering albedo, `E` is the extinction coefficient in the medium, `a` the mean particle radius, and `phi` the filling factor (eq. 30 in Hapke, 2000).

Be wary of the actual physical interpretation of the parameter values.

There are four options for `PhaseFunction`:
- `Isotropic()`, for isotropic scattering
- `Rayleigh()` for Rayleigh scattering
- `HenyeyGreenstein(xi)`, for Henyey-Greenstein scattering, where `xi` is the asymmetry parameter.
- `DoubleHenyeyGreenstein(c, xi)`, for Henyey-Greenstein scattering, where `c` is the weight parameter and `xi` is the asymmetry parameter.

To compute the value of the Hapke BDRF for a given scattering model, call `BDRF(model, mu0, mu, g)` where `model` is a `ScatteringModel` object, `mu0` and `mu` are the cosines of incidence and emergence, and `g` is the phase angle.

