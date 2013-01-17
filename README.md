
# Hapke model

Based on Hapke (2002).

## Usage

Example:

```
using HapkeModel
model = ScatteringModel(HenyeyGreenstein(0.5), 0.1)
r = BRDF(model, 0.0, 0.0, 0.0)
```

To define a scattering model, call `model = ScatteringModel(P, w)`, where P is a `PhaseFunction` and `w` is the single-scattering albedo.

There are four options for `PhaseFunction`:
- `Isotropic()`, for isotropic scattering
- `Rayleigh()` for Rayleigh scattering
- `HenyeyGreenstein(xi)`, for Henyey-Greenstein scattering, where `xi` is the asymmetry parameter.
- `DoubleHenyeyGreenstein(c, xi)`, for Henyey-Greenstein scattering, where `c` is the weight parameter and `xi` is the asymmetry parameter.

To compute the value of the Hapke BDRF for a given scattering model, call `BDRF(model, mu0, mu, g)` where `model` is a `ScatteringModel` object, `mu0` and `mu` are the cosines of incidence and emergence, and `g` is the phase angle.

