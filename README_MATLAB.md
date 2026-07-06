# OpenNTN — MATLAB Adaptation

A MATLAB-native implementation of the OpenNTN stochastic channel model for
Non-Terrestrial Networks (NTN), combining the 3GPP TR 38.811 / TR 38.901
geometry-based stochastic channel model (GSCM) with the ITU-R atmospheric
propagation models (P.676 gases, P.618 scintillation). It is a faithful port of
the core stochastic workflow of the Python [OpenNTN](../README.md) project and
reuses the **same** TR 38.811 parameter tables in `OpenNTN/models/*.json`.

> The original Python project delegates the full polarized-MIMO channel
> coefficient generation to NVIDIA Sionna / TensorFlow. This MATLAB layer is a
> self-contained, dependency-free implementation of the stochastic core: it is
> intended for NTN link-level studies (link budgets, power-delay profiles,
> Doppler, SISO frequency responses) and as a base for progressively porting the
> remaining Sionna-dependent behaviour. See [Scope and simplifications](#scope-and-simplifications).

---

## 1. Quick start

```matlab
% From the repository root:
addpath("matlab");

% 1. Create a scenario: dense_urban | urban | sub_urban, S-band or Ka-band,
%    uplink/downlink, elevation angle in [10, 90] deg.
scenario = openntn.createScenario("urban", 2.0e9, "downlink", 50);

% 2. Place a satellite (BS) and ground users (UTs). Heights in metres.
utLoc = [  0    0  1.5;
          80   30  1.5;
         200 -120  1.5];
bsLoc = [0 0 600e3];                 % LEO satellite at 600 km
scenario = openntn.setTopology(scenario, utLoc, bsLoc, ...
    "UTVelocities", [1 0 0; 0 0 0; 5 0 0], ...
    "Indoor", [false false false], ...
    "LOS", []);                      % [] -> sample LOS from TR 38.811 probability

% 3. Generate a channel realization (32 time samples at 1 kHz).
channel = openntn.generateChannel(scenario, ...
    "NumTimeSamples", 32, "SamplingFrequency", 1e3);

disp(channel.pathloss.total)         % [numBS x numUT] total path loss [dB]
disp(channel.satelliteDopplerHz)     % peak satellite Doppler [Hz]
```

Two runnable scripts are provided:

| Script | Purpose |
| --- | --- |
| [`examples/run_basic_ntn.m`](examples/run_basic_ntn.m) | Minimal end-to-end example with a PDP plot. |
| [`tests/test_openntn.m`](tests/test_openntn.m) | Self-checking sanity/validation suite (no toolboxes required). |

---

## 2. What the model produces

Given a scenario and a topology, the pipeline implements 3GPP TR 38.901 §7.5
(parametrized for NTN by TR 38.811 §6) together with the TR 38.811 §6.6 path
loss:

```
createScenario ─▶ setTopology ─▶ sampleLSP ─▶ sampleRays ─▶ generateChannel
   (tables)        (geometry,     (step 4:      (steps 5-9:    (per-cluster
                    LOS, PL)       7 LSPs)       delays/powers/  SISO CIR +
                                                 angles/XPR)     Doppler)
```

* **Large-scale parameters (LSPs)** — delay spread (DS), azimuth/zenith spreads
  of arrival and departure (ASA, ASD, ZSA, ZSD), shadow fading (SF) and Rician
  K-factor — sampled as a correlated log-normal vector (TR 38.901 step 4).
* **Rays** — cluster delays (step 5), cluster powers (step 6), per-ray
  azimuth/zenith angles of arrival/departure (step 7), random coupling (step 8),
  and cross-polarization power ratios (step 9).
* **Path loss** — free-space loss on the spherical-Earth slant range, shadow
  fading and clutter loss, plus ITU-R gaseous absorption and scintillation.
* **Doppler** — closed-form peak satellite Doppler (TR 38.811 §5.3.4) plus the
  user-motion Doppler projected onto each cluster.

---

## 3. API reference

All functions live in the `+openntn` package (call them as `openntn.<name>`).

### Public workflow

| Function | Summary |
| --- | --- |
| `createScenario(type, fc, direction, elevation, Name,Value)` | Build a scenario struct; validates band/elevation and loads the LOS/NLOS tables. |
| `setTopology(scenario, utLoc, bsLoc, Name,Value)` | Attach geometry, sample LOS, compute path loss. |
| `sampleLSP(scenario)` | Draw one correlated set of large-scale parameters. |
| `sampleRays(scenario, lsp)` | Draw cluster delays, powers, angles and XPR. |
| `generateChannel(scenario, Name,Value)` | Full pipeline → compact SISO CIR with Doppler evolution. |
| `pathloss(scenario)` | Total path loss and its components [dB]. |
| `satelliteDoppler(scenario)` | Peak satellite Doppler shift [Hz]. |

### Helpers (continue/utility)

| Function | Summary |
| --- | --- |
| `generateSingleSectorTopology(scenario, numUT, Name,Value)` | Randomly drop UTs in a sector (MATLAB analogue of the Python `gen_single_sector_topology`). |
| `cirToOFDMChannel(channel, frequencies)` | Convert the cluster CIR to an OFDM frequency response. |
| `defaultEnvironment()` | ITU-R default atmospheric / earth-station parameters. |

### Internal

`selectModelFiles`, `loadModelParameters`, `param`, `gasAttenuation676`,
`scintillationAttenuation` are used internally but documented in their headers.

### `createScenario` options

| Name | Default | Meaning |
| --- | --- | --- |
| `EnablePathloss` | `true` | Apply deterministic path loss to the channel amplitude. |
| `EnableShadowFading` | `true` | Apply the log-normal shadow-fading amplitude. |
| `DopplerEnabled` | `true` | Include satellite + user Doppler in the CIR. |
| `AverageStreetWidth` | `20.0` | Street width [m] (urban geometry parameter). |
| `AverageBuildingHeight` | `5.0` | Building height [m] (urban geometry parameter). |

Supported carrier bands: **S-band** 1.9–4 GHz and **Ka-band** 19–40 GHz, matching
the ranges asserted by the Python `SystemLevelScenario`.

---

## 4. Scope and simplifications

This layer reproduces the **stochastic core** of OpenNTN. The following reference
behaviours are intentionally simplified or omitted; each is a candidate for
future work:

| Area | Reference (Python/Sionna) | This MATLAB layer |
| --- | --- | --- |
| Channel coefficients | Full polarized MIMO antenna-field synthesis (`channel_coefficients.py`) | Per-cluster **SISO** taps with Doppler; XPR/angles are produced but not yet combined with antenna patterns. |
| Spatial LSP correlation | Second Cholesky correlating LSPs across UTs by distance | Cross-LSP correlation only; UTs treated independently. |
| CDL / TDL models | Dummy files present, not parametrized below 50° | Not ported (same status as the reference). |
| Indoor / O2I, HAPS | Stubs kept for future standard updates | LOS forced for indoor UTs; entry loss = 0. |
| Cloud / rain loss | `compute_pathloss_additional` (disabled in reference) | Not included (matches the reference default). |

These choices keep the implementation dependency-free (no toolboxes required) and
focused on NTN link-level research. None of them affect the supported S/Ka-band
satellite scenarios beyond the documented limits.

---

## 5. Correctness notes and parity with the Python reference

The port was validated function-by-function against `lsp.py`, `rays.py`,
`utils.py`, and the scenario classes. Key points a reviewer should know:

* **7-LSP ordering and the 7×7 cross-correlation map** (`sampleLSP.m`) reproduce
  `LSPGenerator._compute_cross_lsp_correlation_matrix` index-for-index.
* **ZoD departure centring** (`sampleRays.m > zenithAngles`) distinguishes the
  LOS first-cluster re-centring from the NLOS `losAngle + zod_offset` form, matching
  `RaysGenerator._zenith_angles`. (An earlier version of the MATLAB port handled
  only the NLOS form and dropped `zod_offset`; this is now fixed. `zod_offset`
  itself is computed in `setTopology.m`.)
* **Ray offset table** intentionally keeps the reference's `-0.1481` entry at
  index 16 (the 3GPP table value is `-1.1481`); this matches the Python/Sionna
  source so realizations are comparable. Change it in `sampleRays.m` if strict
  TR 38.901 Table 7.5-3 values are required.
* **Gas attenuation** (`gasAttenuation676.m`) uses the ITU-R P.676 exponent
  switch `n = 0` (f ≤ 60 GHz) / `n = -15` (f > 60 GHz) in the 60 GHz band, as in
  `utils.compute_pathloss_gas`. This band is outside the supported NTN bands and
  is provided only for completeness.
* **Slant range** uses the spherical-Earth form of TR 38.811 Eq. 6.6-3, driven by
  the satellite height in `bsLoc(:,3)`.

---

## 6. References

* **3GPP TR 38.811** — *Study on New Radio (NR) to support non-terrestrial
  networks* (channel model §6, geometry §6.6, Doppler §5.3.4).
* **3GPP TR 38.901** — *Study on channel model for frequencies from 0.5 to
  100 GHz* (GSCM procedure §7.5, parameter tables §7.5-x).
* **ITU-R P.676** — *Attenuation by atmospheric gases.*
* **ITU-R P.618** — *Propagation data and prediction methods … Earth-space
  telecommunication systems* (scintillation).
* T. Düe, M. Vakilifard, C. Bockelmann, D. Wübben, A. Dekorsy, *OpenNTN: An
  Open-Source Framework for Non-Terrestrial Network Channel Simulations*, WSA 2025.

---

## 7. Repository layout

```
matlab/
├── +openntn/                 % the package (call functions as openntn.<name>)
│   ├── createScenario.m      % scenario construction + validation
│   ├── setTopology.m         % geometry, LOS, ZoD offset, path loss
│   ├── sampleLSP.m           % step 4: correlated large-scale parameters
│   ├── sampleRays.m          % steps 5-9: delays, powers, angles, XPR
│   ├── generateChannel.m     % full pipeline -> SISO CIR + Doppler
│   ├── pathloss.m            % TR 38.811 6.6 link budget
│   ├── gasAttenuation676.m   % ITU-R P.676 gaseous absorption
│   ├── scintillationAttenuation.m  % ITU-R P.618 scintillation
│   ├── satelliteDoppler.m    % TR 38.811 5.3.4 Doppler
│   ├── generateSingleSectorTopology.m % random UT drop helper
│   ├── cirToOFDMChannel.m    % CIR -> OFDM frequency response
│   ├── param.m, selectModelFiles.m, loadModelParameters.m, defaultEnvironment.m
├── examples/run_basic_ntn.m
├── tests/test_openntn.m
└── README_MATLAB.md          % this file
```
