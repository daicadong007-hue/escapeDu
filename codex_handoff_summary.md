# Codex Handoff Summary: Chen 2011 Gear Dynamic Model

## Goal

Reproduce the analytical spur gear dynamic model from `chen2011.pdf` for a cracked gear tooth with constant crack depth along the full tooth width, using and improving the MATLAB code in `code/`.

## Workspace

- Root folder: `C:\Users\Hi Windows 11 23\Desktop\escapedu`
- Paper: `chen2011.pdf`
- Original MATLAB program: `code\geardynamic.m`
- Original stiffness helpers: `code\toothmesh1.m`, `code\toothmesh2.m`
- Bisection helper: `code\bsm2.m`
- New reproduction program: `code\reproduce_chen2011_constant_depth.m`
- Generated outputs: `outputs\chen2011_constant_depth\`

## Important Paper Interpretation

The requested crack model is the constant through-width crack model from Chen/Shao Section 1.1, not the parabolic crack-depth distribution along tooth width from the later sections. The reproduction therefore targets constant crack depth `q0` across the full `20 mm` tooth width and generates analogous dynamic/stiffness results.

## Main Implementation Done

1. Added `code\reproduce_chen2011_constant_depth.m`.
   - Uses Chen 2011 paper-aligned parameters:
     - Pinion/gear teeth: `30/25`
     - Module: `2 mm`
     - Tooth width: `20 mm`
     - Pressure angle: `20 deg`
     - Young modulus: `2e11 Pa`
     - Poisson ratio: `0.3`
     - Pinion/gear speed: `2000/2400 rpm`
     - Load torque: `60 Nm`
     - Bearing stiffness: `6.56e8 N/m`
     - Bearing damping: `1.8e3 Ns/m`
     - Mesh damping: `67 Ns/m`
     - Friction coefficient: `0.06`
   - Computes healthy and cracked single-pair mesh stiffness using `toothmesh1.m` / `toothmesh2.m`.
   - Assembles time-varying mesh stiffness with contact-ratio overlap.
   - Runs a 6-DOF ODE gear dynamic model.
   - Saves figures, CSVs, and `.mat` data.

2. Fixed `code\bsm2.m`.
   - Changed the function declaration from `function Pasokh=bsm(...)` to `function Pasokh=bsm2(...)` so the file/function names match.

3. Updated original `code\geardynamic.m`.
   - User requested changing its ODE from always adding primary + adjacent tooth pair to adding the adjacent pair only during contact-ratio overlap.
   - Added:
     - `params.overlap_start_index = Pb;`
   - In `gear_dyn`, `N_adjacent` now starts at `0`.
   - Adjacent pair is only computed when:
     - `idx > overlap_start_index`
   - This preserves the old program structure while fixing the double-contact logic.

## New Output Files

Under `outputs\chen2011_constant_depth\`:

- `chen2011_constant_depth_results.mat`
- `summary_constant_depth.csv`
- `table3_validation.csv`
- `program_comparison.md`
- `single_pair_stiffness.png`
- `total_mesh_stiffness.png`
- `time_mesh_stiffness_rotational_period.png`
- `angle_mesh_stiffness_rotational_period.png`
- `time_mesh_stiffness_rotational_period.csv`
- `pinion_y_displacement.png`
- `frequency_spectrum.png`
- `rms_kurtosis_change.png`
- `residual_indicators.png`
- `sideband_marker_change.png`

## Validation Results

The new program validates single-pair stiffness against Chen 2011 Table 3 using minimum stiffness over the mesh path:

| Case | Computed Min | Paper Value | Error |
|---|---:|---:|---:|
| Healthy | `1.6255e8 N/m` | `1.52e8 N/m` | `6.94%` |
| Crack No. 1 | `1.5923e8 N/m` | `1.47e8 N/m` | `8.32%` |
| Crack No. 2 | `1.4293e8 N/m` | `1.38e8 N/m` | `3.57%` |

## Important Difference Between Original and New Program

Original `geardynamic.m`:

- Uses pinion/gear teeth `25/30`.
- Uses pinion speed `2400 rpm`.
- Uses bearing stiffness `6.56e9 N/m`.
- Uses damping ratio `zeta = 0.07`.
- Previously always added a healthy adjacent pair in the ODE.

New reproduction program:

- Uses paper-aligned pinion/gear teeth `30/25`.
- Uses pinion/gear speed `2000/2400 rpm`.
- Uses bearing stiffness `6.56e8 N/m`.
- Uses constant mesh damping `C_m = 67 Ns/m`.
- Adds adjacent pair only during contact-ratio overlap.
- Allows cracked tooth to affect the active or entering adjacent pair during a rotational period.
- Saves outputs instead of only showing interactive MATLAB figures.

## Recent User Request Completed

User asked:

> change the original program from "Always adds a primary pair and a healthy adjacent pair" to "Adds the adjacent pair only during the contact-ratio overlap" in the ODE

Completed in `code\geardynamic.m`.

Relevant edited behavior in `gear_dyn`:

```matlab
N_adjacent = 0;
if idx > overlap_start_index
    adjacent_idx = idx - overlap_start_index;
    k_adjacent = interp1(1:PTH, K_healthy, adjacent_idx, 'linear', 'extrap');
    c_adjacent = 2 * zeta * sqrt(k_adjacent * J_e);
    N_adjacent = max(k_adjacent * delta + c_adjacent * dot_delta, 0);
end
```

## MATLAB Verification Notes

- MATLAB is installed at `D:\Matlab\bin\matlab.exe`.
- Running MATLAB commands works.
- MATLAB prints permission warnings related to:
  - `C:\Users\Hi Windows 11 23\AppData\Roaming\MathWorks\MATLAB\R2023a\...`
- These warnings did not prevent script execution or output generation.
- `reproduce_chen2011_constant_depth.m` ran successfully and saved outputs.
- `checkcode('geardynamic.m','-id')` reported style/settings warnings but no syntax failure.

## Suggested Next Steps

- If continuing on the second device, first open:
  - `code\geardynamic.m`
  - `code\reproduce_chen2011_constant_depth.m`
  - `outputs\chen2011_constant_depth\program_comparison.md`
- To rerun the reproduction:

```matlab
cd('C:\Users\Hi Windows 11 23\Desktop\escapedu\code')
reproduce_chen2011_constant_depth
```

- To rerun the original updated program:

```matlab
cd('C:\Users\Hi Windows 11 23\Desktop\escapedu\code')
geardynamic
```
