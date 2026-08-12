# Physics-informed SINDy-PI for lettuce-greenhouse dynamics

This repository contains MATLAB code and preprocessed input data for identifying a four-state lettuce-greenhouse model with a physics-informed parallel implicit sparse identification of nonlinear dynamics (SINDy-PI) workflow.

The code uses greenhouse biophysical knowledge to build a state-specific candidate library, applies sequentially thresholded least squares (STLSQ), selects models using derivative-prediction error on a separate 40-day interval, and evaluates the selected equations on a third 40-day interval. A second script studies sensitivity to additive Gaussian noise.

## Model signals and units

The ordering below follows the MATLAB implementation and must be kept consistent in the manuscript, data files, candidate libraries, and generated equations.

### States

| Index | Symbol | Meaning | Unit used by the code |
|---:|---|---|---|
| 1 | `x1` | Lettuce dry matter | kg m^-2 |
| 2 | `x2` | Indoor CO2 density | kg m^-3 |
| 3 | `x3` | Indoor air temperature | degC |
| 4 | `x4` | Indoor water-vapour density | kg m^-3 |

Some plot labels display `x1` in g m^-2 even though the nominal model state and initial condition are in kg m^-2. Convert by a factor of 1000 if reporting dry matter in g m^-2.

### Control inputs

| Index | Symbol | Meaning | Unit |
|---:|---|---|---|
| 1 | `u1` | CO2 injection rate | mg m^-2 s^-1 |
| 2 | `u2` | Ventilation rate | mm s^-1 |
| 3 | `u3` | Heating input | W m^-2 |

### Weather disturbances

| Index | Symbol | Meaning before/after preprocessing | Unit after preprocessing |
|---:|---|---|---|
| 1 | `d1` | Outdoor global radiation | W m^-2 |
| 2 | `d2` | Outdoor CO2 concentration -> density | kg m^-3 |
| 3 | `d3` | Outdoor temperature | degC |
| 4 | `d4` | Outdoor relative humidity -> vapour density | kg m^-3 |

The code uses `d2` for outdoor CO2 and `d3` for outdoor temperature. Reversing these two signals changes both the nominal dynamics and the physical interpretation of the identified terms.

## Data

The pre-split files contain three consecutive 40-day intervals:

| Role | Interval | Control file and variable | Disturbance file and variable |
|---|---|---|---|
| Training | days 0-40 | `u_0_40Days.mat` / `u_0_40Days` | `d_0_40Days.mat` / `d_0_40Days` |
| Model selection (called `test` in code) | days 40-80 | `u_40_80Days.mat` / `u_40_80Days` | `d_40_80Days.mat` / `d_40_80Days` |
| Final validation | days 80-120 | `u_80_120Days.mat` / `u_80_120Days` | `d_80_120Days.mat` / `d_80_120Days` |

Each control matrix is `3 x 3842`; each disturbance matrix is `4 x 3865`. The scripts transpose and crop them to 3841 samples, corresponding to

```matlab
h = 15*60;                 % 900 s
t = 0:h:40*24*60*60;      % 3841 samples, including both endpoints
```

The control trajectories were produced previously by an MPC simulation and are loaded as fixed inputs. This repository does not contain the MPC implementation that generated them. The disturbance trajectories are based on outdoor weather measurements from a Venlo-type greenhouse in Bleiswijk, the Netherlands. See `bin/data/README.docx` for source-data details.

The raw file `outdoorWeatherWurGlas2014.mat` is provenance material and is not directly loaded by either main experiment script.

## Requirements

The manuscript reports MATLAB R2021b. The static checks for this README were also performed with MATLAB R2024b.

Required MathWorks products:

- MATLAB;
- Symbolic Math Toolbox (`sym`, `vpa`, `vpasolve`, and `matlabFunction`);
- Parallel Computing Toolbox is recommended because the STLSQ normalization loop uses `parfor`. MATLAB may execute this loop serially when no pool is open, depending on the installation.

Project-specific model and conversion functions included under `Functions/`:

- `DefineParameters.m`: returns the nominal greenhouse-model parameter structure, including all fields required by `Functions/GHLettuce_ODE.m`;
- `rh2vaporDens.m`: converts temperature and relative humidity to water-vapour density;
- `vaporDens2rh.m`: inverse conversion, used only by the local output-conversion helper and commented plotting code.

These three functions were checked with MATLAB R2024b. `DefineParameters` contains all 22 fields referenced by the nominal ODE, the humidity conversion round trip had a maximum absolute error of approximately `7.1e-15`, and a 15-minute `ode45` nominal-model smoke test completed with finite states and derivatives. MATLAB Code Analyzer reported no warnings for these three files.

Before running an experiment, verify the installation from the repository root:

```matlab
codeDir = fullfile(pwd, "SINDyPI PIML Lettuce Greenhouse - Fakhira's Thesis");
cd(codeDir)
addpath("Functions")

assert(~isempty(which("DefineParameters")), "Missing DefineParameters.m")
assert(~isempty(which("rh2vaporDens")), "Missing rh2vaporDens.m")
assert(~isempty(which("vaporDens2rh")), "Missing vaporDens2rh.m")
```

## Computational workflow

```mermaid
flowchart LR
    A[Load three input/weather intervals] --> B[Convert outdoor CO2 and RH]
    B --> C[Simulate nominal greenhouse with ode45]
    C --> D[Build state-specific LHS and RHS libraries]
    D --> E[STLSQ over candidate LHS terms and lambda values]
    E --> F[Select by derivative error on days 40-80]
    F --> G[Generate Sindy_ODE_RHS_GHLettuce.m]
    G --> H[Evaluate on days 80-120]
```

### 1. Data generation

`Get_Sim_Data_GHLettuce` advances the four-state nominal model over each 15-minute interval with `ode45`. The input and disturbance are held constant over each interval. Although `ode45` uses adaptive internal steps, outputs are retained at the 15-minute grid.

The function records the state derivative by evaluating the nominal ODE right-hand side at each sampled state. It does **not** compute derivatives by forward finite differences or TV regularization.

### 2. Physics-informed libraries

The candidate variables are restricted state by state:

| State equation | Variables retained by the code |
|---|---|
| `dx1/dt` | `x1`, `x2`, `x3`, `d1` and selected interactions |
| `dx2/dt` | `x1`, `x2`, `x3`, `d1`, `d2`, `u1`, `u2` and selected interactions |
| `dx3/dt` | `x3`, `d1`, `d3`, `u2`, `u3` and selected interactions |
| `dx4/dt` | `x1`, `x3`, `x4`, `d4`, `u2` and selected interactions |

The `*_Test.m` files enumerate broader LHS/RHS combinations. The `*_Best.m` files encode hand-selected candidate terms found during previous searches. These names describe library breadth and do not correspond directly to the training/test/validation data split.

### 3. Sparse regression and model selection

The main noise-free script uses:

```matlab
lam = [1e-20; 1e-19; 1e-18; 1e-17; 1e-16; ...
       1e-15; 1e-14; 1e-13; 1e-12];
N_iter = 20;
NormalizeLib = 1;
```

For each state and LHS guess, `sparsifyDynamics_GHLettuce`:

1. normalizes each RHS-library column by its Euclidean norm;
2. obtains an initial least-squares coefficient vector;
3. repeats thresholding (`abs(Xi) < lambda`) and least-squares refitting 20 times;
4. rescales the coefficients to the original library units; and
5. reconstructs an explicit symbolic state equation with `vpasolve`.

Candidate models are selected on days 40-80 using

```matlab
relative_error = norm(dData - dData_est) / norm(dData);
```

This is a one-step/right-hand-side derivative error, not a free-run state-trajectory error. The selected equations are then evaluated on days 80-120.

## Running the experiments

Always start MATLAB in the code directory because the scripts use relative paths and generate a MATLAB function in the current directory.

### Noise-free library/LHS search

```matlab
cd("SINDyPI PIML Lettuce Greenhouse - Fakhira's Thesis")
SINDyPI_Greenhouselettuce_PhysML_BruteforceSearch
```

Important workspace outputs include:

| Variable | Description |
|---|---|
| `ODE_Best` | Selected symbolic equation for each of the four states |
| `BEST_MODEL` | LHS guess, library description, lambda, validation score, and equation |
| `Score` | Model-selection error indexed by state, LHS guess, and lambda |
| `dScore_all_Best` | Four validation derivative errors |
| `Data_True_val`, `Data_Es_val` | Nominal and identified validation state trajectories |
| `dData_True_val`, `dData_Es_val` | Nominal and identified validation derivatives |

The script creates figures but does not save them automatically. It also creates a `Results` directory but does not write results into it.

### Noise-sensitivity experiment

The distributed script currently evaluates only `Noise = 1e-7`. To reproduce the intended sweep, replace that assignment with:

```matlab
Noise = [0; 1e-12; 1e-11; 1e-10; 1e-9; 1e-8; 1e-7];
rng(1, "twister");  % add a documented seed before any simulations
SINDyPI_Greenhouselettuce_PhysML_Noise
```

The noise script fixes `lambda = 1e-16` for all four state equations. Its final summary variable is `dScore_X_Noise`, with rows corresponding to `x1` through `x4` and columns corresponding to the entries of `Noise`.

In the current implementation, `Noise` is the absolute standard deviation of zero-mean Gaussian noise added separately to every state and derivative channel in their native units:

```matlab
Data(:,i)  = x_list(:,i)  + Noise*randn(...);
dData(:,i) = dx_list(:,i) + Noise*randn(...);
```

It is not a percentage, signal-to-noise ratio, or variance normalized by the scale of each variable. This matters because the four states and their derivatives have very different units and magnitudes.

## Generated files and path precedence

`Generate_ODE_RHS_GHLettuce` calls `matlabFunction` and overwrites

```text
./Sindy_ODE_RHS_GHLettuce.m
```

for every candidate evaluated during model selection and again for the final model. Because the repository root is MATLAB's current directory, this root copy takes precedence over the duplicate copy in `Functions/`.

Do not manually edit the generated file as a source of truth. Preserve final equations separately (for example, a versioned MAT/JSON/text result file) before rerunning the search.

## Mapping code to manuscript evidence

| Manuscript item | Available code path | Status |
|---|---|---|
| Nominal simulated datasets | `GHLettuce_ODE.m`, `DefineParameters.m`, conversion functions, `Get_Sim_Data_GHLettuce.m` | Dependencies present; nominal ODE and 15-minute integration smoke-tested |
| Physics-informed library search | `...BruteforceSearch.m`, `GuessLib_*`, `SINDyLib_*` | Present |
| Lambda grid and STLSQ iterations | Main scripts and `sparsifyDynamics_GHLettuce.m` | Present |
| Validation predictions and derivative errors | Main scripts and generated RHS | Present |
| Noise sensitivity | `...Noise.m` | Present, but full noise vector is commented out and no RNG seed is set |
| Conventional SINDy-PI versus physics-informed SINDy-PI comparison | Not found in the supplied source tree | Not reproducible from the current repository |
| MPC generation of control signals | Not found; only precomputed `u_*.mat` files are supplied | Inputs usable, generation not reproducible |

## Known issues and pre-release checklist

The following items should be resolved before citing this repository as a fully reproducible research artifact.

1. **Make primary function names match file names.** MATLAB's Code Analyzer reports mismatches in:
   - `Generate_ODE_RHS_GHLettuce.m` (`Generate_ODE_RHS_DW` in the declaration);
   - `sparsifyDynamics_GHLettuce.m` (`sparsifyDynamicsDW` in the declaration);
   - `GuessLib_GHLettuce_PhysML_Best.m` (`GuessLib_GHTomato_PhysML_Best` in the declaration).
2. **Resolve the duplicate generated RHS.** Keep one generated `Sindy_ODE_RHS_GHLettuce.m`, or write generated models to a dedicated output directory with unique names.
3. **Add deterministic noise control.** Call `rng` with a documented seed and save it with every result. At present, repeated noise experiments need not agree.
4. **Define noise mathematically.** State that the current `eta` is an absolute standard deviation in native units, or preferably normalize every channel and define a dimensionless relative noise level/SNR.
5. **Do not add independent noise to both sides of a validation comparison.** The current noise script generates new random perturbations for the nominal and identified validation simulations. Evaluate both models at the same clean validation states/inputs, or use the same explicitly stored noisy observations.
6. **Fix the shuffle branch.** `Get_Sim_Data_GHLettuce` refers to undefined `Data_Es` when `Shuffle == 1`; leave `Shuffle = 0` until this is corrected.
7. **Document the derivative source consistently.** The code evaluates exact model derivatives; it does not use the forward-difference expression described in the manuscript.
8. **Add the conventional SINDy-PI baseline script.** The supplied source does not reproduce the conventional-baseline results or timing comparison reported in the revised manuscript.
9. **Save complete artifacts automatically.** At minimum save configuration, library terms, coefficients, selected LHS, lambda, random seed, errors, MATLAB release, runtime, and figures under a timestamped result directory.
10. **Add automated tests.** The three newly included dependency files pass manual smoke tests. Add permanent tests for unit-conversion round trips, dimensions/units of all MAT inputs, one-step nominal-model evaluation, library column/symbol alignment, and recovery of the polynomial indoor-temperature equation.
11. **Add a license and citation metadata.** No `LICENSE` or `CITATION.cff` file is included in the supplied folder.

## Suggested reproducibility record

For every published run, archive at least the following:

```matlab
run_info.matlab_release = version('-release');
run_info.lambda_grid = lam;
run_info.stlsq_iterations = N_iter;
run_info.normalize_library = NormalizeLib;
run_info.noise = Noise;
run_info.rng = rng;
run_info.sample_time_s = ops.h;
run_info.initial_states = [x0, x0_test, x0_val];

save(fullfile('Results','run_info.mat'), 'run_info', 'ODE_Best', ...
     'Score', 'dScore_all_Best', 'LHS_Sym_all', 'SINDy_Struct_all')
```

Also record the Git commit hash and the exact data-file hashes used for the run.

## Data and model limitations

- The state trajectories used for identification, model selection, and validation are generated by the nominal greenhouse model; they are not measured greenhouse state trajectories.
- Outdoor disturbances originate from measured weather, while control inputs are precomputed MPC trajectories.
- Splitting consecutive time intervals provides temporal holdout validation, but it is not by itself cross-scenario validation under different years, greenhouse configurations, crop conditions, or control strategies.
- The identified models for dry matter, CO2, and humidity are low-order approximations of dynamics containing rational and exponential terms. Their coefficients should not be interpreted as direct estimates of nominal physical parameters unless algebraic equivalence and identifiability are established.
- No closed-loop control validation is provided by this code.

## Citation

If you use this code, cite the accompanying manuscript and the original SINDy and SINDy-PI methods. Add the final journal citation and DOI here after publication.

Public repository URL:

<https://github.com/xcheng20/SINDyPI-PIML-Lettuce-Greenhouse>
