# phytometr 0.0.0.9000

## Initial Development Release – Host–Environment–Holobiont Power Framework

This is the first functional prototype of **phytometr**, providing a complete 
simulation and power-analysis system for clonal host experiments and 
environmentally driven holobiont responses (Phase 1 of package development).

### New Core Features

#### Environmental Simulation Framework
* Added **three independent modes** for environmental definition:
  - `define_env_space()` – gradient-based environmental ranges for experiments  
  - `define_env_locations()` – user-specified site/condition definitions  
  - `define_env_grid()` – grid of possible experimental conditions for factorial or multivariate setups  
* All modes support:
  - Environmental shapes (`linear`, `nonlinear`, `random`, `patchy`)  
  - Realistic **within-group variability**  
  - Environmental **precision control** via `target_sd`  
* Added `simulate_individual_env()` to expand group-level environments to individual-level conditions.

#### Host Model Framework
* Added `define_host_model()` with two modes:
  - `"clone"` – no host genetic variance  
  - `"natural"` – placeholder structure for Phase 2 (population genetics model)  
* Stores host variance, population structure metadata, and simulation notes.

#### Holobiont Component Specification & Simulation
* Added `simulate_holobiont_components()` to simulate:
  - User-defined holobiont components (bacteria, fungi, herbivores, etc.)
  - Environment-driven effect sizes (`r2_env`)
  - Optional negative effects (e.g., stressors)
  - Group-level component means and variability
* Added `simulate_env_response()` to generate individual-level holobiont responses combining:
  - Environmental effects  
  - Host model variance  
  - Residual variance  

#### Power Analysis Engine
* Added `simulate_holobiont_power()`:
  - Computes power across combinations of:
    - Number of groups  
    - Individuals per group  
  - Automatically uses all defined groups in `env_space`, `env_locations`, or `env_grid`
  - Estimates component-wise and overall environmental response power
  - Returns tidy power table for downstream analysis

### New Visualisation Suite

#### Environmental Plots
* `plot_env_space()` – planned vs. realised environmental gradients  
* `plot_env_individual()` – individual-level environmental distributions  

#### Holobiont Diagnostics
* `plot_holobiont_components()` –  
  - PCA of components (with group-based colouring and convex hulls)  
  - Component correlation heatmap (Pearson r with p-values)  
  - Environment–component correlation heatmap (Pearson r with p-values)

#### Power & Variance Plots
* `plot_power_landscape()` – heatmap of sample-size × group-number power  
* `plot_variance_decomposition()` – variance partitioning for each holobiont component  

#### Design Diagnostics
* `diagnose_phytomet_design()` – quick diagnostics summary  
* `plot_phytomet_diagnostics()` – publication-ready diagnostic visualisation

### Utility Improvements
* Fully consistent ID + group tracking across all simulation steps  
* Harmonised scaling of effect sizes  
* Better defaults for gradient noise  
* Standardised return objects (`env_space`, `host_model`, `holo_components`)  
* Harmonised colour scales and aesthetics for all heatmaps

### Notes
* API may change slightly as functionality stabilises.

---

