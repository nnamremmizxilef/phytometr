# phytometr 0.0.0.9001

## Development Changes
* Added full Phase 1 framework for clonal hosts, environmental simulations, and holobiont-response power analysis.
* Added environmental gradient comparison:
  - `define_env_gradient_range()`
  - `define_env_locations()`
  - `generate_env_grid()`
* Added individual-level environmental simulation:
  - `simulate_individual_env()`
* Added system-wide holobiont response simulation:
  - `simulate_env_response()`
  - `simulate_holobiont_components()`
* Added full N × environmental-range power analysis:
  - `simulate_holobiont_power()`
  - `plot_power_landscape()`
  - `plot_variance_decomposition()`
* Added extensive diagnostic tools:
  - `plot_env_space()`
  - `plot_env_individual()`
  - `plot_holobiont_components()`
  - `plot_phytomet_diagnostics()`
* Improved heatmaps (correlations + env–holobiont):  
  - better colour scales (–1 to 1 diverging)  
  - r-values printed on tiles  
  - p-values included  
* PCA now uses grouped colours, point ellipses, and cluster outlines.

---

# phytometr 0.0.0.9000

## Initial Development Version
* Added basic environmental space simulation (`define_env_space()`).
* Added early clonal-host modelling (`define_host_model()`).
* Added initial holobiont component simulator.
* Added prototype power analysis for simple designs.
* Set up project structure, documentation, and GitHub integration.
