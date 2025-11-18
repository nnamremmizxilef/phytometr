# phytometR

<div style="display: flex; align-items: center;">
  <img src="inst/hex/logo.png" alt="phytometr logo" height="160" style="margin-right: 25px;">
  <p style="margin: 0;">
    <strong>phytometR</strong> provides simulation-based power analysis and experimental design tools for host–environment–holobiont systems. It simulates host genotypes (clonal or natural populations), environmental gradients (field or controlled), and holobiont responses (diversity, composition, herbivory, etc.) to determine required sample sizes, environmental range, and replication.
  </p>
</div>

---

## Overview

The package supports both:

### Phase 1 — **Clonal host systems (common garden / controlled experiments)**
- Host genetic variation is negligible or absent  
- Environmental variation is the main driver  
- Goal: Determine number of conditions and replicates needed  
- Use cases: tree clones, grafted plants, model organisms  

### Phase 2 — **Natural host populations**
- Realistic host genomic variation  
- Optional population structure (panmictic, IBD, discrete, admixture)  
- Environment + host + G×E effects on holobiont  
- Goal: Plan spatial field sampling and holobiont characterization  

phytometR allows:
- realistic simulation of environmental variables  
- host model definition (clone vs. natural)  
- holobiont component simulation  
- response generation  
- power estimation across N × environmental gradient combinations  
- diagnostics for gradient quality, response sensitivity, and design performance  
- optional integration of real environmental or holobiont data  

---

# Quick Start

## 1. Define environmental gradients

    env_grad <- define_env_gradient(
      variable = "temperature",
      range = c(15, 30),
      shape = "linear"
    )

    env_space <- define_env_space(
      variables = list(
        temperature = list(range = c(15, 30), shape = "linear"),
        moisture    = list(range = c(20, 60), shape = "random")
      ),
      n_groups = 6,
      context = "field"
    )

## 2. Choose host model

Clonal host system (Phase 1)

    host_clone <- define_host_model("clone")

Natural population (Phase 2)

    host_nat <- define_host_model(
      type = "natural",
      pop_structure = "isolation_by_distance",
      var_host = 0.15
    )

## 3. Define holobiont components

    components <- simulate_holobiont_components(
      n_components = 3,
      names = c("fungal_diversity", "bacterial_PC1", "herbivore_load"),
      env_effects = c(0.25, 0.30, 0.15),
      host_effects = c(0.10, 0.05, 0.00),
      gxe_effects = c(0.05, 0.10, 0.00)
    )

## 4. Simulate individual environmental conditions

    ind_env <- simulate_env_conditions(env_space, n_per_group = 20)

## 5. Run power analysis

    power <- simulate_holobiont_power(
      env_space = env_space,
      host_model = host_clone,
      components = components,
      n_groups_use = c(4, 6, 8),
      n_per_group = c(10, 20, 30),
      n_replicates = 200
    )

## 6. Visualize power results

    plot_power_landscape(power)
    summary(power)

---

# Diagnostics

Holobiont PCA & correlations

    plot_holobiont_components(
      components_df = simulated_components,
      env = ind_env
    )

Environmental distributions

    plot_env_individual(ind_env)
    plot_env_space(env_space)

Variance decomposition

    plot_variance_decomposition(power, components)

---

# Real Data Integration

    my_env <- read.csv("environment.csv")
    my_holo <- read.csv("holobiont_data.csv")

    env_real <- define_env_space(
      variables = list(
        temperature = list(range = range(my_env$temp), shape = "real")
      ),
      n_groups = length(unique(my_env$site_id)),
      context = "field"
    )

    env_real$realised <- my_env[, c("site_id", "temp")]

    plot_holobiont_components(
      components_df = my_holo[, c("fungi_div", "bact_PC1")],
      env = my_env$temp
    )

---

# Experimental Design Applications

### Field Gradient Studies:
- Number of sites  
- Spacing of sites  
- Required environmental gradient width  
- Optimal number of replicates per site  

### Controlled Experiments:
- Number of experimental conditions  
- Environmental precision (e.g., ±0.2°C vs ±2°C)  
- Replication needed per condition  
- Designs balancing N vs gradient strength  

### Clonal Systems:
- Testing environmental effects on holobiont diversity/composition  
- Finding minimal N needed for reliable response detection  
- Designing multi-factor experiments  

### Natural Populations:
- Identify environmental variation needed to detect signals  
- Model population structure  
- Account for genetic, environmental, and G×E contributions  

---

# Features

- Simulation of multi-dimensional environmental gradients  
- Host genotype simulation (clone or natural)  
- Flexible holobiont modules (root fungi, bacteria, herbivores, pathogens…)  
- User-defined effect sizes from literature  
- Power surface estimation for N × gradient width  
- Diagnostics for environmental noise, collinearity, and response sensitivity  
- Real data calibration to refine simulations  

---

# Planned Extensions

- Integration of spatial coordinates and real GIS gradients  
- Modular addition of holobiont domains (bacteria/mycorrhiza/herbivores)  
- Link functions for abundance, composition, and functional diversity  
- Hierarchical experimental designs  

---

# License

MIT © Felix Zimmermann

# Contact

- Author: Felix Zimmermann  
- Institution: WSL – Swiss Federal Institute for Forest, Snow and Landscape Research  
- Email: felix.zimmermann@wsl.ch  
- GitHub: https://github.com/nnamremmizxilef  
- Issues: https://github.com/nnamremmizxilef/phytometr/issues
