# phytometR

<div style="display: flex; align-items: center;">
  <img src="inst/hex/logo.png" alt="phytometr logo" height="160" style="margin-right: 25px;">
  <p style="margin: 0;">
    <strong>phytometR</strong> provides simulation-based power analysis and experimental design tools for host–environment–holobiont systems. It simulates host genotypes (clonal or natural populations), environmental gradients (field, controlled, location-based, grid-based), and holobiont responses (diversity, composition, herbivory, etc.) to determine required sample sizes, environmental range, and replication.
  </p>
</div>

---

# Overview

The package supports both:

## Application 1 — **Clonal host systems (common garden / controlled experiments)**
- Negligible host genetic variation  
- Environmental variation drives holobiont response  
- Goal: Determine **number of conditions and replicates** needed  
- Examples: tree clones, grafted plants, model systems  

## Application 2 — **Natural host populations**
- Includes realistic host genomic variation  
- Optional population structure (panmictic, IBD, discrete, admixture)  
- Environment + host + G×E effects  
- Goal: **plan spatial field sampling and holobiont characterization**

phytometR provides:
- environmental gradient simulation  
- location‑based or grid‑based environmental designs  
- host model definition (clone or natural)  
- holobiont component modelling  
- holobiont response simulation  
- power estimation for N × environmental gradient  
- correlation, PCA, and variance decomposition diagnostics  
- integration of real environmental or holobiont data  

---

# Quick Start

## 1. Define environmental conditions

### Field/controlled gradient
```r
env_space <- define_env_space(
  variables = list(
    temperature = list(range = c(15, 30), shape = "linear"),
    moisture    = list(range = c(20, 60), shape = "random")
  ),
  n_groups = 6,
  context = "field"
)
```

### Location-based definition (new)
```r
locations <- data.frame(
  site = c("A","B","C"),
  temperature = c(10, 18, 25),
  moisture = c(35, 55, 80)
)

env_loc <- define_env_locations(
  locations       = locations,
  within_group_sd = c(temperature = 1.2, moisture = 4)
)
```

### Environmental grid (new)
```r
env_grid <- define_env_grid(
  variables = list(
    temperature = seq(10, 25, length.out = 5),
    moisture    = seq(40, 80, length.out = 5)
  )
)
```

---

## 2. Simulate individual-level environmental values

```r
ind_env <- simulate_individual_env(env_space, n_per_group = 20)
```

---

## 3. Define host model

### Clone system
```r
host_clone <- define_host_model("clone", var_host = 0)
```

### Natural population
```r
host_nat <- define_host_model(
  type = "natural",
  pop_structure = "panmictic",
  var_host = 0.25
)
```

---

## 4. Define holobiont components

Each component receives an **R² partition** into environmental drivers.

```r
components <- list(
  root_bacteria = list(r2_env = c(temperature = 0.30)),
  root_fungi    = list(r2_env = c(temperature = 0.25, moisture = 0.10)),
  herbivores    = list(r2_env = c(moisture = 0.40))
)
```

---

## 5. Simulate holobiont components for individuals

```r
components_sim <- simulate_holobiont_components(
  components = components,
  env_data   = ind_env
)
```

What this produces:
- A value per individual per component  
- Values are generated so that correlations with environmental drivers **match the R² you specified**  
- Results include:
  - component values  
  - “true” correlations  
  - variance explained  

---

## 6. Produce diagnostics

### PCA + clustering by group
```r
plot_holobiont_components(components_sim, env = ind_env)
```

This includes:
- PCA with group cluster outlines  
- Component‑component correlation heatmap  
- Environment–component correlation heatmap  
- All heatmaps include  
  - r values  
  - p values  
  - full −1 to 1 diverging colour scales  

---

## 7. Power analysis

```r
power <- simulate_holobiont_power(
  env_space    = env_space,
  host_model   = host_clone,
  components   = components,
  n_groups_use = c(3,5,7),
  n_per_group  = c(10,20,40),
  n_replicates = 200
)
```

---

## 8. Plot design performance

```r
plot_power_landscape(power)
plot_variance_decomposition(power, components)
plot_phytomet_diagnostics(power)
```

---

# How phytometR Calculates Things

## Environmental simulation
- Group means follow your defined ranges  
- Individual values = group mean + noise  
- Noise depends on **within_group_sd**, reflecting:
  - field-year unpredictability  
  - chamber precision  
- Continuous, grid, and location-based designs all supported

## Holobiont simulation
For each component and each environmental variable:

```
y = β_env * X_env  + β_host * X_host + β_resid
```

Where:
- β_env values are chosen such that **R² = effect strength you specified**
- Host contribution is added for natural models
- Residual noise is scaled so total variance = 1

This ensures:
- Achieved correlations ≈ expected R²
- PCA structure is meaningful
- Diagnostics reflect real simulation behaviour

---

# Understanding Power in phytometR

Power = probability of detecting an environment–holobiont relationship under your design.

### Common thresholds:
| Expected biological signal | Suggested minimum power |
|----------------------------|--------------------------|
| Strong relationships (R² > 0.3) | 0.80 – 0.90 |
| Moderate effects (0.1–0.3) | 0.70 – 0.80 |
| Weak ecological effects (~0.05) | 0.50 – 0.60 (if biologically relevant) |

phytometR provides power for:
- **environment-only tests**  
- each holobiont component individually  
- global multivariate tests  

---

# Real Data Integration

```r
my_env <- read.csv("env.csv")
my_holo <- read.csv("holo.csv")

env_space_real <- define_env_locations(my_env)
components_real <- my_holo

plot_holobiont_components(components_real, env = my_env)
```

This allows:
- calibrating simulation to real field data  
- comparing expected vs observed correlations  
- refining design before sampling  

---

# Experimental Design Applications

### Field gradients
- Select number of sites  
- Evaluate if gradient width is sufficient  
- Test sensitivity to environmental noise

### Controlled experiments
- Decide number of chambers / climate settings  
- Quantify precision vs sample size trade‑offs  
- Test multi‑factor designs

### Natural populations
- Include population structure  
- Integrate genomic variation  
- Conduct feasibility checks for G×E analyses  

---

# License

MIT © Felix Zimmermann

---

# Contact

Felix Zimmermann  
WSL – Swiss Federal Institute for Forest, Snow and Landscape Research  
**Email:** felix.zimmermann@wsl.ch  
**GitHub:** https://github.com/nnamremmizxilef/phytometr  
