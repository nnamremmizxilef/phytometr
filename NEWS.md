# phytometR 0.0.0.9000

## Initial Development Release — Core Framework for Host–Environment–Holobiont Simulation

This is the very first development version of **phytometR**, introducing a unified framework
to simulate environmental gradients, host variation, holobiont responses, and experimental
power for both **clonal systems (Phase 1)** and **natural populations (Phase 2)**.

The package is designed to support experimental planning in ecology, microbiome research,
plant–microbe interactions, stress experiments, and host–holobiont biology.

---

## Key Features Introduced

### 1. Environmental Gradient Simulation
**Functions:**
- `define_env_gradient()`  
- `simulate_environmental_conditions()`

Supports:
- continuous, categorical, or multivariate environmental variables  
- spatial vs. controlled experimental gradients  
- controllable environmental noise (field vs. greenhouse precision)  
- user-defined environmental ranges based on literature or pilot data  

This establishes the foundation for asking:
> “How much environmental variation is needed to detect an effect?”

---

### 2. Host System Modeling (Clone vs. Natural Population)
**Function:**
- `define_host_model()`

Two modes:
- **Clonal host system (Phase 1)**  
  No genetic variation; all individuals share one genotype.

- **Natural population (Phase 2)**  
  With optional genetic variance and population structure:
  - panmictic  
  - isolation-by-distance  
  - discrete populations  
  - structured allele frequency variance  

This enables experiments across the continuum from controlled clonal setups to landscape-level sampling.

---

### 3. Holobiont Component Simulation
**Function:**
- `simulate_holobiont_components()`

Allows simulation of:
- microbiome diversity  
- microbial functional traits  
- herbivore/pathogen load  
- community composition proxies (e.g., PC1, richness, evenness)  

Each component can have:
- environmental effect size  
- host genetic contribution  
- environmental noise  
- optional G×E interaction  

Users can define as many components as needed, reflecting literature-based effect sizes.

---

### 4. Experimental Power Simulation
**Function:**
- `simulate_power()` (initial version)

Allows exploration of:
- how sample size interacts with environmental differences  
- trade-offs between replication, environmental range, and effect detectability  
- differences between clonal vs. natural host systems  

This is the core engine to support planning field trials or controlled experiments.

---

### 5. Diagnostic Tools (Initial Version)
**Function:**
- `plot_diagnostics()`

Provides visual summaries of:
- environmental variance  
- host genetic variance (if applicable)  
- holobiont response distributions  
- signal-to-noise structure  

These diagnostics guide users toward realistic and statistically supported experimental designs.

---

## Goals of This Release

This initial version establishes the **foundation** for the package:

- Simulation-first workflow  
- Modular design (environment → host → holobiont → power)  
- Flexible enough for any holobiont system  
- Supports both field gradients and controlled experiments  
- Prepares the package for later integration of:
  - real data ingestion  
  - more advanced power metrics  
  - genetic coalescent simulations (Phase 2 extensions)

---

## Functions Included in 0.0.0.9000

| Category | Functions |
|---------|-----------|
| Environment simulation | `define_env_gradient()`, `simulate_environmental_conditions()` |
| Host modeling | `define_host_model()` |
| Holobiont responses | `simulate_holobiont_components()` |
| Power analysis | `simulate_power()` |
| Diagnostics | `plot_diagnostics()` |

---

## Notes

This is a development-stage release. API changes are expected before the first stable version.

---

