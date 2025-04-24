# Airlift Scheduling Problem Using 2-Stage Linear Programming Methods ✈️ 

This project implements and compares **Fat** and **Slim** Benders Decomposition approaches for a two-stage stochastic airlift scheduling problem. It includes full support for **SMPS file generation**, scenario sampling via **log-normal distributions**, and scalable experimentation on small and large instances.

---

## Problem Description

We study a two-stage stochastic linear program (SLP) for airlift scheduling:

- **First Stage:** Decide how many flight hours to allocate from each aircraft to each route.
- **Second Stage (Recourse):** Adjust allocations and slack to meet uncertain demand via scenario-specific reassignments.

The model incorporates:
- Aircraft capacity constraints
- Route-specific demand (uncertain)
- Aircraft re-routing and associated switching costs
- Penalties for unmet demand

---

## Methods

### Fat Benders
Solves the full extensive form for all scenarios in a single model.
- Fast and accurate on small instances
- Memory-intensive for large scenario counts

### Slim Benders
Decomposes the problem into a master and scenario-wise subproblems.
- Adds **aggregated Benders cuts** (not per-scenario!)
- Scales to large problem sizes with many scenarios
- Supports **feasibility and optimality cuts**

---

## Experiments

We include tests on both:
- `.sto.first`: Jointly sampled log-normal demand (fully correlated)

Performance is reported across:
- Runtime
- Convergence speed (iterations)
- Objective optimality
- Cuts added

We did not run tests using `.sto.second` due to the memory deficiency on Gavin's brick of a computer.

## Generators

Use `write_smps(...)` to generate synthetic SMPS instances:

- Configurable number of aircraft (I), routes (J), and scenarios (S)
- Supports both **joint** and **independent** log-normal scenario generation
- Automatically writes `.cor`, `.tim`, and `.sto` files

---

## Repository Structure

```
.
├── benders.py                # Fat + Slim Benders implementations
├── solver.py                 # Full extensive form solver
├── solve_generated_instances.py  # Script to run & log generated instances
├── solve_provided_instances.py   # Script to run & log provided instances
├── generator.py              # SMPS generator + log-normal sampling
├── reader.py                 # SMPS file parser + matrix builder
├── airlift/                  # Provided instance files (.cor, .tim, .sto)
├── requirements.txt          # Dependency text
└── README.md                 # Project overview and documentation
```

## Dependencies

run pip install -r requirements.txt

## Authors

- Gavin Daves, Rice University
- Beck Edwards, Rice University
- Matthew Cihlar, Rice University
- A Little Help From Our Friend/TA ChatGPT 🤖
