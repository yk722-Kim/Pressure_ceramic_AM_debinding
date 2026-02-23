# Debinding Gas-Pressure Model for Ceramic VPP

This repository provides a physics-based Python implementation for estimating
the **internal gas pressure evolution** during polymer-binder debinding in ceramic VPP,
using a **TGA-derived kinetic model** and an **interfacial transport (slit) pressure-relaxation model**.

The script computes:
- Binder conversion **α(t)** via **Coats–Redfern** kinetic fitting (F1 model)
- Decomposition rate **dα/dt**
- Interfacial transport timescale **τ_diff(α, T)**
- Internal pressure evolution **p_g(t)** and normalized pressure **p_g/p0**
- CSV export and optional plots for downstream analysis

---

## Key Features

### 1) Coats–Redfern kinetic fitting
This implementation uses the **first-order reaction model**:
- f(α) = 1 − α
- g(α) = −ln(1 − α)

From the Coats–Redfern linearization,
the code estimates the activation energy **E** and pre-exponential factor **A**
using the user-specified heating rate **β (K/min)**.

### 2) Interfacial pressure accumulation and relaxation model
The pressure calculation is carried out at the **polymer-dominated interlayer** scale,
modeled as a thin cylindrical control region with specimen radius **r** and cohesive interlayer thickness **ℓ**:

- Interlayer volume:  V_layer = π r² ℓ
- Initial void volume: V_void0 = ϕ V_layer
- Initial binder volume: V_binder0 = (1 − ϕ) V_layer
- Evolving pore volume: V_p(α) = V_void0 + α V_binder0

Gas generation is driven by the decomposition rate:
- n0 = (ρ_binder V_binder0) / M_gas
- dn/dt |_gen = n0 · dα/dt

Pressure relaxation is modeled through a slit-like pathway that evolves with conversion:
- δ(α) = α ℓ
- k_frac(α) = δ(α)² / 12
- D_p(α, T) = k_frac(α) · p0 / μ(T)
- L(α) = max(α, α_min) · r
- τ_diff(α, T) = L(α)² / D_p(α, T)

The internal gas amount n(t) is updated by balancing generation and venting toward equilibrium:
- n_eq(t) = (p0 V_p(α(t))) / (R T(t))

Finally, pressure is computed using the ideal gas law in the evolving pore volume:
- p_g(t) = n(t) R T(t) / V_p(α(t))

---

## Repository Structure

```
.
├── run.py # Main script
├── example_TGA.csv # Example input TGA dataset
├── output_results.csv # Generated output (created after running)
└── README.md
```


---

## Input CSV Format (TGA Data)

Your input CSV must contain the following **exact column names**:

| Column name         | Description                     |
|---------------------|---------------------------------|
| `Temperature (°C)`  | Temperature history in Celsius  |
| `Time (min)`        | Elapsed time in minutes         |
| `Weight (%)`        | Remaining mass (%) from TGA     |

Rows can be unordered; the script internally sorts by time.

---

## How to Run

### 1) Install dependencies
```bash
pip install numpy pandas matplotlib


### 2) Run

```bash
python run.py
```

### User Parameters to Edit (in run.py)

At the top of run.py, you may edit:

- INPUT_CSV, OUTPUT_CSV
- Geometry:
- R_layer : specimen radius r [m]
- ell_layer: cohesive interlayer thickness ℓ [m]
- Material / gas:
- rho_binder [kg/m³]
- M_gas [kg/mol] (default CO2)
- phi_void initial void fraction ϕ
- Coats–Redfern heating rate:
- beta [K/min]

### 2) Output

A CSV file (e.g., `output_results.csv`) containing:

- time_min  
- T_C  
- T_K  
- alpha_kinetic  
- dalpha_dt_kin  
- P_over_P0  
- P_Pa  

---

## Example Plots

- **Pressure vs. Temperature (P/P0)**
- **Binder conversion α vs. Temperature**

Plots help identify:
- gas accumulation stages,
- early pressure peaks,
- delamination risk conditions,
- transport-limited vs kinetics-limited regimes.

---

## Model Assumptions & Notes

- Ideal gas behavior: pV = nRT
- CO₂ is used as a representative dominant gaseous product (via M_gas)
- Interfacial transport uses a parallel-plate permeability approximation: k_frac = δ²/12
- Gas viscosity μ(T) is computed using a Sutherland-type correlation for CO₂
- α_min = 1e−3 is used only to avoid L(α) = 0 at α → 0 for numerical stability

These constants are scientific model parameters,  
**not user settings** and must remain unchanged.

---

## Citation

If you use or modify this code, please cite:

```
Kim, Yuseok & Sobhani, Sadaf. Additive Manufacturing, 2026, 105130, https://doi.org/10.1016/j.addma.2026.105130
Mechanistic Insights into Debinding-induced Defects in VPP-printed Ceramics
GitHub Repository.
```

---

