# Debinding Gas-Pressure Model (Slit-Based Transport)

This repository provides a physics-based Python implementation for estimating  
**internal gas pressure during polymer-binder debinding** in ceramic additive manufacturing (VPP/SLA/DLP).  
The model computes:

- Binder conversion α(t) using **Coats–Redfern** kinetic fitting  
- Gas generation rate dα/dt  
- Slit-based interfacial transport and effective venting time scale τ  
- Internal pressure evolution p_g(t) and normalized pressure p_g / p₀  
- Visualization and CSV export for downstream analysis

This code is intended for researchers studying debinding defects  
(delamination, blistering, Z-cracks, XY-surface cracks, binder modulus effects).

---

## ✨ Key Features

### ✔ Coats–Redfern kinetic fitting  
Supports solid-state reaction models:
- D1 (1D diffusion)
- D2 (Jander)
- D3 (Crank)
- F1 (First-order reaction)

Automatically estimates activation energy **E** and pre-exponential factor **A**.

---

### ✔ Slit-based venting model for interfacial binder layers  
The pressure model assumes gas escapes primarily through a **developing interfacial slit**,  
representing a delamination gap or transport channel:

- Gap height grows: `g(α) = α · h_layer`
- Slit permeability: `k = g² / 12`
- Effective venting diffusivity: `D_eff = k P₀ / μ(T)`
- Flow length increases as debinding progresses: `L_flow = max(α, α_min) · R_layer`
- Venting timescale: `τ = L_flow² / D_eff`

This formulation captures early-stage gas accumulation due to:
- small initial permeability,
- strong viscosity dependence μ(T),
- limited transport length at low α.

---

## 📂 Repository Structure

```
.
├── debinding_model.py        # Main Python script
├── example_data/             # Folder for example input CSV (user-provided)
├── output/                   # Generated time-history results (optional)
└── README.md
```

---

## 📄 Input CSV Format (TGA Data)

Your input CSV must contain the following **exact column names**:

| Column name         | Description                              | Example |
|---------------------|-------------------------------------------|---------|
| `Temperature (°C)`  | Temperature history in Celsius            | 25 → 600 |
| `Time (min)`        | Timestamp in minutes                      | 0 → 60 |
| `Weight (%)`        | Remaining mass (%) from TGA               | 100 → 0 |

**Rows can be in any order** — the script sorts internally by time and temperature.

Place the file in:

```
example_data/your_file.csv
```

Then set in the Python script:

```python
INPUT_CSV = "example_data/your_file.csv"
```

---

## ▶️ How to Run

### 1️⃣ Install dependencies

```bash
pip install numpy pandas matplotlib
```

### 2️⃣ Run the model

```bash
python debinding_model.py
```

### 3️⃣ Output

A CSV file (e.g., `output_results.csv`) containing:

- time_min  
- T_C  
- T_K  
- alpha_kinetic  
- dalpha_dt_kin  
- P_over_P0  
- P_Pa  

---

## 📊 Example Plots

- **Pressure vs. Temperature (P/P0)**
- **Binder conversion α vs. Temperature**

Plots help identify:
- gas accumulation stages,
- early pressure peaks (before DTG maximum),
- delamination risk conditions,
- transport-limited vs kinetics-limited regimes.

---

## 📘 Model Assumptions & Notes

- Ideal gas behavior (pV = nRT)  
- CO₂ as dominant pyrolysis product  
- Slit-based permeability (parallel-plate approximation)  
- Debinding proceeds from the outer radius inward  
- α_min (10⁻³) prevents division-by-zero early in the simulation  
- Viscosity via Sutherland's law (CO₂ parameters)

These constants are scientific model parameters,  
**not user settings** and must remain unchanged.

---

## 📎 Citation

If you use or modify this code, please cite:

```
Kim, Y. & Collaborators (2025).
Slit-Based Gas Pressure Modeling for Binder Debinding in Ceramic Additive Manufacturing.
GitHub Repository.
```

---

## 📬 Contact

If you have questions or wish to contribute improvements,  
please open an issue or submit a pull request.

---

# ❓ FAQ

### Q. *Can I upload my own example dataset?*  
Yes — simply place your CSV file into:

```
example_data/
```

and point the script to it via:

```python
INPUT_CSV = "example_data/your_dataset.csv"
```

---

