"""
Debinding Gas-Pressure Model (Slit-Based Transport)
Author: Your Name
License: MIT

This script computes:
 - Binder conversion α(t) via Coats–Redfern kinetic fitting
 - Gas generation rate dα/dt
 - Internal pressure p_g(t) using a slit-based interfacial transport model
 - Output CSV with all time-history variables

INPUT CSV FORMAT (TGA dataset)
-----------------------------------
The input CSV must contain the following columns:
 - "Temperature (°C)" : Temperature history in Celsius
 - "Time (min)"       : Elapsed time in minutes
 - "Weight (%)"       : Relative sample mass (0–100%)

Rows must be ordered arbitrarily; the script internally sorts by time.
All values must be numeric.
"""

import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt

# ============================================================
# 0. USER PARAMETERS (edit these for your dataset)
# ============================================================

INPUT_CSV  = "example_TGA.csv"      # Input TGA file (replace with your path)
OUTPUT_CSV = "output_results.csv" # Output results file

# Select Coats–Redfern model: one of {"D1", "D2", "D3", "F1"}
MODEL = "F1"

# Geometry of interfacial binder layer
R_layer = 12.5e-3      # [m] specimen radius
h_layer = 0.5e-6       # [m] interfacial binder thickness

# Binder + gas properties
rho_binder = 1.27 * 1000.0   # [kg/m^3] binder density
M_gas = 44e-3                # [kg/mol] CO2 molecular weight
R_gas = 8.314                # [J/mol/K]
P0 = 101325.0                # [Pa] ambient pressure

# Initial void fraction in the interfacial control volume
eps0_void = 0.02             # Example: 2% initial porosity/roughness volume

# ============================================================
# INTERNAL CONSTANTS (fixed model parameters, not user-editable)
# ============================================================

# Minimum transport porosity used only to avoid L_flow = 0 at α=0.
# A tiny lower bound (1e-3) stabilizes the flow length and prevents division by zero.
_ALPHA_MIN_TRANSPORT = 1e-3

# Very small numerical epsilon used for clipping α to avoid log(0) or division by zero.
_EPS_ALPHA = 1e-12


# ============================================================
# 1. LOAD TGA DATA
# ============================================================

df_raw = pd.read_csv(INPUT_CSV)

T_C_raw      = pd.to_numeric(df_raw["Temperature (°C)"], errors="coerce").to_numpy()
time_min_raw = pd.to_numeric(df_raw["Time (min)"], errors="coerce").to_numpy()
W_raw        = pd.to_numeric(df_raw["Weight (%)"], errors="coerce").to_numpy()

mask = np.isfinite(T_C_raw) & np.isfinite(time_min_raw) & np.isfinite(W_raw)
T_C_raw      = T_C_raw[mask]
time_min_raw = time_min_raw[mask]
W_raw        = W_raw[mask]

# Convert units
T_K_raw = T_C_raw + 273.15
t_s_raw = time_min_raw * 60.0


# ============================================================
# 2. MODEL FUNCTIONS: f(α), g(α) FOR COATS–REDFERN FITTING
# ============================================================

def get_model_funcs(model: str):
    """
    Returns model-specific functions f(α) and g(α) used in Coats–Redfern fitting.
    These match standard solid-state kinetic models.
    """
    eps = _EPS_ALPHA
    def clip(a):
        return np.clip(a, eps, 1.0 - eps)

    m = model.upper()

    if m == "D1":
        # One-dimensional diffusion model
        def f(a):  return 0.5 / clip(a)
        def g(a):  return clip(a)**2

    elif m == "D2":
        # Jander diffusion model
        def f(a):
            a = clip(a)
            return 2 * (1-a)**(2/3) / (1 - (1-a)**(1/3))
        def g(a):
            a = clip(a)
            return (1 - (1-a)**(1/3))**2

    elif m == "D3":
        # Crank diffusion model
        def f(a):
            a = clip(a)
            return 1.5 / ((1-a)**(-1/3) - 1)
        def g(a):
            a = clip(a)
            return 1 - (2/3)*a - (1-a)**(2/3)

    elif m == "F1":
        # First-order reaction model
        def f(a):
            return 1 - clip(a)
        def g(a):
            a = clip(a)
            return -np.log(1 - a)

    else:
        raise ValueError(f"Unknown model: {model}")

    return np.vectorize(f), np.vectorize(g)

f_alpha, g_alpha = get_model_funcs(MODEL)


# ============================================================
# 3. COATS–REDFERN LINEAR FIT FOR A AND E
# ============================================================

R = R_gas

order_T = np.argsort(T_K_raw)
T_K_CR = T_K_raw[order_T]
W_CR   = W_raw[order_T]

N = len(W_CR)
k_end = max(5, N//20)     # Use 5% of data at each end for stable W0, W∞ averaging
W0    = W_CR[:k_end].mean()
W_inf = W_CR[-k_end:].mean()

denom = W0 - W_inf
alpha_CR = (W0 - W_CR) / denom
alpha_CR = np.clip(alpha_CR, 0.0, 0.999999)

# Fit only the mid-range (0.1 < α < 0.9)
mask_fit = (alpha_CR > 0.1) & (alpha_CR < 0.9)
T_fit = T_K_CR[mask_fit]
alpha_fit = alpha_CR[mask_fit]

g_fit = g_alpha(alpha_fit)
g_fit = np.where(g_fit <= 0, 1e-12, g_fit)

y = np.log(g_fit / (T_fit**2))
x = 1.0 / T_fit

A_mat = np.vstack([np.ones_like(x), x]).T
coeffs, *_ = np.linalg.lstsq(A_mat, y, rcond=None)
a, b = coeffs

# Heating rate (K/min) — This is a user-independent assumption because CR formula requires β.
# Replace with your true TGA heating rate if known.
beta = 2.0

E = -b * R
A_pre = (beta * E / R) * math.exp(a)

# ============================================================
# 4. COMPUTE α(t) USING ACTUAL TEMPERATURE PROFILE
# ============================================================

order_t = np.argsort(t_s_raw)
T_C = T_C_raw[order_t]
T_K = T_K_raw[order_t]
t_s = t_s_raw[order_t]

A_s = A_pre / 60.0   # convert A from 1/min → 1/s

alpha_kin = np.zeros_like(t_s)
alpha_kin[0] = 1e-4       # Small initial α for numerical stability

stop_index = len(t_s)-1

for i in range(len(t_s)-1):
    dt = max(t_s[i+1] - t_s[i], 1e-6)
    T_i = T_K[i]
    k_i = A_s * math.exp(-E/(R*T_i))

    a_i = max(alpha_kin[i], 1e-4)
    dalpha_dt_i = f_alpha(a_i) * k_i

    new_alpha = alpha_kin[i] + dalpha_dt_i*dt

    if new_alpha >= 0.999:
        alpha_kin[i+1:] = 0.999
        stop_index = i+1
        break
    alpha_kin[i+1] = new_alpha

# Compute kinetic generation rate dα/dt(t)
dalpha_dt_kin = np.zeros_like(t_s)
for i in range(len(t_s)):
    if i >= stop_index:
        continue
    T_i = T_K[i]
    k_i = A_s * math.exp(-E/(R*T_i))
    a_i = max(alpha_kin[i], 1e-4)
    dalpha_dt_kin[i] = f_alpha(a_i) * k_i


# ============================================================
# 5. PRESSURE MODEL — SLIT TRANSPORT
# ============================================================

V_layer_tot = math.pi * R_layer**2 * h_layer
V_void0   = eps0_void * V_layer_tot
V_binder0 = (1 - eps0_void) * V_layer_tot

m0    = rho_binder * V_binder0
n_max = m0 / M_gas      # Max possible CO2 moles if α=1

def Vp_of_alpha(a):
    """Pore volume: existing voids + gas-filled former binder volume."""
    a = float(np.clip(a, 0.0, 0.999999))
    return V_void0 + a * V_binder0

# ===== Gas viscosity model (Sutherland law) =====
def mu_CO2_sutherland(T_K_val):
    """
    CO₂ dynamic viscosity using Sutherland's correlation.
    Constants are typical literature values, not user-specific.
    """
    mu0 = 1.37e-5   # Pa·s at 273 K
    T0  = 273.0
    S   = 222.0
    return mu0 * (T_K_val/T0)**1.5 * (T0+S)/(T_K_val+S)

def L_flow_alpha(a):
    """
    Flow length grows from edge inward as debinding progresses.
    To avoid zero, enforce a ≥ α_min_transport.
    """
    a_eff = max(a, _ALPHA_MIN_TRANSPORT)
    return a_eff * R_layer

def tau_slit(a, T_K_val):
    """
    Characteristic diffusion/venting timescale inside a slit:
       g = a*h_layer       : slit gap
       k = g²/12           : permeability of parallel-plate slit
       D_eff = k * P0 / μ
       τ = L_flow² / D_eff
    """
    a_eff = max(a, _ALPHA_MIN_TRANSPORT)
    g = a_eff * h_layer
    k = g**2 / 12.0
    L_eff = L_flow_alpha(a_eff)
    mu = mu_CO2_sutherland(T_K_val)

    D_eff = k * P0 / mu
    if D_eff <= 0:
        return 1e12
    return L_eff**2 / D_eff


# ===== Integrate pressure =====
n = np.zeros_like(t_s)
P = np.zeros_like(t_s)

a0 = max(alpha_kin[0], 1e-4)
Vp0 = Vp_of_alpha(a0)
n_eq0 = P0 * Vp0 / (R_gas * T_K[0])

n[0] = n_eq0
P[0] = P0

for i in range(len(t_s)-1):
    dt = max(t_s[i+1] - t_s[i], 1e-6)

    a_i  = max(alpha_kin[i], 1e-4)
    Vp_i = Vp_of_alpha(a_i)
    n_eq_i = P0 * Vp_i / (R_gas*T_K[i])

    tau_i = tau_slit(alpha_kin[i], T_K[i])
    dn_gen_dt = n_max * dalpha_dt_kin[i]

    # Analytic update of: dn/dt = generation – venting
    if tau_i < 1e-9:
        n_new = n_eq_i
    else:
        expfac = math.exp(-dt/tau_i)
        n_new = n_eq_i + (n[i] - n_eq_i)*expfac + dn_gen_dt*tau_i*(1-expfac)

    n[i+1] = max(n_new, 0.0)

    a_ip1 = max(alpha_kin[i+1], 1e-4)
    P[i+1] = n[i+1] * R_gas * T_K[i+1] / Vp_of_alpha(a_ip1)

P_over_P0 = P / P0


# ============================================================
# 6. SAVE RESULTS
# ============================================================

out_df = pd.DataFrame({
    "time_min": t_s/60,
    "T_C": T_C,
    "T_K": T_K,
    "alpha_kinetic": alpha_kin,
    "dalpha_dt_kin": dalpha_dt_kin,
    "P_over_P0": P_over_P0,
    "P_Pa": P
})
out_df.to_csv(OUTPUT_CSV, index=False)
print(f"Saved results to {OUTPUT_CSV}")


# ============================================================
# 7. OPTIONAL PLOTS
# ============================================================

plt.figure()
plt.plot(T_C, P_over_P0)
plt.xlabel("Temperature (°C)")
plt.ylabel("P/P0")
plt.grid(True)
plt.title("Pressure ratio vs Temperature")

plt.figure()
plt.plot(T_C, alpha_kin)
plt.xlabel("Temperature (°C)")
plt.ylabel("Alpha (kinetic)")
plt.grid(True)
plt.title("Kinetic conversion vs Temperature")

plt.show()
