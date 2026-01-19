"""
Debinding Gas-Pressure Model
Author: Yuseok Kim and Sadaf Sobhani
License: MIT

This script computes:
 - Binder conversion α(t) via Coats–Redfern kinetic fitting
 - Gas generation rate dα/dt
 - Internal pressure p_g(t) using an interfacial transport model
 - Output CSV with all time-history variables

INPUT CSV FORMAT (TGA dataset)
-----------------------------------
The input CSV must contain the following columns:
 - "Temperature (°C)" : Temperature history in Celsius
 - "Time (min)"       : Elapsed time in minutes
 - "Weight (%)"       : Relative sample mass (0–100%)

Rows can be unordered; the script internally sorts by time.
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
OUTPUT_CSV = "output_results.csv"   # Output results file

# Geometry of interfacial binder layer
R_layer  = 12.5e-3      # [m] specimen radius
ell_layer = 0.5e-6      # [m] cohesive interlayer thickness ℓ

# Binder + gas properties
rho_binder = 1.27 * 1000.0   # [kg/m^3] binder density
M_gas = 44e-3                # [kg/mol] CO2 molecular weight
R_gas = 8.314                # [J/mol/K]
p0 = 101325.0                # [Pa] ambient pressure  (CHANGED: P0 -> p0)

# Initial void fraction (ϕ) in the interfacial control volume
phi_void = 0.02              # Example: 2% initial porosity/roughness volume

# Heating rate used in Coats–Redfern fitting (K/min)
beta = 2.0

# ============================================================
# INTERNAL CONSTANTS (fixed model parameters, not user-editable)
# ============================================================

# Minimum transport conversion used only to avoid L = 0 at α=0.
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
# 2. MODEL FUNCTIONS: First-order reaction
# ============================================================

def _clip_alpha(a):
    return np.clip(a, _EPS_ALPHA, 1.0 - _EPS_ALPHA)

def f_alpha(a):
    """First-order reaction model: f(α) = 1 - α"""
    return 1.0 - _clip_alpha(a)

def g_alpha(a):
    """Integral form for Coats–Redfern: g(α) = -ln(1 - α)"""
    a = _clip_alpha(a)
    return -np.log(1.0 - a)


# ============================================================
# 3. COATS–REDFERN LINEAR FIT FOR A AND E
# ============================================================

R = R_gas

order_T = np.argsort(T_K_raw)
T_K_CR = T_K_raw[order_T]
W_CR   = W_raw[order_T]

N = len(W_CR)
k_end = max(5, N // 20)   # ~5% of data at each end (min 5 points)
W0    = W_CR[:k_end].mean()
W_inf = W_CR[-k_end:].mean()

denom = W0 - W_inf
if abs(denom) < 1e-12:
    raise ValueError("W0 - W_inf is too small; cannot normalize α from weight history.")

alpha_CR = (W0 - W_CR) / denom
alpha_CR = np.clip(alpha_CR, 0.0, 1.0 - 1e-12)

# Fit only the mid-range (0.1 < α < 0.9)
mask_fit = (alpha_CR > 0.1) & (alpha_CR < 0.9)
T_fit = T_K_CR[mask_fit]
alpha_fit = alpha_CR[mask_fit]

g_fit = g_alpha(alpha_fit)
g_fit = np.where(g_fit <= 0, 1e-12, g_fit)

# Coats–Redfern linearization: y = ln(g(α)/T^2), x = 1/T
y = np.log(g_fit / (T_fit**2))
x = 1.0 / T_fit

A_mat = np.vstack([np.ones_like(x), x]).T
coeffs, *_ = np.linalg.lstsq(A_mat, y, rcond=None)
a, b = coeffs

E = -b * R
A_pre = (beta * E / R) * math.exp(a)   # A in 1/min (consistent with beta in K/min)

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

stop_index = len(t_s) - 1

for i in range(len(t_s) - 1):
    dt = max(t_s[i + 1] - t_s[i], 1e-6)
    T_i = T_K[i]
    k_i = A_s * math.exp(-E / (R * T_i))

    a_i = max(alpha_kin[i], 1e-4)
    dalpha_dt_i = float(f_alpha(a_i)) * k_i

    new_alpha = alpha_kin[i] + dalpha_dt_i * dt

    if new_alpha >= 0.999:
        alpha_kin[i + 1:] = 0.999
        stop_index = i + 1
        break

    alpha_kin[i + 1] = new_alpha

# Compute kinetic generation rate dα/dt(t)
dalpha_dt_kin = np.zeros_like(t_s)
for i in range(len(t_s)):
    if i >= stop_index:
        continue
    T_i = T_K[i]
    k_i = A_s * math.exp(-E / (R * T_i))
    a_i = max(alpha_kin[i], 1e-4)
    dalpha_dt_kin[i] = float(f_alpha(a_i)) * k_i


# ============================================================
# 5. PRESSURE MODEL
# ============================================================

V_layer_tot = math.pi * R_layer**2 * ell_layer
V_void0   = phi_void * V_layer_tot
V_binder0 = (1.0 - phi_void) * V_layer_tot

m0 = rho_binder * V_binder0
n0 = m0 / M_gas        # max possible gas moles if α=1 (paper notation n0)

def Vp_of_alpha(a):
    """Pore volume: existing voids + gas-filled former binder volume."""
    a = float(np.clip(a, 0.0, 1.0 - 1e-12))
    return V_void0 + a * V_binder0

# ===== Gas viscosity model (Sutherland law) =====
def mu_CO2_sutherland(T_K_val):
    """CO₂ dynamic viscosity using Sutherland's correlation (typical constants)."""
    mu0 = 1.37e-5   # Pa·s at 273 K
    T0  = 273.0
    S   = 222.0
    return mu0 * (T_K_val / T0) ** 1.5 * (T0 + S) / (T_K_val + S)

def L(a):
    """
    L(α): characteristic flow length.
    (CHANGED: L_flow_alpha -> L)
    """
    a_eff = max(float(a), _ALPHA_MIN_TRANSPORT)
    return a_eff * R_layer

def tau_slit(a, T_K_val):
    """
    Characteristic diffusion/venting timescale inside a slit:
       δ(α) = α*ℓ               : slit gap
       k_frac(α) = δ(α)²/12     : parallel-plate permeability
       D_p = k_frac * p0 / μ
       τ_diff = L(α)² / D_p
    """
    a_eff = max(float(a), _ALPHA_MIN_TRANSPORT)

    delta = a_eff * ell_layer
    k_frac = (delta ** 2) / 12.0

    mu = mu_CO2_sutherland(T_K_val)

    D_p = k_frac * p0 / mu
    if D_p <= 0.0:
        return 1e12

    return (L(a_eff) ** 2) / D_p


# ===== Integrate pressure =====
n = np.zeros_like(t_s)
p_g = np.zeros_like(t_s)

a0 = max(alpha_kin[0], 1e-4)
Vp0 = Vp_of_alpha(a0)
n_eq0 = p0 * Vp0 / (R_gas * T_K[0])

n[0] = n_eq0
p_g[0] = p0

for i in range(len(t_s) - 1):
    dt = max(t_s[i + 1] - t_s[i], 1e-6)

    a_i = max(alpha_kin[i], 1e-4)
    Vp_i = Vp_of_alpha(a_i)
    n_eq_i = p0 * Vp_i / (R_gas * T_K[i])

    tau_i = tau_slit(alpha_kin[i], T_K[i])
    dn_gen_dt = n0 * dalpha_dt_kin[i]   # generation: n0 * dα/dt

    # Analytic update of: dn/dt = generation – venting
    if tau_i < 1e-9:
        n_new = n_eq_i
    else:
        expfac = math.exp(-dt / tau_i)
        n_new = n_eq_i + (n[i] - n_eq_i) * expfac + dn_gen_dt * tau_i * (1.0 - expfac)

    n[i + 1] = max(n_new, 0.0)

    a_ip1 = max(alpha_kin[i + 1], 1e-4)
    p_g[i + 1] = n[i + 1] * R_gas * T_K[i + 1] / Vp_of_alpha(a_ip1)

p_g_over_p0 = p_g / p0


# ============================================================
# 6. SAVE RESULTS
# ============================================================

out_df = pd.DataFrame({
    "time_min": t_s / 60.0,
    "T_C": T_C,
    "T_K": T_K,
    "alpha_kinetic": alpha_kin,
    "dalpha_dt_kin": dalpha_dt_kin,
    "p_g_over_p0": p_g_over_p0,
    "p_g_Pa": p_g
})
out_df.to_csv(OUTPUT_CSV, index=False)
print(f"Saved results to {OUTPUT_CSV}")


# ============================================================
# 7. OPTIONAL PLOTS
# ============================================================

plt.figure()
plt.plot(T_C, p_g_over_p0)
plt.xlabel("Temperature (°C)")
plt.ylabel("p_g/p0")
plt.grid(True)
plt.title("Pressure ratio vs Temperature")

plt.figure()
plt.plot(T_C, alpha_kin)
plt.xlabel("Temperature (°C)")
plt.ylabel("Alpha (kinetic)")
plt.grid(True)
plt.title("Kinetic conversion vs Temperature")

plt.show()
