"""
Rocket Nozzle Design & Analysis (1D isentropic model)
-----------------------------------------------------
Inputs: chamber conditions, propellant gas properties, throat size, ambient pressure.
Outputs: exit Mach, area ratio, thrust, Isp, Mach/pressure/temperature along the nozzle,
         and a parametric study of how expansion ratio and altitude affect performance.
Requires: numpy, matplotlib
"""
import os
print(os.getcwd())
import numpy as np
import matplotlib.pyplot as plt

g0 = 9.80665

# ---------- Inputs (edit these) ----------
P0 = 5.0e6        # chamber pressure [Pa]
T0 = 3000.0       # chamber temperature [K]
gamma = 1.2       # specific heat ratio of exhaust gas
R = 350.0         # gas constant of exhaust [J/kg/K]
rt = 0.05         # throat radius [m]
Pa = 101325.0     # ambient pressure [Pa] (sea level)
half_div = np.radians(15)   # divergent cone half-angle
half_conv = np.radians(30)  # convergent half-angle
# -----------------------------------------

At = np.pi * rt**2

def area_ratio(M, g=gamma):
    """A/A* as a function of Mach number."""
    return (1/M) * ((2/(g+1)) * (1 + (g-1)/2 * M**2)) ** ((g+1)/(2*(g-1)))

def mach_from_area(AR, supersonic=True, g=gamma):
    """Invert A/A* by bisection."""
    lo, hi = (1.0, 50.0) if supersonic else (1e-4, 1.0)
    for _ in range(100):
        mid = 0.5*(lo+hi)
        f = area_ratio(mid, g) - AR
        if (f > 0) == supersonic:
            hi = mid
        else:
            lo = mid
    return 0.5*(lo+hi)

def p_ratio(M, g=gamma):   # p/p0
    return (1 + (g-1)/2*M**2) ** (-g/(g-1))

def t_ratio(M, g=gamma):   # T/T0
    return 1 / (1 + (g-1)/2*M**2)

def performance(P0, T0, Pe, Pa, At, g=gamma, R=R):
    Te = T0 * (Pe/P0) ** ((g-1)/g)
    ve = np.sqrt(2*g/(g-1) * R * T0 * (1 - (Pe/P0)**((g-1)/g)))
    Me = ve / np.sqrt(g*R*Te)
    AR = area_ratio(Me, g)
    Ae = AR * At
    mdot = P0 * At * np.sqrt(g/(R*T0)) * (2/(g+1)) ** ((g+1)/(2*(g-1)))
    F = mdot*ve + (Pe - Pa)*Ae
    return dict(Me=Me, AR=AR, Ae=Ae, ve=ve, mdot=mdot, F=F, Isp=F/(mdot*g0), Te=Te)

# ---------- 1. Design point: nozzle fully expanded to ambient pressure ----------
d = performance(P0, T0, Pa, Pa, At)
print("DESIGN POINT (fully expanded at sea level)")
print(f"  Exit Mach        : {d['Me']:.2f}")
print(f"  Area ratio Ae/At : {d['AR']:.2f}")
print(f"  Exit velocity    : {d['ve']:.0f} m/s")
print(f"  Mass flow        : {d['mdot']:.2f} kg/s")
print(f"  Thrust           : {d['F']/1e3:.1f} kN")
print(f"  Isp              : {d['Isp']:.0f} s")

# ---------- 2. Flow properties along a conical nozzle ----------
re = np.sqrt(d['Ae']/np.pi)
x_conv = (1.5*rt - rt) / np.tan(half_conv)         # converging section length
x_div = (re - rt) / np.tan(half_div)               # diverging section length
x = np.linspace(-x_conv, x_div, 400)
r = np.where(x < 0, rt + (-x)*np.tan(half_conv), rt + x*np.tan(half_div))
AR_x = (r/rt)**2
M = np.array([mach_from_area(a, supersonic=(xi > 0)) if a > 1.0001 else 1.0
              for a, xi in zip(AR_x, x)])
P = P0 * p_ratio(M)
T = T0 * t_ratio(M)

fig, ax = plt.subplots(2, 2, figsize=(11, 8))
ax[0,0].plot(x*100, r*100, 'k'); ax[0,0].plot(x*100, -r*100, 'k')
ax[0,0].set_title("Nozzle contour"); ax[0,0].set_xlabel("x [cm]"); ax[0,0].set_ylabel("radius [cm]"); ax[0,0].axis('equal')
ax[0,1].plot(x*100, M); ax[0,1].set_title("Mach number"); ax[0,1].set_xlabel("x [cm]")
ax[1,0].plot(x*100, P/1e6); ax[1,0].set_title("Static pressure [MPa]"); ax[1,0].set_xlabel("x [cm]")
ax[1,1].plot(x*100, T); ax[1,1].set_title("Static temperature [K]"); ax[1,1].set_xlabel("x [cm]")
for a in ax.flat: a.grid(alpha=.3)
plt.tight_layout(); plt.savefig("nozzle_flow_profiles.png", dpi=150); plt.close()

# ---------- 3. Parametric study: fixed nozzle (design area ratio) vs altitude ----------
alt = np.linspace(0, 40e3, 80)
Pa_alt = 101325 * np.exp(-alt/7200)                # simple exponential atmosphere
Me_fixed = d['Me']; Pe_fixed = P0 * p_ratio(Me_fixed)
F_alt = []
for pa in Pa_alt:
    F_alt.append(d['mdot']*d['ve'] + (Pe_fixed - pa)*d['Ae'])
F_alt = np.array(F_alt)
plt.figure(figsize=(6,4))
plt.plot(alt/1e3, F_alt/1e3); plt.xlabel("Altitude [km]"); plt.ylabel("Thrust [kN]")
plt.title("Thrust vs altitude (fixed nozzle)"); plt.grid(alpha=.3)
plt.tight_layout(); plt.savefig("thrust_vs_altitude.png", dpi=150); plt.close()

# ---------- 4. Parametric study: Isp vs expansion ratio in vacuum ----------
ARs = np.linspace(3, 100, 60)
Isp_vac = []
for AR in ARs:
    Me = mach_from_area(AR)
    Pe = P0 * p_ratio(Me)
    p = performance(P0, T0, Pe, 0.0, At)
    Isp_vac.append(p['Isp'])
plt.figure(figsize=(6,4))
plt.plot(ARs, Isp_vac); plt.xlabel("Expansion ratio Ae/At"); plt.ylabel("Vacuum Isp [s]")
plt.title("Vacuum Isp vs expansion ratio"); plt.grid(alpha=.3)
plt.tight_layout(); plt.savefig("isp_vs_expansion_ratio.png", dpi=150); plt.close()
print("\nSaved: nozzle_flow_profiles.png, thrust_vs_altitude.png, isp_vs_expansion_ratio.png")
