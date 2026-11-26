"""
Quasi-1D Euler solver for a converging-diverging nozzle (finite volume, HLL flux, RK2 time marching).
Compares the converged numerical solution with the exact isentropic solution and runs a mesh refinement study.
Run:  python nozzle_cfd1d.py
"""
import os
print(os.getcwd())
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from nozzle_lib import *

# ---------------- inputs ----------------
g, R = 1.2, 350.0          # exhaust gas properties
P0, T0 = 5.0e6, 3000.0     # chamber conditions
rt, Pa = 0.05, 101325.0    # throat radius, design ambient pressure
cp = g*R/(g-1)

d = design(P0, T0, g, R, rt, Pa)
x_geo, r_geo = conical_profile(rt, d['re'])
x_in, x_ex = x_geo[0], x_geo[-1]
L = x_ex - x_in
At = np.pi*rt**2

def radius(x):
    hc, hd = np.radians(30), np.radians(15)
    return np.where(x < 0, rt + (-x)*np.tan(hc), rt + x*np.tan(hd))

def exact_mach(x):
    AR = (radius(x)/rt)**2
    return np.array([mach_from_area(a, g, supersonic=(xi > 0)) if a > 1 + 1e-12 else 1.0
                     for a, xi in zip(AR, x)])

def hll(rL, uL, pL, rR, uR, pR):
    aL, aR = np.sqrt(g*pL/rL), np.sqrt(g*pR/rR)
    EL = pL/((g-1)*rL) + 0.5*uL**2
    ER = pR/((g-1)*rR) + 0.5*uR**2
    UL = np.array([rL, rL*uL, rL*EL]); UR = np.array([rR, rR*uR, rR*ER])
    FL = np.array([rL*uL, rL*uL**2 + pL, (rL*EL + pL)*uL])
    FR = np.array([rR*uR, rR*uR**2 + pR, (rR*ER + pR)*uR])
    SL = np.minimum(uL - aL, uR - aR); SR = np.maximum(uL + aL, uR + aR)
    mid = (SR*FL - SL*FR + SL*SR*(UR - UL))/(SR - SL)
    return np.where(SL >= 0, FL, np.where(SR <= 0, FR, mid))

def solve(N, n_flow_through=30, cfl=0.4):
    xf = np.linspace(x_in, x_ex, N+1); xc = 0.5*(xf[1:] + xf[:-1]); dx = xf[1] - xf[0]
    Af, Ac = np.pi*radius(xf)**2, np.pi*radius(xc)**2
    # start-up: cold gas at rest at ambient pressure, chamber gas pushed in through the inlet
    rho = np.full(N, Pa/(R*300.0)); u = np.zeros(N); p = np.full(N, Pa)
    q = np.array([rho*Ac, rho*u*Ac, (p/(g-1) + 0.5*rho*u**2)*Ac])

    def rhs(q):
        rho = q[0]/Ac; u = q[1]/q[0]; p = (g-1)*(q[2]/Ac - 0.5*rho*u**2)
        p = np.maximum(p, 1.0); rho = np.maximum(rho, 1e-6)
        # subsonic inlet from stagnation conditions, supersonic outlet by extrapolation
        ui = max(u[0], 0.0); Ti = T0 - ui**2/(2*cp)
        pi = P0*(Ti/T0)**(g/(g-1)); ri = pi/(R*Ti)
        rl = np.concatenate(([ri], rho)); ul = np.concatenate(([ui], u)); pl = np.concatenate(([pi], p))
        rr = np.concatenate((rho, [rho[-1]])); ur = np.concatenate((u, [u[-1]])); pr = np.concatenate((p, [p[-1]]))
        F = hll(rl, ul, pl, rr, ur, pr)*Af
        dq = -(F[:, 1:] - F[:, :-1])/dx
        dq[1] += p*(Af[1:] - Af[:-1])/dx
        a = np.sqrt(g*p/rho)
        return dq, np.max(np.abs(u) + a)

    t, t_end = 0.0, n_flow_through*L/d['ve']
    while t < t_end:
        k1, smax = rhs(q); dt = cfl*dx/smax
        k2, _ = rhs(q + dt*k1)
        q = q + 0.5*dt*(k1 + k2); t += dt
    rho = q[0]/Ac; u = q[1]/q[0]; p = (g-1)*(q[2]/Ac - 0.5*rho*u**2)
    M = u/np.sqrt(g*p/rho)
    return xc, M, p, rho, u, Ac

# ---------------- mesh refinement study ----------------
results, sols = [], {}
mdot_exact = choked_mdot(P0, T0, At, g, R)
for N in (50, 100, 200, 400):
    xc, M, p, rho, u, Ac = solve(N)
    Mex = exact_mach(xc)
    e2 = np.sqrt(np.mean(((M - Mex)/Mex)**2))*100           # relative L2 error in Mach [%]
    mdot = rho*u*Ac
    md_err = (np.mean(mdot[1:]) - mdot_exact)/mdot_exact*100
    Me_err = (M[-1] - Mex[-1])/Mex[-1]*100
    F_num = rho[-1]*u[-1]**2*Ac[-1] + (p[-1] - Pa)*Ac[-1]
    F_err = (F_num - d['F'])/d['F']*100
    results.append(dict(N=N, L2_Mach_pct=e2, exit_Mach_err_pct=Me_err, mdot_err_pct=md_err, thrust_err_pct=F_err))
    sols[N] = (xc, M, p)
    print(f"N={N:4d}  Mach L2 error {e2:6.3f}%  exit-Mach error {Me_err:+6.3f}%  mass-flow error {md_err:+6.3f}%  thrust error {F_err:+6.3f}%")

orders = [float(np.log2(results[i]["L2_Mach_pct"]/results[i+1]["L2_Mach_pct"])) for i in range(len(results)-1)]
print("Observed order of convergence between successive meshes:", [round(o, 2) for o in orders])
json.dump(dict(design={k: float(v) for k, v in d.items()}, mesh_study=results, orders=orders),
          open("results.json", "w"), indent=2)

# ---------------- plots ----------------
xc, M, p = sols[400]; xe = np.linspace(x_in, x_ex, 800); Mex = exact_mach(xe)
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(xe*100, Mex, 'k-', label="exact isentropic"); ax[0].plot(xc*100, M, 'r--', label="finite volume (N=400)")
ax[0].set_xlabel("x [cm]"); ax[0].set_ylabel("Mach number"); ax[0].legend(); ax[0].grid(alpha=.3); ax[0].set_title("Mach number along the nozzle")
pex = P0*p_ratio(Mex, g)/1e6
ax[1].plot(xe*100, pex, 'k-', label="exact isentropic"); ax[1].plot(xc*100, p/1e6, 'r--', label="finite volume (N=400)")
ax[1].set_xlabel("x [cm]"); ax[1].set_ylabel("Pressure [MPa]"); ax[1].legend(); ax[1].grid(alpha=.3); ax[1].set_title("Static pressure along the nozzle")
plt.tight_layout(); plt.savefig("cfd1d_vs_exact.png", dpi=150); plt.close()

Ns = [r['N'] for r in results]; errs = [r['L2_Mach_pct'] for r in results]
plt.figure(figsize=(5.5, 4)); plt.loglog(Ns, errs, 'o-', label="Mach L2 error")
plt.loglog(Ns, [errs[0]*Ns[0]/n for n in Ns], 'k--', label="first-order slope")
plt.xlabel("Number of cells"); plt.ylabel("Relative L2 error [%]"); plt.grid(alpha=.3, which="both"); plt.legend()
plt.title("Mesh refinement study"); plt.tight_layout(); plt.savefig("mesh_convergence.png", dpi=150); plt.close()

# ---------------- geometry export for CAD / meshing ----------------
with open("nozzle_profile.csv", "w") as f:
    f.write("x_mm,r_mm\n")
    for xi, ri in zip(x_geo[::8], r_geo[::8]):
        f.write(f"{xi*1000:.3f},{ri*1000:.3f}\n")
print("\nDesign: Me=%.3f  Ae/At=%.3f  exit radius=%.1f mm  convergent length=%.1f mm  divergent length=%.1f mm"
      % (d['Me'], d['AR'], d['re']*1000, -x_in*1000, x_ex*1000))
