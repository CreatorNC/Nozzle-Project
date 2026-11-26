"""Shared isentropic-flow relations and nozzle geometry (used by the other scripts)."""
import numpy as np

g0 = 9.80665

def area_ratio(M, g):
    """A/A* for a given Mach number."""
    return (1/M) * ((2/(g+1)) * (1 + (g-1)/2*M**2)) ** ((g+1)/(2*(g-1)))

def mach_from_area(AR, g, supersonic=True):
    """Invert A/A* by bisection (subsonic or supersonic branch)."""
    lo, hi = (1.0, 60.0) if supersonic else (1e-5, 1.0)
    for _ in range(200):
        mid = 0.5*(lo+hi)
        f = area_ratio(mid, g) - AR
        if (f > 0) == supersonic:
            hi = mid
        else:
            lo = mid
    return 0.5*(lo+hi)

def p_ratio(M, g):   # p/p0
    return (1 + (g-1)/2*M**2) ** (-g/(g-1))

def t_ratio(M, g):   # T/T0
    return 1/(1 + (g-1)/2*M**2)

def mach_from_pressure_ratio(pr, g):
    """Mach number for a given p/p0 (bisection)."""
    lo, hi = 1e-5, 60.0
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if p_ratio(mid, g) > pr: lo = mid
        else: hi = mid
    return 0.5*(lo+hi)

def choked_mdot(P0, T0, At, g, R):
    return P0*At*np.sqrt(g/(R*T0))*(2/(g+1))**((g+1)/(2*(g-1)))

def design(P0, T0, g, R, rt, Pa):
    """Nozzle fully expanded to Pa: returns dict of exit conditions and performance."""
    Me = mach_from_pressure_ratio(Pa/P0, g)
    AR = area_ratio(Me, g)
    At = np.pi*rt**2
    Te = T0*t_ratio(Me, g)
    ve = Me*np.sqrt(g*R*Te)
    mdot = choked_mdot(P0, T0, At, g, R)
    F = mdot*ve                      # fully expanded: pressure term is zero
    return dict(Me=Me, AR=AR, At=At, Ae=AR*At, re=rt*np.sqrt(AR), Te=Te, ve=ve,
                mdot=mdot, F=F, Isp=F/(mdot*g0))

def conical_profile(rt, re, half_conv_deg=30.0, half_div_deg=15.0, r_in_factor=1.5, n=400):
    """Axisymmetric conical nozzle: converging section, sharp throat, diverging cone.
    Returns x (throat at x=0) and radius r."""
    hc, hd = np.radians(half_conv_deg), np.radians(half_div_deg)
    r_in = r_in_factor*rt
    x_in = -(r_in - rt)/np.tan(hc)
    x_ex = (re - rt)/np.tan(hd)
    x = np.linspace(x_in, x_ex, n)
    r = np.where(x < 0, rt + (-x)*np.tan(hc), rt + x*np.tan(hd))
    return x, r
