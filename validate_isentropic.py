import os
print(os.getcwd())
"""Validate the isentropic relations against standard gas-dynamics table values (gamma = 1.4)."""
from nozzle_lib import area_ratio, p_ratio, t_ratio, mach_from_area

g = 1.4
# Standard isentropic-flow table values: M, A/A*, p/p0, T/T0
TABLE = [(1.5, 1.1762, 0.2724, 0.6897),
         (2.0, 1.6875, 0.1278, 0.5556),
         (3.0, 4.2346, 0.02722, 0.3571)]

print(f"{'M':>4} | {'A/A* calc':>10} {'table':>8} {'err %':>7} | {'p/p0 calc':>10} {'table':>8} {'err %':>7} | {'T/T0 calc':>10} {'table':>8} {'err %':>7}")
worst = 0.0
for M, A, P, T in TABLE:
    a, p, t = area_ratio(M, g), p_ratio(M, g), t_ratio(M, g)
    ea, ep, et = [abs(c-r)/r*100 for c, r in ((a, A), (p, P), (t, T))]
    worst = max(worst, ea, ep, et)
    print(f"{M:4.1f} | {a:10.4f} {A:8.4f} {ea:7.3f} | {p:10.5f} {P:8.5f} {ep:7.3f} | {t:10.4f} {T:8.4f} {et:7.3f}")

# Round trip: Mach -> A/A* -> Mach
rt_err = max(abs(mach_from_area(area_ratio(M, g), g) - M) for M, *_ in TABLE)
print(f"\nWorst error against table: {worst:.3f} %  (table values are rounded to 4 significant digits)")
print(f"Round-trip Mach error (M -> A/A* -> M): {rt_err:.2e}")
