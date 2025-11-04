import math

def T0_over_T(M, gamma):
    return 1 + (gamma - 1) / 2 * M ** 2

def P0_over_P(M, gamma):
    return T0_over_T(M, gamma) ** (gamma / (gamma - 1))

def rho0_over_rho(M, gamma):
    return T0_over_T(M, gamma) ** (1 / (gamma - 1))

def A_over_Astar(M, gamma):
    term = (2 / (gamma + 1)) * T0_over_T(M, gamma)
    exponent = (gamma + 1) / (2 * (gamma - 1))
    return (1 / M) * term ** exponent

def solve_bisection(f, target, low, high, tol=1e-10, max_iter=200):
    a, b = low, high
    fa, fb = f(a) - target, f(b) - target
    if fa * fb > 0:
        return float('nan')
    for _ in range(max_iter):
        mid = 0.5 * (a + b)
        fm = f(mid) - target
        if abs(fm) < tol:
            return mid
        if fa * fm <= 0:
            b, fb = mid, fm
        else:
            a, fa = mid, fm
    return 0.5 * (a + b)

def mach_from_T0T(T0T, gamma):
    if T0T < 1:
        return float('nan')
    return math.sqrt((2 / (gamma - 1)) * (T0T - 1))

def mach_from_P0P(P0P, gamma):
    f = lambda M: P0_over_P(M, gamma)
    return solve_bisection(f, P0P, 1e-8, 50)

def mach_from_rho0rho(rho0rho, gamma):
    f = lambda M: rho0_over_rho(M, gamma)
    return solve_bisection(f, rho0rho, 1e-8, 50)

def mach_from_AAstar(Aa, gamma, branch="subsonic"):
    if Aa < 1:
        return float('nan')
    f = lambda M: A_over_Astar(M, gamma)
    if branch == "subsonic":
        return solve_bisection(f, Aa, 1e-8, 1 - 1e-8)
    else:
        return solve_bisection(f, Aa, 1 + 1e-8, 50)

def isentropic_calculator(input_type, input_value, gamma=1.4, branch="subsonic"):
    if input_type == "M":
        M = input_value
    elif input_type == "T0T":
        M = mach_from_T0T(input_value, gamma)
    elif input_type == "P0P":
        M = mach_from_P0P(input_value, gamma)
    elif input_type == "rho0rho":
        M = mach_from_rho0rho(input_value, gamma)
    elif input_type == "AoverAstar":
        M = mach_from_AAstar(input_value, gamma, branch)
    else:
        return None

    if not math.isfinite(M):
        return None

    results = {
        'Mach': M,
        'T0/T': T0_over_T(M, gamma),
        'P0/P': P0_over_P(M, gamma),
        'rho0/rho': rho0_over_rho(M, gamma),
        'A/A*': A_over_Astar(M, gamma)
    }
    return results

if __name__ == "__main__":
    print("Isentropic Flow Calculator (γ adjustable)")
    gamma = float(input("Enter gamma (default 1.4): ") or 1.4)
    print("Choose input type: M, T0T, P0P, rho0rho, AoverAstar")
    t = input("Type: ")
    val = float(input("Value: "))
    branch = "subsonic"
    if t == "AoverAstar":
        branch = input("Branch (subsonic/supersonic): ") or "subsonic"
    res = isentropic_calculator(t, val, gamma, branch)
    if res:
        for k, v in res.items():
            print(f"{k}: {v:.6f}")
    else:
        print("Invalid input.")
















