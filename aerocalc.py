import math
import numpy as np
import streamlit as st


# === Core Functions ===
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
    if input_type == "Mach":
        M = input_value
    elif input_type == "T0/T":
        M = mach_from_T0T(input_value, gamma)
    elif input_type == "P0/P":
        M = mach_from_P0P(input_value, gamma)
    elif input_type == "ρ0/ρ":
        M = mach_from_rho0rho(input_value, gamma)
    elif input_type == "A/A*":
        M = mach_from_AAstar(input_value, gamma, branch)
    else:
        return None

    if not math.isfinite(M):
        return None

    return {
        'Mach': M,
        'T0/T': T0_over_T(M, gamma),
        'P0/P': P0_over_P(M, gamma),
        'ρ0/ρ': rho0_over_rho(M, gamma),
        'A/A*': A_over_Astar(M, gamma)
    }

# === Streamlit UI ===
st.set_page_config(page_title="Isentropic Flow Calculator", page_icon="🌬️", layout="centered")

# Custom styling
st.markdown("""
    <style>
        :root {
            --osu-orange: #FF6600;
        }
        .stApp {
            background-color: #121212;
            color: white;
        }
        h1, h2, h3 {
            color: var(--osu-orange);
        }
        .stButton button {
            background-color: var(--osu-orange);
            color: white;
            border: none;
            border-radius: 0.5rem;
            padding: 0.5rem 1rem;
            font-weight: bold;
        }
        .stButton button:hover {
            background-color: #e65c00;
        }
    </style>
""", unsafe_allow_html=True)

st.title("🌬️ Isentropic Flow Calculator")
st.caption("A compressible flow tool built with Python — adjustable γ, subsonic/supersonic modes, and OSU styling.")

gamma = st.number_input("Ratio of Specific Heats (γ)", value=1.4, min_value=1.0, step=0.01)
input_type = st.selectbox("Select Input Type", ["Mach", "T0/T", "P0/P", "ρ0/ρ", "A/A*"])
input_value = st.number_input(f"Enter value for {input_type}", value=2.0, min_value=0.0)

branch = "subsonic"
if input_type == "A/A*":
    branch = st.radio("Flow Regime", ["subsonic", "supersonic"], horizontal=True)

if st.button("Calculate"):
    result = isentropic_calculator(input_type, input_value, gamma, branch)
    if result:
        st.success("Calculation Successful ✅")
        st.write("### Results:")
        for key, value in result.items():
            st.write(f"**{key}:** {value:.6f}")

        # Plot section
        st.write("### 📈 Property Variation with Mach Number")
        M_vals = np.linspace(0.1, 5, 200)
        T_vals = [T0_over_T(M, gamma) for M in M_vals]
        P_vals = [P0_over_P(M, gamma) for M in M_vals]
        rho_vals = [rho0_over_rho(M, gamma) for M in M_vals]
        A_vals = [A_over_Astar(M, gamma) for M in M_vals]

        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(M_vals, T_vals, label="T₀/T")
        ax.plot(M_vals, P_vals, label="P₀/P")
        ax.plot(M_vals, rho_vals, label="ρ₀/ρ")
        ax.plot(M_vals, A_vals, label="A/A*")
        ax.axvline(result['Mach'], color="#FF6600", linestyle="--", label="Current Mach")
        ax.set_xlabel("Mach Number (M)")
        ax.set_ylabel("Ratio Value")
        ax.set_title("Isentropic Flow Ratios vs Mach Number")
        ax.legend()
        ax.grid(True, linestyle="--", alpha=0.4)
        st.pyplot(fig)
    else:
        st.error("Invalid input or no solution found.")

