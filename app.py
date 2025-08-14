from flask import Flask, render_template, request
import json, math, os

app = Flask(__name__)

# --- Load shapes database (if present) ---
SHAPES_PATH = os.path.join(os.path.dirname(__file__), 'shapes.json')
try:
    with open(SHAPES_PATH, 'r', encoding='utf-8') as f:
        SHAPES = json.load(f)
except Exception:
    SHAPES = {}  # empty; you can create shapes.json using the sample below

# ---------------- Beam helpers (UDL on simple span) ----------------
# Units: w [kN/m], L [m], Fy [MPa], Zx [mm^3], phi_b [-]
# Returns kN·m and kN, plus % margin

def beam_calc(w_kN_per_m: float, L_m: float, Fy_MPa: float, Zx_mm3: float, phi_b: float = 0.9):
    # Max bending moment and shear for UDL on simply supported beam
    M_u_kNm = w_kN_per_m * (L_m ** 2) / 8.0  # kN·m
    V_max_kN = w_kN_per_m * L_m / 2.0        # kN

    # Nominal moment resistance (very simplified): MR = φ * Fy * Zx
    # Convert Fy[MPa]=N/mm^2 and Zx[mm^3] to N·mm, then to kN·m
    MR_kNm = phi_b * Fy_MPa * Zx_mm3 / 1e6   # (N/mm^2 * mm^3 = N·mm) /1e6 -> kN·m

    ratio_percent = (MR_kNm / M_u_kNm) * 100.0 if M_u_kNm > 0 else float('inf')
    return {
        'Mu_kNm': round(M_u_kNm, 3),
        'Vmax_kN': round(V_max_kN, 3),
        'MR_kNm': round(MR_kNm, 3),
        'MR_over_Mu_percent': round(ratio_percent, 1)
    }

# ---------------- Load combination helper (Pu) ----------------
# Inputs: D[kN], L[kN]
# Combinations considered: 1.4D ; 1.25D + 1.5L ; 0.9L (per user spec)

def pu_combo(D_kN: float, L_kN: float):
    combos = {
        '1.4D': 1.4 * D_kN,
        '1.25D+1.5L': 1.25 * D_kN + 1.5 * L_kN,
        '0.9L': 0.9 * L_kN,
    }
    gov_name, Pu = max(combos.items(), key=lambda kv: kv[1])
    return {
        'Pu_gov': round(Pu, 2),
        'governing': gov_name,
        'combos': {k: round(v, 2) for k, v in combos.items()}
    }

# ---------------- Column helpers (Class 3 sections) ----------------
# Inputs: Pu[kN], L[m], Kx, Ky, A[mm^2], rx[mm], ry[mm], E[MPa], Fy[MPa], phi_c, n
# Outputs: capacities for Yielding, Euler (x & y), and Column Curve (real), plus governing mode

def col_calc(Pu_kN: float, L_m: float, Kx: float, Ky: float,
             A_mm2: float, rx_mm: float, ry_mm: float,
             E_MPa: float, Fy_MPa: float, phi_c: float = 0.9, n_exp: float = 1.34):
    # Conversions
    L_mm = L_m * 1000.0

    # 1) Yielding (ideal): P_y = φ * A * Fy  (N) -> kN
    P_y_kN = phi_c * A_mm2 * Fy_MPa / 1000.0

    # 2) Global (flexural) buckling (Euler) about x and y using I = A r^2
    #    P_cr = φ * π^2 * E * I / (K L)^2 = φ * π^2 * E * A * r^2 / (K L)^2
    def pcr(axis_r_mm: float, K: float):
        if axis_r_mm <= 0 or K <= 0 or L_mm <= 0:
            return float('nan')
        return phi_c * (math.pi ** 2) * E_MPa * A_mm2 * (axis_r_mm ** 2) / ((K * L_mm) ** 2) / 1000.0

    P_crx_kN = pcr(rx_mm, Kx)
    P_cry_kN = pcr(ry_mm, Ky)

    # Corresponding elastic buckling stresses Fe_x, Fe_y (MPa)
    # Fe = π^2 E / (KL/r)^2 = π^2 E * r^2 / (K L)^2
    def fe(axis_r_mm: float, K: float):
        if axis_r_mm <= 0 or K <= 0 or L_mm <= 0:
            return float('nan')
        return (math.pi ** 2) * E_MPa * (axis_r_mm ** 2) / ((K * L_mm) ** 2)

    Fe_x = fe(rx_mm, Kx)
    Fe_y = fe(ry_mm, Ky)
    Fe_min = min(Fe_x, Fe_y)

    # 3) Column curve (real column) using generalized Rankine with exponent n:
    #     1 / F_r^n = 1 / Fy^n + 1 / Fe_min^n   ->  F_r = ( Fy^{-n} + Fe^{-n} ) ^ (-1/n)
    # Note: valid for Class 3 sections. For singly/asymmetric sections, use Fe_min from appropriate buckling check.
    if Fe_min and Fy_MPa and Fe_min > 0 and Fy_MPa > 0:
        Fr_MPa = ( (Fy_MPa ** (-n_exp)) + (Fe_min ** (-n_exp)) ) ** (-1.0 / n_exp)
    else:
        Fr_MPa = float('nan')

    P_curve_kN = phi_c * A_mm2 * Fr_MPa / 1000.0

    # Slenderness ratios
    lam_x = (Kx * L_mm / rx_mm) if rx_mm > 0 else float('inf')
    lam_y = (Ky * L_mm / ry_mm) if ry_mm > 0 else float('inf')

    # Pick governing nominal resistance
    candidates = [
        ('Yielding', P_y_kN),
        ('Euler-x', P_crx_kN),
        ('Euler-y', P_cry_kN),
        ('Column Curve', P_curve_kN)
    ]
    # Filter out NaNs
    candidates = [(name, val) for name, val in candidates if not (val is None or isinstance(val, float) and math.isnan(val))]

    governing = min(candidates, key=lambda t: t[1]) if candidates else ('N/A', float('nan'))

    # Utilization/margins
    def util(Pn):
        return (Pu_kN / Pn * 100.0) if Pn and Pn > 0 else float('inf')

    result = {
        'P_y_kN': round(P_y_kN, 2),
        'P_crx_kN': round(P_crx_kN, 2),
        'P_cry_kN': round(P_cry_kN, 2),
        'P_curve_kN': round(P_curve_kN, 2),
        'util_yield_percent': round(util(P_y_kN), 1),
        'util_eulerx_percent': round(util(P_crx_kN), 1),
        'util_eulery_percent': round(util(P_cry_kN), 1),
        'util_curve_percent': round(util(P_curve_kN), 1),
        'Fe_x_MPa': round(Fe_x, 1),
        'Fe_y_MPa': round(Fe_y, 1),
        'Fe_min_MPa': round(Fe_min, 1),
        'Fr_MPa': round(Fr_MPa, 1),
        'lambda_x': round(lam_x, 1),
        'lambda_y': round(lam_y, 1),
        'governing_mode': governing[0],
        'Pn_governing_kN': round(governing[1], 2),
        'Pu_over_Pn_percent': round(util(governing[1]), 1)
    }
    return result

@app.route('/', methods=['GET', 'POST'])
def index():
    beam_out = None
    col_out = None
    pu_out = None

    if request.method == 'POST':
        form_type = request.form.get('form_type')

        if form_type == 'beam':
            w = float(request.form.get('w', 0))
            L = float(request.form.get('L', 0))
            Fy = float(request.form.get('Fy', 350))
            Zx = float(request.form.get('Zx', 0))
            phi_b = float(request.form.get('phi_b', 0.9))
            beam_out = beam_calc(w, L, Fy, Zx, phi_b)

        elif form_type == 'pu_combo':
            D = float(request.form.get('D', 0))
            Lvl = float(request.form.get('L', 0))
            pu_out = pu_combo(D, Lvl)

        elif form_type == 'column':
            Pu = float(request.form.get('Pu', 0))
            L = float(request.form.get('L', 0))
            Kx = float(request.form.get('Kx', 1.0))
            Ky = float(request.form.get('Ky', 1.0))
            A = float(request.form.get('A', 0))
            rx = float(request.form.get('rx', 0))
            ry = float(request.form.get('ry', 0))
            E = float(request.form.get('E', 200000))
            Fy = float(request.form.get('Fy', 350))
            phi_c = float(request.form.get('phi_c', 0.9))
            n_exp = float(request.form.get('n_exp', 1.34))
            col_out = col_calc(Pu, L, Kx, Ky, A, rx, ry, E, Fy, phi_c, n_exp)

    # pass shapes and a friendly list of names for the selector
    shape_names = sorted(list(SHAPES.keys()))
    return render_template(
        'index.html',
        beam_out=beam_out,
        col_out=col_out,
        pu_out=pu_out,
        shapes=SHAPES,
        shape_names=shape_names,
    )

if __name__ == '__main__':
    app.run(debug=True)