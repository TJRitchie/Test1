import tkinter as tk
from tkinter import ttk
from fractions import Fraction
import math

# ---------- Parsing for inch inputs (supports fractions) ----------
def parse_dim(s: str) -> float:
    s = s.strip().replace("-", " ")
    parts = s.split()
    if len(parts) == 1:
        return float(Fraction(parts[0])) if "/" in parts[0] else float(parts[0])
    elif len(parts) == 2:
        whole = int(parts[0])
        frac = float(Fraction(parts[1]))
        return whole + frac if whole >= 0 else whole - frac
    else:
        raise ValueError(f"Can't parse dimension: {s}")

# ---------- Nominal calculator ----------
def calculate_nominal():
    try:
        low = parse_dim(lower_entry.get())
        high = parse_dim(upper_entry.get())
        if high < low:
            raise ValueError("Upper must be ≥ lower")
        nom = (low + high) / 2
        total_tol = high - low
        pm_tol = total_tol / 2
        result_var.set(f"Nominal: {nom:.6f}\n±Tol: {pm_tol:.6f}  (Total: {total_tol:.6f})")
    except Exception as e:
        result_var.set(f"Error: {e}")

# ---------- Unit converter (in, mm, cm, m, km, °C, °F) ----------
def try_set(var: tk.StringVar, value: str):
    if var.get() != value:
        var.set(value)

_updating = {"flag": False}  # simple re-entrancy guard across all converter fields

def _set_all_from_inches(inches: float):
    """Set all length unit fields based on an inches value."""
    _updating["flag"] = True
    try:
        mm = inches * 25.4
        cm = mm / 10.0
        m  = cm / 100.0
        km = m / 1000.0

        try_set(in_var,   f"{inches:.6f}")
        try_set(mm_var,   f"{mm:.4f}")
        try_set(cm_var,   f"{cm:.4f}")
        try_set(m_var,    f"{m:.6f}")
        try_set(km_var,   f"{km:.9f}")
        conv_err_var.set("")
    finally:
        _updating["flag"] = False

def from_in_entry(*_):
    if _updating["flag"]: return
    try:
        txt = in_var.get().strip()
        if txt == "":
            _updating["flag"] = True
            try:
                try_set(mm_var, ""); try_set(cm_var, ""); try_set(m_var, ""); try_set(km_var, "")
            finally:
                _updating["flag"] = False
            conv_err_var.set("")
            return
        inches = parse_dim(txt)
        _set_all_from_inches(inches)
    except Exception as e:
        conv_err_var.set(f"Inches error: {e}")

def from_mm_entry(*_):
    if _updating["flag"]: return
    try:
        txt = mm_var.get().strip()
        if txt == "":
            _updating["flag"] = True
            try:
                try_set(in_var, ""); try_set(cm_var, ""); try_set(m_var, ""); try_set(km_var, "")
            finally:
                _updating["flag"] = False
            conv_err_var.set("")
            return
        mm = float(txt)
        _set_all_from_inches(mm / 25.4)
    except Exception as e:
        conv_err_var.set(f"mm error: {e}")

def from_cm_entry(*_):
    if _updating["flag"]: return
    try:
        txt = cm_var.get().strip()
        if txt == "":
            _updating["flag"] = True
            try:
                try_set(in_var, ""); try_set(mm_var, ""); try_set(m_var, ""); try_set(km_var, "")
            finally:
                _updating["flag"] = False
            conv_err_var.set("")
            return
        cm = float(txt)
        _set_all_from_inches((cm * 10.0) / 25.4)
    except Exception as e:
        conv_err_var.set(f"cm error: {e}")

def from_m_entry(*_):
    if _updating["flag"]: return
    try:
        txt = m_var.get().strip()
        if txt == "":
            _updating["flag"] = True
            try:
                try_set(in_var, ""); try_set(mm_var, ""); try_set(cm_var, ""); try_set(km_var, "")
            finally:
                _updating["flag"] = False
            conv_err_var.set("")
            return
        m = float(txt)
        _set_all_from_inches(m / 0.0254)
    except Exception as e:
        conv_err_var.set(f"m error: {e}")

def from_km_entry(*_):
    if _updating["flag"]: return
    try:
        txt = km_var.get().strip()
        if txt == "":
            _updating["flag"] = True
            try:
                try_set(in_var, ""); try_set(mm_var, ""); try_set(cm_var, ""); try_set(m_var, "")
            finally:
                _updating["flag"] = False
            conv_err_var.set("")
            return
        km = float(txt)
        _set_all_from_inches(km / 0.0000254)
    except Exception as e:
        conv_err_var.set(f"km error: {e}")

# --- Temperature converters (°C ↔ °F) ---
def from_c_entry(*_):
    if _updating["flag"]: return
    try:
        txt = c_var.get().strip()
        if txt == "":
            _updating["flag"] = True
            try:
                try_set(f_var, "")
            finally:
                _updating["flag"] = False
            return
        c = float(txt)
        f = c * 9.0/5.0 + 32.0
        _updating["flag"] = True
        try:
            try_set(f_var, f"{f:.2f}")
        finally:
            _updating["flag"] = False
    except Exception as e:
        conv_err_var.set(f"Celsius error: {e}")

def from_f_entry(*_):
    if _updating["flag"]: return
    try:
        txt = f_var.get().strip()
        if txt == "":
            _updating["flag"] = True
            try:
                try_set(c_var, "")
            finally:
                _updating["flag"] = False
            return
        f = float(txt)
        c = (f - 32.0) * 5.0/9.0
        _updating["flag"] = True
        try:
            try_set(c_var, f"{c:.2f}")
        finally:
            _updating["flag"] = False
    except Exception as e:
        conv_err_var.set(f"Fahrenheit error: {e}")

# ---------- Speeds & Feeds (Lathe Turning) ----------
MATERIAL_DATA = {
    "Aluminum (6061)":      {"carbide_sfm": 800, "hss_sfm": 200, "ipr_rough": (0.010, 0.020), "ipr_finish": (0.004, 0.010)},
    "Low carbon steel (1018)": {"carbide_sfm": 500, "hss_sfm": 100, "ipr_rough": (0.010, 0.018), "ipr_finish": (0.004, 0.010)},
    "Alloy steel (4140)":   {"carbide_sfm": 350, "hss_sfm": 80,  "ipr_rough": (0.008, 0.016), "ipr_finish": (0.003, 0.008)},
    "Stainless (304/316)":  {"carbide_sfm": 250, "hss_sfm": 60,  "ipr_rough": (0.006, 0.012), "ipr_finish": (0.003, 0.006)},
    "Brass / Bronze":       {"carbide_sfm": 600, "hss_sfm": 150, "ipr_rough": (0.008, 0.016), "ipr_finish": (0.003, 0.008)},
    "Cast iron (gray)":     {"carbide_sfm": 400, "hss_sfm": 90,  "ipr_rough": (0.010, 0.020), "ipr_finish": (0.004, 0.010)},
    "Plastics (acetal/nylon)": {"carbide_sfm": 600, "hss_sfm": 200, "ipr_rough": (0.010, 0.020), "ipr_finish": (0.005, 0.015)},
}

def lathe_calc():
    try:
        mat = material_cb.get()
        tool = tool_cb.get()
        op   = op_cb.get()
        d_txt = dia_entry.get().strip()
        if not mat or not d_txt:
            sf_result_var.set("Enter diameter and choose material/tool/op.")
            return
        d_in = parse_dim(d_txt)
        if d_in <= 0:
            raise ValueError("Diameter must be > 0")

        rec = MATERIAL_DATA[mat]
        sfm = rec["carbide_sfm"] if tool == "Carbide" else rec["hss_sfm"]
        ipr_lo, ipr_hi = rec["ipr_rough"] if op == "Rough" else rec["ipr_finish"]
        ipr = (ipr_lo + ipr_hi) / 2.0

        rpm = (sfm * 12.0) / (math.pi * d_in)
        ipm = rpm * ipr

        sf_result_var.set(
            f"SFM: {sfm:.0f}\nRPM: {rpm:.0f}\nFeed: {ipr:.4f} in/rev  (≈ {ipm:.2f} in/min)"
        )
    except Exception as e:
        sf_result_var.set(f"Error: {e}")

# ---------- Fastener torque calculator ----------
PROOF_IMP_PSI = {
    "Grade 2": 55000,
    "Grade 5": 85000,
    "Grade 8": 120000,
    "Custom (psi)": None,
}
PROOF_MET_MPA = {
    "8.8": 600,
    "10.9": 830,
    "12.9": 970,
    "Custom (MPa)": None,
}
K_PRESETS = {
    "Dry steel (0.20)": 0.20,
    "Zinc plated (0.22)": 0.22,
    "Light oil (0.18)": 0.18,
    "Moly lube (0.10)": 0.10,
    "Custom": None,
}

def stress_area_imperial(d_in: float, tpi: float) -> float:
    return 0.7854 * (d_in - 0.9743 / tpi) ** 2  # in^2

def stress_area_metric(d_mm: float, pitch_mm: float) -> float:
    return (math.pi / 4.0) * (d_mm - 0.9382 * pitch_mm) ** 2  # mm^2

def calc_torque():
    try:
        system = units_var.get()
        pct = float(pct_entry.get())
        if not (0 < pct <= 100):
            raise ValueError("% clamp must be between 0 and 100")
        k = float(k_entry.get())

        if system == "imperial":
            d_in = parse_dim(diam_entry.get().strip())
            tpi = float(thread_entry.get())
            if d_in <= 0 or tpi <= 0:
                raise ValueError("Diameter and TPI must be > 0")
            At = stress_area_imperial(d_in, tpi)  # in^2

            grade = grade_cb.get()
            if grade == "Custom (psi)":
                proof = float(proof_entry.get())
            else:
                proof = PROOF_IMP_PSI[grade]
            if proof is None or proof <= 0:
                raise ValueError("Enter a valid proof (psi)")

            F = (pct / 100.0) * proof * At  # lbf
            Tinlb = k * F * d_in            # in-lbf
            TNm   = Tinlb * 0.112984829     # N·m
            Tftlb = Tinlb / 12.0

            torque_result_var.set(
                f"Stress area A_t: {At:.5f} in²\n"
                f"Preload F: {F:,.0f} lbf\n"
                f"T = K·F·d  →  {Tinlb:.1f} in·lbf   |   {Tftlb:.2f} ft·lbf   |   {TNm:.2f} N·m"
            )

        else:
            d_mm = float(diam_entry.get())
            pitch = float(thread_entry.get())
            if d_mm <= 0 or pitch <= 0:
                raise ValueError("Diameter and pitch must be > 0")
            As = stress_area_metric(d_mm, pitch)  # mm^2

            grade = grade_cb.get()
            if grade == "Custom (MPa)":
                proof_mpa = float(proof_entry.get())
            else:
                proof_mpa = PROOF_MET_MPA[grade]
            if proof_mpa is None or proof_mpa <= 0:
                raise ValueError("Enter a valid proof (MPa)")

            F_N = (pct / 100.0) * proof_mpa * As  # N
            d_m = d_mm / 1000.0
            TNm = k * F_N * d_m                   # N·m
            Tinlb = TNm * 8.850745767
            Tftlb = Tinlb / 12.0

            torque_result_var.set(
                f"Stress area A_s: {As:.1f} mm²\n"
                f"Preload F: {F_N:,.0f} N\n"
                f"T = K·F·d  →  {TNm:.2f} N·m   |   {Tftlb:.2f} ft·lbf   |   {Tinlb:.1f} in·lbf"
            )

        torque_err_var.set("")
    except Exception as e:
        torque_err_var.set(f"Error: {e}")
        torque_result_var.set("")

def on_units_change():
    system = units_var.get()
    if system == "imperial":
        diam_label.config(text="Diameter (in)")
        thread_label.config(text="Thread (TPI)")
        grade_cb.config(values=list(PROOF_IMP_PSI.keys()))
        grade_cb.set("Grade 5")
        proof_label.config(text="Custom proof (psi)")
        if grade_cb.get() != "Custom (psi)":
            proof_entry.delete(0, tk.END)
    else:
        diam_label.config(text="Diameter (mm)")
        thread_label.config(text="Pitch (mm)")
        grade_cb.config(values=list(PROOF_MET_MPA.keys()))
        grade_cb.set("10.9")
        proof_label.config(text="Custom proof (MPa)")
        if grade_cb.get() != "Custom (MPa)":
            proof_entry.delete(0, tk.END)
    torque_result_var.set("")
    torque_err_var.set("")

def on_grade_change(_=None):
    g = grade_cb.get()
    if "Custom" in g:
        proof_entry.config(state="normal")
    else:
        proof_entry.config(state="disabled")
        proof_entry.delete(0, tk.END)

def on_kpreset_change(_=None):
    name = kpreset_cb.get()
    val = K_PRESETS.get(name)
    if val is not None:
        k_entry.delete(0, tk.END)
        k_entry.insert(0, f"{val}")

# ---------- GUI setup ----------
root = tk.Tk()
root.title("Nominal Dimension Calculator")

# Nominal UI
tk.Label(root, text="Lower Limit:").grid(row=0, column=0, sticky="e", padx=5, pady=5)
lower_entry = tk.Entry(root, width=20)
lower_entry.grid(row=0, column=1, padx=5, pady=5)

tk.Label(root, text="Upper Limit:").grid(row=1, column=0, sticky="e", padx=5, pady=5)
upper_entry = tk.Entry(root, width=20)
upper_entry.grid(row=1, column=1, padx=5, pady=5)

calc_button = tk.Button(root, text="Calculate Nominal", command=calculate_nominal)
calc_button.grid(row=2, column=0, columnspan=2, pady=10)

result_var = tk.StringVar()
result_label = tk.Label(root, textvariable=result_var, font=("Arial", 12), justify="left")
result_label.grid(row=3, column=0, columnspan=2, pady=10)

# Converter UI
tk.Label(root, text="Unit Converter (type in any field)").grid(row=4, column=0, columnspan=2, pady=(10, 0))

in_var = tk.StringVar()
mm_var = tk.StringVar()
cm_var = tk.StringVar()
m_var  = tk.StringVar()
km_var = tk.StringVar()
c_var  = tk.StringVar()
f_var  = tk.StringVar()
conv_err_var = tk.StringVar()

tk.Label(root, text="Inches (in):").grid(row=5, column=0, sticky="e", padx=5, pady=3)
in_entry = tk.Entry(root, textvariable=in_var, width=20)
in_entry.grid(row=5, column=1, padx=5, pady=3)

tk.Label(root, text="Millimeters (mm):").grid(row=6, column=0, sticky="e", padx=5, pady=3)
mm_entry = tk.Entry(root, textvariable=mm_var, width=20)
mm_entry.grid(row=6, column=1, padx=5, pady=3)

tk.Label(root, text="Centimeters (cm):").grid(row=7, column=0, sticky="e", padx=5, pady=3)
cm_entry = tk.Entry(root, textvariable=cm_var, width=20)
cm_entry.grid(row=7, column=1, padx=5, pady=3)

tk.Label(root, text="Meters (m):").grid(row=8, column=0, sticky="e", padx=5, pady=3)
m_entry = tk.Entry(root, textvariable=m_var, width=20)
m_entry.grid(row=8, column=1, padx=5, pady=3)

tk.Label(root, text="Kilometers (km):").grid(row=9, column=0, sticky="e", padx=5, pady=3)
km_entry = tk.Entry(root, textvariable=km_var, width=20)
km_entry.grid(row=9, column=1, padx=5, pady=3)

# NEW: Temperature fields
tk.Label(root, text="Celsius (°C):").grid(row=10, column=0, sticky="e", padx=5, pady=3)
c_entry = tk.Entry(root, textvariable=c_var, width=20)
c_entry.grid(row=10, column=1, padx=5, pady=3)

tk.Label(root, text="Fahrenheit (°F):").grid(row=11, column=0, sticky="e", padx=5, pady=3)
f_entry = tk.Entry(root, textvariable=f_var, width=20)
f_entry.grid(row=11, column=1, padx=5, pady=3)

conv_err_label = tk.Label(root, textvariable=conv_err_var, fg="#b00020", justify="left")
conv_err_label.grid(row=12, column=0, columnspan=2, padx=5, pady=(2, 8), sticky="w")

# Live updates (bidirectional). Inch field accepts fractions.
in_var.trace_add("write", from_in_entry)
mm_var.trace_add("write", from_mm_entry)
cm_var.trace_add("write", from_cm_entry)
m_var.trace_add("write",  from_m_entry)
km_var.trace_add("write", from_km_entry)
c_var.trace_add("write",  from_c_entry)
f_var.trace_add("write",  from_f_entry)

# Enter key normalizes current field
in_entry.bind("<Return>", from_in_entry)
mm_entry.bind("<Return>", from_mm_entry)
cm_entry.bind("<Return>", from_cm_entry)
m_entry.bind("<Return>",  from_m_entry)
km_entry.bind("<Return>", from_km_entry)
c_entry.bind("<Return>",  from_c_entry)
f_entry.bind("<Return>",  from_f_entry)

# Speeds & Feeds UI (Lathe Turning)  (rows shifted by +2)
ttk.Separator(root, orient="horizontal").grid(row=13, column=0, columnspan=2, sticky="ew", padx=5, pady=(6,6))
tk.Label(root, text="Lathe Speeds & Feeds (Turning)").grid(row=14, column=0, columnspan=2, pady=(0, 4))

tk.Label(root, text="Material:").grid(row=15, column=0, sticky="e", padx=5, pady=3)
material_cb = ttk.Combobox(root, values=list(MATERIAL_DATA.keys()), state="readonly", width=26)
material_cb.grid(row=15, column=1, padx=5, pady=3)
material_cb.set("Low carbon steel (1018)")

tk.Label(root, text="Tool:").grid(row=16, column=0, sticky="e", padx=5, pady=3)
tool_cb = ttk.Combobox(root, values=["Carbide","HSS"], state="readonly", width=12)
tool_cb.grid(row=16, column=1, padx=5, pady=3, sticky="w")
tool_cb.set("Carbide")

tk.Label(root, text="Operation:").grid(row=17, column=0, sticky="e", padx=5, pady=3)
op_cb = ttk.Combobox(root, values=["Rough","Finish"], state="readonly", width=12)
op_cb.grid(row=17, column=1, padx=5, pady=3, sticky="w")
op_cb.set("Rough")

tk.Label(root, text="Diameter (in):").grid(row=18, column=0, sticky="e", padx=5, pady=3)
dia_entry = tk.Entry(root, width=20)
dia_entry.grid(row=18, column=1, padx=5, pady=3)

sf_btn = tk.Button(root, text="Calculate S&F", command=lathe_calc)
sf_btn.grid(row=19, column=0, columnspan=2, pady=6)

sf_result_var = tk.StringVar()
sf_result_lbl = tk.Label(root, textvariable=sf_result_var, justify="left", font=("Arial", 11))
sf_result_lbl.grid(row=20, column=0, columnspan=2, padx=5, pady=(0,10), sticky="w")

# Enter key runs S&F too
dia_entry.bind("<Return>", lambda e: lathe_calc())
material_cb.bind("<<ComboboxSelected>>", lambda e: lathe_calc())
tool_cb.bind("<<ComboboxSelected>>", lambda e: lathe_calc())
op_cb.bind("<<ComboboxSelected>>", lambda e: lathe_calc())

# ---------- Fastener Torque UI ----------  (rows shifted by +2)
ttk.Separator(root, orient="horizontal").grid(row=21, column=0, columnspan=2, sticky="ew", padx=5, pady=(6,6))
tk.Label(root, text="Fastener Torque (T = K · F · d)").grid(row=22, column=0, columnspan=2, pady=(0, 6))

units_var = tk.StringVar(value="imperial")
units_frame = ttk.Frame(root)
units_frame.grid(row=23, column=0, columnspan=2, sticky="w", padx=5)
ttk.Radiobutton(units_frame, text="Imperial (in, TPI)", variable=units_var, value="imperial", command=on_units_change).pack(side="left")
ttk.Radiobutton(units_frame, text="Metric (mm, pitch)", variable=units_var, value="metric", command=on_units_change).pack(side="left", padx=(10,0))

diam_label = tk.Label(root, text="Diameter (in)")
diam_label.grid(row=24, column=0, sticky="e", padx=5, pady=3)
diam_entry = tk.Entry(root, width=20)
diam_entry.grid(row=24, column=1, padx=5, pady=3)

thread_label = tk.Label(root, text="Thread (TPI)")
thread_label.grid(row=25, column=0, sticky="e", padx=5, pady=3)
thread_entry = tk.Entry(root, width=20)
thread_entry.grid(row=25, column=1, padx=5, pady=3)

tk.Label(root, text="% Clamp load:").grid(row=26, column=0, sticky="e", padx=5, pady=3)
pct_entry = tk.Entry(root, width=20)
pct_entry.insert(0, "75")
pct_entry.grid(row=26, column=1, padx=5, pady=3)

tk.Label(root, text="K factor:").grid(row=27, column=0, sticky="e", padx=5, pady=3)
k_entry = tk.Entry(root, width=20)
k_entry.insert(0, "0.20")
k_entry.grid(row=27, column=1, padx=5, pady=3)

tk.Label(root, text="K presets:").grid(row=28, column=0, sticky="e", padx=5, pady=3)
kpreset_cb = ttk.Combobox(root, values=list(K_PRESETS.keys()), state="readonly", width=20)
kpreset_cb.grid(row=28, column=1, padx=5, pady=3, sticky="w")
kpreset_cb.set("Dry steel (0.20)")
kpreset_cb.bind("<<ComboboxSelected>>", on_kpreset_change)

tk.Label(root, text="Grade / class:").grid(row=29, column=0, sticky="e", padx=5, pady=3)
grade_cb = ttk.Combobox(root, values=list(PROOF_IMP_PSI.keys()), state="readonly", width=20)
grade_cb.grid(row=29, column=1, padx=5, pady=3, sticky="w")
grade_cb.set("Grade 5")
grade_cb.bind("<<ComboboxSelected>>", on_grade_change)

proof_label = tk.Label(root, text="Custom proof (psi)")
proof_label.grid(row=30, column=0, sticky="e", padx=5, pady=3)
proof_entry = tk.Entry(root, width=20, state="disabled")
proof_entry.grid(row=30, column=1, padx=5, pady=3)

torque_btn = tk.Button(root, text="Calculate Torque", command=calc_torque)
torque_btn.grid(row=31, column=0, columnspan=2, pady=6)

torque_err_var = tk.StringVar()
tk.Label(root, textvariable=torque_err_var, fg="#b00020", justify="left").grid(row=32, column=0, columnspan=2, sticky="w", padx=5)

torque_result_var = tk.StringVar()
tk.Label(root, textvariable=torque_result_var, justify="left", font=("Arial", 11)).grid(row=33, column=0, columnspan=2, sticky="w", padx=5, pady=(0,10))

# Enter runs torque calc
for w in (diam_entry, thread_entry, pct_entry, k_entry, proof_entry):
    w.bind("<Return>", lambda e: calc_torque())

# Initialize unit-dependent UI
on_units_change()

root.mainloop()
