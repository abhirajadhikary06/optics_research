# double_slit_analysis.py
# IEEE Project: Double-Slit Diffraction & Interference
# Enhanced with 4 graphs + 3 tables + uncertainty + visibility + envelope analysis

import numpy as np
import matplotlib.pyplot as plt
import cv2
import pandas as pd
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import os

# ===============================
# 1. SIMULATION OF INTERFERENCE
# ===============================
def simulate_double_slit(wavelength=632.8e-9, d_list=[0.125e-3, 0.250e-3], L=1.50,
                        screen_width=0.05, pixels=2048, slit_width=0.04e-3):
    results = {}
    x = np.linspace(-screen_width/2, screen_width/2, pixels)
    
    for d in d_list:
        # Single-slit envelope
        beta = (np.pi * slit_width * x) / (wavelength * L)
        envelope = (np.sinc(beta / np.pi))**2
        
        # Double-slit interference
        delta = (np.pi * d * x) / (wavelength * L)
        interference = np.cos(delta)**2
        
        I = envelope * interference
        I /= I.max()
        results[d] = (x * 1e3, I, envelope)  # x in mm, return envelope too
    return results

# ===============================
# 2. IMAGE ANALYSIS
# ===============================
def analyze_image(image_path, crop_roi=None, pixel_to_mm=0.025):
    img = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
    if img is None:
        raise FileNotFoundError(f"Image not found: {image_path}")
    if crop_roi:
        x, y, w, h = crop_roi
        img = img[y:y+h, x:x+w]
    
    profile = np.mean(img, axis=0).astype(float)
    profile /= profile.max()
    
    x_pixels = np.arange(len(profile))
    x_mm = x_pixels * pixel_to_mm
    
    return x_mm, profile, img

# ===============================
# 3. FRINGE SPACING & VISIBILITY
# ===============================
def find_fringe_spacing(x, intensity, min_prominence=0.1):
    peaks, props = find_peaks(intensity, prominence=min_prominence)
    if len(peaks) < 3:
        return None, None, None, None
    diffs = np.diff(x[peaks])
    spacing = np.mean(diffs)
    std = np.std(diffs)
    # Visibility = (Imax - Imin) / (Imax + Imin) ~ average over central region
    central = slice(peaks[0], peaks[-1])
    I_max = np.max(intensity[central])
    I_min = np.min(intensity[central])
    visibility = (I_max - I_min) / (I_max + I_min) if (I_max + I_min) > 0 else 0
    return spacing, std, visibility, peaks

# ===============================
# 4. ENVELOPE FIT (Single-Slit Diffraction)
# ===============================
def sinc_envelope(x, A, w):
    beta = (np.pi * w * 1e-3 * x) / (632.8e-9 * 1.50)
    return A * (np.sinc(beta / np.pi))**2

def fit_envelope(x, I, p0=[1.0, 0.04]):
    try:
        popt, _ = curve_fit(sinc_envelope, x, I, p0=p0, bounds=(0, [2, 0.1]))
        return popt[1]  # fitted slit width in mm
    except:
        return None

# ===============================
# 5. COMPARATIVE TABLE 1: Fringe Spacing
# ===============================
def make_fringe_table(sim_data, real_data=None, wavelength=632.8, L=1.50):
    rows = []
    for d_mm in sim_data.keys():
        d = d_mm * 1e-3
        theo = (wavelength * 1e9 * L) / d
        
        x_sim, I_sim, _ = sim_data[d_mm]
        sim_spacing, sim_std, sim_vis, _ = find_fringe_spacing(x_sim, I_sim)
        
        real_spacing = real_std = real_vis = None
        if real_data and d_mm in real_data:
            x_real, I_real, _ = real_data[d_mm]
            real_spacing, real_std, real_vis, _ = find_fringe_spacing(x_real, I_real)
        
        rows.append({
            'Slit Sep. $d$ (mm)': d_mm,
            'Theoretical $\\Delta y$ (mm)': f"{theo:.2f}",
            'Simulated $\\Delta y$ (mm)': f"{sim_spacing:.2f} $\\pm$ {sim_std:.2f}" if sim_spacing else "—",
            'Measured $\\Delta y$ (mm)': f"{real_spacing:.2f} $\\pm$ {real_std:.2f}" if real_spacing else "—",
            '% Error (Sim)': f"{abs(sim_spacing - theo)/theo*100:.1f}\%" if sim_spacing else "—",
            '% Error (Exp)': f"{abs(real_spacing - theo)/theo*100:.1f}\%" if real_spacing else "—"
        })
    return pd.DataFrame(rows)

# ===============================
# 6. TABLE 2: Fringe Visibility
# ===============================
def make_visibility_table(sim_data, real_data=None):
    rows = []
    for d_mm in sim_data.keys():
        x_sim, I_sim, _ = sim_data[d_mm]
        _, _, sim_vis, _ = find_fringe_spacing(x_sim, I_sim)
        
        real_vis = None
        if real_data and d_mm in real_data:
            x_real, I_real, _ = real_data[d_mm]
            _, _, real_vis, _ = find_fringe_spacing(x_real, I_real)
        
        rows.append({
            'Slit Sep. $d$ (mm)': d_mm,
            'Simulated Visibility': f"{sim_vis:.3f}" if sim_vis is not None else "—",
            'Measured Visibility': f"{real_vis:.3f}" if real_vis is not None else "—"
        })
    return pd.DataFrame(rows)

# ===============================
# 7. TABLE 3: Envelope Fit
# ===============================
def make_envelope_table(sim_data, real_data=None):
    rows = []
    true_width = 0.040  # mm
    for d_mm in sim_data.keys():
        x_sim, I_sim, envelope = sim_data[d_mm]
        fitted_width = fit_envelope(x_sim, envelope)
        
        real_width = None
        if real_data and d_mm in real_data:
            x_real, I_real, _ = real_data[d_mm]
            real_width = fit_envelope(x_real, I_real)
        
        rows.append({
            'Slit Sep. $d$ (mm)': d_mm,
            'True Slit Width (mm)': true_width,
            'Fitted (Sim) (mm)': f"{fitted_width:.3f}" if fitted_width else "—",
            'Fitted (Exp) (mm)': f"{real_width:.3f}" if real_width else "—",
            '% Error (Sim)': f"{abs(fitted_width - true_width)/true_width*100:.1f}\%" if fitted_width else "—"
        })
    return pd.DataFrame(rows)

# ===============================
# 8. PLOTTING FUNCTIONS
# ===============================
def plot_interference_patterns(sim_data, real_data=None, save_path="figs/"):
    os.makedirs(save_path, exist_ok=True)
    plt.figure(figsize=(6.5, 4))
    colors = ['#1f77b4', '#d62728']
    for i, (d_mm, (x, I, env)) in enumerate(sim_data.items()):
        plt.plot(x, I, label=f'Sim $d={d_mm}$ mm', color=colors[i], lw=1.5)
        if real_data and d_mm in real_data:
            xr, Ir, _ = real_data[d_mm]
            plt.plot(xr, Ir + 0.05, '--', label=f'Exp $d={d_mm}$ mm', color=colors[i], alpha=0.8)
    
    plt.xlabel('Position on Screen (mm)')
    plt.ylabel('Normalized Intensity')
    plt.title('Double-Slit Interference Patterns')
    plt.legend(fontsize=9)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{save_path}interference_comparison.pdf", dpi=300)
    plt.savefig(f"{save_path}interference_comparison.png", dpi=300)
    plt.close()

def plot_2d_pattern(d=0.125e-3, wavelength=632.8e-9, L=1.50, pixels=512):
    x = np.linspace(-0.025, 0.025, pixels)
    y = np.linspace(-0.015, 0.015, pixels)
    X, Y = np.meshgrid(x, y)
    r1 = np.sqrt(L**2 + (X - d/2)**2 + Y**2)
    r2 = np.sqrt(L**2 + (X + d/2)**2 + Y**2)
    delta = (2*np.pi / wavelength) * (r2 - r1)
    I = np.cos(delta/2)**2
    I /= I.max()
    
    plt.figure(figsize=(5, 3.5))
    plt.imshow(I, extent=[x.min()*1e3, x.max()*1e3, y.min()*1e3, y.max()*1e3],
               cmap='hot', aspect='auto', origin='lower')
    plt.colorbar(label='Normalized Intensity')
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    plt.title(f'2D Interference Pattern ($d={d*1e3:.3f}$ mm)')
    plt.tight_layout()
    plt.savefig("figs/pattern_2d.pdf", dpi=300)
    plt.savefig("figs/pattern_2d.png", dpi=300)
    plt.close()

def plot_visibility_vs_d(sim_data, real_data=None, save_path="figs/"):
    os.makedirs(save_path, exist_ok=True)
    d_list = []
    vis_sim = []
    vis_real = []
    for d_mm in sim_data.keys():
        x, I, _ = sim_data[d_mm]
        _, _, vis, _ = find_fringe_spacing(x, I)
        d_list.append(d_mm)
        vis_sim.append(vis)
        if real_data and d_mm in real_data:
            xr, Ir, _ = real_data[d_mm]
            _, _, vis_r, _ = find_fringe_spacing(xr, Ir)
            vis_real.append(vis_r)
        else:
            vis_real.append(None)
    
    plt.figure(figsize=(5, 3.5))
    plt.plot(d_list, vis_sim, 'o-', label='Simulated', color='blue')
    if any(v is not None for v in vis_real):
        plt.plot(d_list, [v for v in vis_real if v is not None], 's--', label='Measured', color='red')
    plt.xlabel('Slit Separation $d$ (mm)')
    plt.ylabel('Fringe Visibility')
    plt.title('Fringe Visibility vs Slit Separation')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{save_path}visibility_vs_d.pdf", dpi=300)
    plt.savefig(f"{save_path}visibility_vs_d.png", dpi=300)
    plt.close()

def plot_envelope_analysis(sim_data, real_data=None, save_path="figs/"):
    os.makedirs(save_path, exist_ok=True)
    d_mm = min(sim_data.keys())  # pick smallest d for clarity
    x, I, envelope = sim_data[d_mm]
    
    plt.figure(figsize=(6, 3.5))
    plt.plot(x, I, label='Total Intensity', color='black')
    plt.plot(x, envelope, '--', label='Diffraction Envelope', color='orange', lw=2)
    
    if real_data and d_mm in real_data:
        xr, Ir, _ = real_data[d_mm]
        plt.plot(xr, Ir, ':', label='Measured', color='green', alpha=0.7)
    
    # Fit and plot
    fitted_w = fit_envelope(x, envelope)
    if fitted_w:
        x_fit = np.linspace(x.min(), x.max(), 1000)
        y_fit = sinc_envelope(x_fit, 1.0, fitted_w)
        plt.plot(x_fit, y_fit, '-', label=f'Fit ($w={fitted_w:.3f}$ mm)', color='red')
    
    plt.xlabel('Position on Screen (mm)')
    plt.ylabel('Normalized Intensity')
    plt.title(f'Envelope Analysis ($d={d_mm}$ mm)')
    plt.legend(fontsize=8)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{save_path}envelope_analysis.pdf", dpi=300)
    plt.savefig(f"{save_path}envelope_analysis.png", dpi=300)
    plt.close()

# ===============================
# 9. MAIN EXECUTION
# ===============================
if __name__ == "__main__":
    # Parameters
    wavelength_nm = 632.8
    L_distances = 1.50
    d_list_mm = [0.125, 0.250]
    slit_width_mm = 0.040

    # 1. Simulate
    sim_data = simulate_double_slit(
        wavelength=wavelength_nm*1e-9,
        d_list=[d*1e-3 for d in d_list_mm],
        L=L_distances,
        slit_width=slit_width_mm*1e-3
    )

    # 2. (Optional) Load real images
    real_data = {}
    # Example:
    # real_data[0.125] = analyze_image("photos/d_0.125mm.jpg", crop_roi=(100,200,1600,400), pixel_to_mm=0.026)
    # real_data[0.250] = analyze_image("photos/d_0.250mm.jpg", crop_roi=(120,180,1550,380), pixel_to_mm=0.026)

    # 3. Generate Tables
    table1 = make_fringe_table(sim_data, real_data, wavelength_nm, L_distances)
    table2 = make_visibility_table(sim_data, real_data)
    table3 = make_envelope_table(sim_data, real_data)

    # Save Tables
    os.makedirs("tables", exist_ok=True)
    table1.to_latex("tables/fringe_spacing.tex", index=False, escape=False, column_format="lccccc")
    table2.to_latex("tables/visibility.tex", index=False, escape=False, column_format="lcc")
    table3.to_latex("tables/envelope.tex", index=False, escape=False, column_format="lcccc")
    table1.to_csv("tables/fringe_spacing.csv", index=False)
    table2.to_csv("tables/visibility.csv", index=False)
    table3.to_csv("tables/envelope.csv", index=False)

    # 4. Generate Plots
    os.makedirs("figs", exist_ok=True)
    plot_interference_patterns(sim_data, real_data)
    plot_2d_pattern(d=0.125e-3)
    plot_visibility_vs_d(sim_data, real_data)
    plot_envelope_analysis(sim_data, real_data)

    # 5. Print Summary
    print("\n" + "="*60)
    print("ALL FILES GENERATED SUCCESSFULLY")
    print("="*60)
    print("TABLES:")
    print("  • tables/fringe_spacing.tex    → Table I")
    print("  • tables/visibility.tex        → Table II")
    print("  • tables/envelope.tex          → Table III")
    print("\nFIGURES:")
    print("  • figs/interference_comparison.pdf  → Fig. 1")
    print("  • figs/pattern_2d.pdf               → Fig. 2")
    print("  • figs/visibility_vs_d.pdf          → Fig. 3")
    print("  • figs/envelope_analysis.pdf        → Fig. 4")
    print("\nInsert into IEEE paper using:")
    print("  \\input{tables/fringe_spacing.tex}")
    print("  \\includegraphics{figs/interference_comparison.pdf}")
    print("="*60)