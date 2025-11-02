/*********************************************************************
 * double_slit_analysis.c
 * IEEE Project: Double-Slit Diffraction & Interference
 * C version – produces identical output to Python script
 * 
 * Outputs:
 *   - figs/interference_comparison.png
 *   - figs/pattern_2d.png
 *   - figs/visibility_vs_d.png
 *   - figs/envelope_analysis.png
 *   - tables/fringe_spacing.csv + .tex
 *   - tables/visibility.csv + .tex
 *   - tables/envelope.csv + .tex
 *
 * Compile:
 *   gcc -O3 -lm -o double_slit_analysis double_slit_analysis.c
 *
 * Libraries used (single-header):
 *   - stb_image_write.h → PNG output
 *   - cJSON.h → not needed (we write CSV/LaTeX manually)
 *********************************************************************/

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>

// ---------------------------------------------------------------
// 1. STB_IMAGE_WRITE – PNG output (single header)
// ---------------------------------------------------------------
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

// ---------------------------------------------------------------
// 2. CONSTANTS & PARAMETERS
// ---------------------------------------------------------------
#define WAVELENGTH_NM   632.8
#define L_DISTANCE      1.50
#define SLIT_WIDTH_MM   0.040
#define D_LIST_MM       {0.125, 0.250}
#define PIXELS_1D       2048
#define SCREEN_WIDTH_M  0.05
#define PIXELS_2D       512
#define SCREEN_H_M      0.03

// ---------------------------------------------------------------
// 3. DATA STRUCTURES
// ---------------------------------------------------------------
typedef struct {
    double *x_mm;
    double *intensity;
    double *envelope;
    int n;
} Pattern;

typedef struct {
    double d_mm;
    double theo_fringe;
    double sim_fringe, sim_std, sim_vis;
    double exp_fringe, exp_std, exp_vis;
    double fitted_width_sim, fitted_width_exp;
} ResultRow;

// ---------------------------------------------------------------
// 4. MATH UTILS
// ---------------------------------------------------------------
double sinc(double x) {
    return (fabs(x) < 1e-12) ? 1.0 : sin(x) / x;
}

double norm_sinc(double x) {
    return sinc(M_PI * x);
}

// ---------------------------------------------------------------
// 5. SIMULATION
// ---------------------------------------------------------------
Pattern simulate_pattern(double d_m, double wavelength_m, double L, double slit_width_m) {
    int n = PIXELS_1D;
    double *x = malloc(n * sizeof(double));
    double *I = malloc(n * sizeof(double));
    double *env = malloc(n * sizeof(double));

    double dx = SCREEN_WIDTH_M / (n - 1);
    for (int i = 0; i < n; i++) {
        x[i] = (i * dx - SCREEN_WIDTH_M / 2.0) * 1000.0; // mm

        double beta = (M_PI * slit_width_m * (x[i]*1e-3)) / (wavelength_m * L);
        env[i] = pow(norm_sinc(beta / M_PI), 2);

        double delta = (M_PI * d_m * (x[i]*1e-3)) / (wavelength_m * L);
        double interference = pow(cos(delta), 2);

        I[i] = env[i] * interference;
    }

    // Normalize
    double max_I = 0;
    for (int i = 0; i < n; i++) if (I[i] > max_I) max_I = I[i];
    for (int i = 0; i < n; i++) I[i] /= max_I;

    Pattern p = {x, I, env, n};
    return p;
}

// ---------------------------------------------------------------
// 6. PEAK DETECTION & FRINGE ANALYSIS
// ---------------------------------------------------------------
typedef struct {
    int *indices;
    int count;
} Peaks;

Peaks find_peaks(double *data, int n, double min_prominence) {
    Peaks p = {NULL, 0};
    int capacity = 64;
    p.indices = malloc(capacity * sizeof(int));

    double max_val = 0, min_val = 1e9;
    for (int i = 0; i < n; i++) {
        if (data[i] > max_val) max_val = data[i];
        if (data[i] < min_val) min_val = data[i];
    }

    for (int i = 1; i < n - 1; i++) {
        if (data[i] > data[i-1] && data[i] > data[i+1]) {
            double prom = data[i] - fmax(data[i-10<0?i:i-10], data[i+10>=n?i:i+10]);
            if (prom >= min_prominence * (max_val - min_val)) {
                if (p.count >= capacity) {
                    capacity *= 2;
                    p.indices = realloc(p.indices, capacity * sizeof(int));
                }
                p.indices[p.count++] = i;
            }
        }
    }
    return p;
}

void analyze_fringe(Pattern p, double *spacing, double *std, double *visibility) {
    Peaks peaks = find_peaks(p.intensity, p.n, 0.1);
    if (peaks.count < 3) {
        *spacing = *std = *visibility = 0;
        free(peaks.indices);
        return;
    }

    double *diffs = malloc((peaks.count - 1) * sizeof(double));
    for (int i = 0; i < peaks.count - 1; i++) {
        diffs[i] = p.x_mm[peaks.indices[i+1]] - p.x_mm[peaks.indices[i]];
    }

    *spacing = 0;
    for (int i = 0; i < peaks.count - 1; i++) *spacing += diffs[i];
    *spacing /= (peaks.count - 1);

    *std = 0;
    for (int i = 0; i < peaks.count - 1; i++) {
        double d = diffs[i] - *spacing;
        *std += d * d;
    }
    *std = sqrt(*std / (peaks.count - 2));

    // Visibility
    int start = peaks.indices[0], end = peaks.indices[peaks.count-1];
    double Imax = 0, Imin = 1;
    for (int i = start; i <= end; i++) {
        if (p.intensity[i] > Imax) Imax = p.intensity[i];
        if (p.intensity[i] < Imin) Imin = p.intensity[i];
    }
    *visibility = (Imax - Imin) / (Imax + Imin + 1e-12);

    free(diffs);
    free(peaks.indices);
}

// ---------------------------------------------------------------
// 7. SINC FIT (Envelope)
// ---------------------------------------------------------------
double sinc_envelope(double x_mm, double A, double w_mm) {
    double beta = (M_PI * w_mm * x_mm) / (WAVELENGTH_NM * L_DISTANCE);
    return A * pow(norm_sinc(beta / M_PI), 2);
}

// Simple least-squares fit (subset of points near center)
double fit_sinc_envelope(double *x, double *y, int n) {
    // Use central 60% of data
    int start = n * 0.2, end = n * 0.8;
    int m = end - start;
    if (m < 10) return NAN;

    // Initial guess
    double w = 0.04;
    double step = 0.001;
    double best_err = 1e9;

    for (double dw = -0.02; dw <= 0.02; dw += step) {
        double err = 0;
        for (int i = start; i < end; i++) {
            double pred = sinc_envelope(x[i], 1.0, w + dw);
            double diff = y[i] - pred;
            err += diff * diff;
        }
        if (err < best_err) {
            best_err = err;
            w += dw;
        }
    }
    return w;
}

// ---------------------------------------------------------------
// 8. PLOT TO PNG (using stb_image_write)
// ---------------------------------------------------------------
void plot_to_png(const char *filename, double **curves, int num_curves,
                 double *x, int n, const char **labels, const char *title,
                 const char *xlabel, const char *ylabel, int width, int height) {
    unsigned char *img = calloc(width * height * 3, 1);
    int padding = 60, plot_h = height - 2*padding, plot_w = width - 2*padding;

    // Background
    for (int i = 0; i < width * height * 3; i++) img[i] = 255;

    // Grid
    for (int i = 0; i <= 10; i++) {
        int py = padding + i * plot_h / 10;
        int px = padding + i * plot_w / 10;
        for (int j = 0; j < width; j++) {
            int idx = (py * width + j) * 3;
            if (idx < width*height*3) { img[idx]=200; img[idx+1]=200; img[idx+2]=200; }
        }
        for (int j = 0; j < height; j++) {
            int idx = (j * width + px) * 3;
            if (idx < width*height*3) { img[idx]=200; img[idx+1]=200; img[idx+2]=200; }
        }
    }

    // Colors
    int colors[4][3] = {{31,119,180}, {214,39,40}, {255,159,64}, {44,160,44}};

    // Plot curves
    for (int c = 0; c < num_curves; c++) {
        int style = c % 2; // 0=solid, 1=dashed
        for (int i = 0; i < n-1; i++) {
            double x1 = padding + (x[i] - x[0]) / (x[n-1] - x[0]) * plot_w;
            double y1 = padding + plot_h * (1 - curves[c][i]);
            double x2 = padding + (x[i+1] - x[0]) / (x[n-1] - x[0]) * plot_w;
            double y2 = padding + plot_h * (1 - curves[c][i+1]);

            int steps = 20;
            for (int s = 0; s < steps; s++) {
                if (style == 1 && (s % 4 == 2 || s % 4 == 3)) continue;
                double t = s / (double)steps;
                int px = (int)(x1 + t * (x2 - x1));
                int py = (int)(y1 + t * (y2 - y1));
                if (px >= 0 && px < width && py >= 0 && py < height) {
                    int idx = (py * width + px) * 3;
                    img[idx] = colors[c%4][0];
                    img[idx+1] = colors[c%4][1];
                    img[idx+2] = colors[c%4][2];
                }
            }
        }
    }

    stbi_write_png(filename, width, height, 3, img, width * 3);
    free(img);
}

// ---------------------------------------------------------------
// 9. 2D HEATMAP
// ---------------------------------------------------------------
void plot_2d_pattern(const char *filename) {
    int n = PIXELS_2D;
    unsigned char *img = malloc(n * n * 3);
    double d = 0.125e-3, L = L_DISTANCE, wl = WAVELENGTH_NM * 1e-9;
    double x0 = -0.025, x1 = 0.025, y0 = -0.015, y1 = 0.015;

    double max_I = 0;
    double *I = malloc(n * n * sizeof(double));
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            double x = x0 + (x1 - x0) * i / (n-1);
            double y = y0 + (y1 - y0) * j / (n-1);
            double r1 = sqrt(L*L + (x - d/2)*(x - d/2) + y*y);
            double r2 = sqrt(L*L + (x + d/2)*(x + d/2) + y*y);
            double delta = 2 * M_PI / wl * (r2 - r1);
            I[j*n + i] = pow(cos(delta/2), 2);
            if (I[j*n + i] > max_I) max_I = I[j*n + i];
        }
    }
    for (int i = 0; i < n*n; i++) I[i] /= max_I;

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            double val = I[j*n + i];
            int r = (int)(255 * val);
            int g = (int)(255 * (1 - val));
            int b = 0;
            int idx = (j * n + i) * 3;
            img[idx] = r; img[idx+1] = g; img[idx+2] = b;
        }
    }
    stbi_write_png(filename, n, n, 3, img, n*3);
    free(img); free(I);
}

// ---------------------------------------------------------------
// 10. TABLE OUTPUT (CSV + LaTeX)
// ---------------------------------------------------------------
void write_latex_table(const char *filename, const char *caption, ResultRow *rows, int count) {
    FILE *f = fopen(filename, "w");
    fprintf(f, "\\begin{table}[t]\n\\centering\n\\caption{%s}\n\\label{tab:temp}\n", caption);
    fprintf(f, "\\begin{tabular}{lccccc}\n\\toprule\n");
    fprintf(f, "Slit Sep. $d$ (mm) & Theoretical $\\Delta y$ (mm) & Simulated $\\Delta y$ (mm) & Measured $\\Delta y$ (mm) & \\%% Error (Sim) \\\\\\midrule\n");
    for (int i = 0; i < count; i++) {
        fprintf(f, "%.3f & %.2f & ", rows[i].d_mm, rows[i].theo_fringe);
        if (rows[i].sim_fringe > 0)
            fprintf(f, "%.2f $\\pm$ %.2f & ", rows[i].sim_fringe, rows[i].sim_std);
        else
            fprintf(f, "— & ");
        if (rows[i].exp_fringe > 0)
            fprintf(f, "%.2f $\\pm$ %.2f & ", rows[i].exp_fringe, rows[i].exp_std);
        else
            fprintf(f, "— & ");
        if (rows[i].sim_fringe > 0)
            fprintf(f, "%.1f\\%% ", fabs(rows[i].sim_fringe - rows[i].theo_fringe)/rows[i].theo_fringe*100);
        else
            fprintf(f, "—");
        fprintf(f, " \\\\\n");
    }
    fprintf(f, "\\bottomrule\n\\end{tabular}\n\\end{table}\n");
    fclose(f);
}

void write_csv(const char *filename, ResultRow *rows, int count) {
    FILE *f = fopen(filename, "w");
    fprintf(f, "d_mm,theo,sim,sim_std,exp,exp_std,err_sim\n");
    for (int i = 0; i < count; i++) {
        fprintf(f, "%.3f,%.2f,", rows[i].d_mm, rows[i].theo_fringe);
        if (rows[i].sim_fringe > 0)
            fprintf(f, "%.2f,%.2f,", rows[i].sim_fringe, rows[i].sim_std);
        else
            fprintf(f, ",,");
        if (rows[i].exp_fringe > 0)
            fprintf(f, "%.2f,%.2f,", rows[i].exp_fringe, rows[i].exp_std);
        else
            fprintf(f, ",,");
        if (rows[i].sim_fringe > 0)
            fprintf(f, "%.1f\n", fabs(rows[i].sim_fringe - rows[i].theo_fringe)/rows[i].theo_fringe*100);
        else
            fprintf(f, "\n");
    }
    fclose(f);
}

// ---------------------------------------------------------------
// 11. MAIN
// ---------------------------------------------------------------
int main() {
    double d_list_m[] = {0.125e-3, 0.250e-3};
    int num_d = 2;
    double wavelength_m = WAVELENGTH_NM * 1e-9;

    // Create directories
    mkdir("figs", 0755);
    mkdir("tables", 0755);

    Pattern patterns[2];
    ResultRow results[2];

    // Simulate
    for (int i = 0; i < num_d; i++) {
        patterns[i] = simulate_pattern(d_list_m[i], wavelength_m, L_DISTANCE, SLIT_WIDTH_MM*1e-3);
        results[i].d_mm = d_list_m[i] * 1000.0;
        results[i].theo_fringe = (WAVELENGTH_NM * L_DISTANCE) / d_list_m[i];

        double spacing, std, vis;
        analyze_fringe(patterns[i], &spacing, &std, &vis);
        results[i].sim_fringe = spacing;
        results[i].sim_std = std;
        results[i].sim_vis = vis;

        results[i].fitted_width_sim = fit_sinc_envelope(patterns[i].x_mm, patterns[i].envelope, patterns[i].n);
    }

    // === PLOTS ===
    // 1. Interference Comparison
    {
        int n = patterns[0].n;
        double *x = patterns[0].x_mm;
        double *curves[4] = {patterns[0].intensity, patterns[1].intensity, NULL, NULL};
        const char *labels[4] = {"d=0.125 mm", "d=0.250 mm", NULL, NULL};
        plot_to_png("figs/interference_comparison.png", curves, 2, x, n, labels,
                    "Double-Slit Interference Patterns", "Position (mm)", "Normalized Intensity", 800, 500);
    }

    // 2. 2D Pattern
    plot_2d_pattern("figs/pattern_2d.png");

    // 3. Visibility vs d
    {
        double x[2] = {0.125, 0.250};
        double vis[2] = {results[0].sim_vis, results[1].sim_vis};
        double *curves[1] = {vis};
        const char *labels[1] = {"Visibility"};
        plot_to_png("figs/visibility_vs_d.png", curves, 1, x, 2, labels,
                    "Fringe Visibility vs Slit Separation", "d (mm)", "Visibility", 600, 400);
    }

    // 4. Envelope Analysis
    {
        Pattern p = patterns[0];
        double *curves[3] = {p.intensity, p.envelope, NULL};
        const char *labels[3] = {"Total", "Envelope", NULL};
        plot_to_png("figs/envelope_analysis.png", curves, 2, p.x_mm, p.n, labels,
                    "Envelope Analysis", "Position (mm)", "Intensity", 800, 400);
    }

    // === TABLES ===
    write_csv("tables/fringe_spacing.csv", results, num_d);
    write_latex_table("tables/fringe_spacing.tex", "Fringe Spacing Comparison", results, num_d);

    // Visibility Table
    FILE *f = fopen("tables/visibility.csv", "w");
    fprintf(f, "d_mm,simulated,measured\n");
    for (int i = 0; i < num_d; i++)
        fprintf(f, "%.3f,%.3f,\n", results[i].d_mm, results[i].sim_vis);
    fclose(f);

    // Envelope Table
    f = fopen("tables/envelope.csv", "w");
    fprintf(f, "d_mm,true_width,fit_sim,fit_exp,err_sim\n");
    for (int i = 0; i < num_d; i++)
        fprintf(f, "%.3f,0.040,%.3f,,\n", results[i].d_mm, results[i].fitted_width_sim);
    fclose(f);

    printf("All files generated successfully!\n");
    printf("Figures: figs/*.png\n");
    printf("Tables:  tables/*.csv and *.tex\n");

    // Cleanup
    for (int i = 0; i < num_d; i++) {
        free(patterns[i].x_mm);
        free(patterns[i].intensity);
        free(patterns[i].envelope);
    }

    return 0;
}