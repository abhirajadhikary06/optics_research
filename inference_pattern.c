Pattern simulate_pattern(double d_m, double wavelength_m, double L, double slit_width_m) {
    int n = PIXELS_1D;
    double *x = malloc(n * sizeof(double));
    double *I = malloc(n * sizeof(double));
    double *env = malloc(n * sizeof(double));
    double dx = SCREEN_WIDTH_M / (n - 1);
    for (int i = 0; i < n; i++) {
        x[i] = (i * dx - SCREEN_WIDTH_M / 2.0) * 1000.0; // mm
        double beta = (M_PI * slit_width_m * (x[i]*1e-3)) / (wavelength_m * L);
        env[i] = pow(sinc(beta / M_PI), 2);
        double delta = (M_PI * d_m * (x[i]*1e-3)) / (wavelength_m * L);
        I[i] = env[i] * pow(cos(delta), 2);
    }
    double max_I = 0;
    for (int i = 0; i < n; i++) if (I[i] > max_I) max_I = I[i];
    for (int i = 0; i < n; i++) I[i] /= max_I;
    Pattern p = {x, I, env, n};
    return p;
}