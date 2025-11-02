void plot_to_png(const char *filename, double **curves, int num_curves,
                 double *x, int n, const char **labels, const char *title,
                 const char *xlabel, const char *ylabel, int width, int height) {
    unsigned char *img = calloc(width * height * 3, 1);
    int padding = 60, plot_h = height - 2*padding, plot_w = width - 2*padding;
    for (int i = 0; i < width * height * 3; i++) img[i] = 255; // White background
    // Grid and curve plotting code follows (omitted for brevity)
    stbi_write_png(filename, width, height, 3, img, width * 3);
    free(img);
}