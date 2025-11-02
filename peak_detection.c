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