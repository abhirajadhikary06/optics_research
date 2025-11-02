for (int i = 0; i < N; i++) {
    x = (i - N/2) * dx;                     // position in mm
    beta  = M_PI * a * x / (lambda * L);
    delta = M_PI * d * x / (lambda * L);
    envelope = pow(sinc(beta / M_PI), 2);
    I[i] = envelope * pow(cos(delta), 2);
}
I[i] /= I_max;  // normalize