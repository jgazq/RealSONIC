PROCEDURE update() {
        A_t = Adrive * stimon  : time-dependent amplitude for lookups
        Vm = V(A_t, v, q1, f1) : + V1(A_t, v, q1, f1) * cos (2 * PI * Fdrive/1000 * t + psi1(A_t, v, q1, f1)) : charge density dependent voltage
        y = v                  : charge density dependent rates
}