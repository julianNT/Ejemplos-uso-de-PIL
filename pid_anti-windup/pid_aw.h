#ifndef PID_AW_H
#define PID_AW_H
/*
 * PID paralelo con derivada filtrada y anti-windup por back-calculation,
 * discretizado con Forward Euler (s -> (z-1)/Ts).
 *
 *   C(s) = Kp + Ki/s + Kd*s/(Tf*s + 1),   dI/dt = Ki*e + Kb*(u_sat - v),   Kb = 1/Tt
 *
 * Vector de parametros p (8 elementos):
 *   p[0] Kp   p[1] Ki   p[2] Kd   p[3] Tf   p[4] Kb   p[5] Ts   p[6] umin   p[7] umax
 * Kb = 0 desactiva el anti-windup.
 */
void   pid_aw_init(void);                                     /* estados en cero: I = 0, xf = 0                */
double pid_aw_step(double r, double y, const double *p,
                   double *v);                                /* devuelve u_sat[k]; *v = u[k] sin saturar      */
#endif
