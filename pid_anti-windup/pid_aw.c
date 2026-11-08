#include "pid_aw.h"

/* Estados del controlador (una sola instancia: ModelReferenceNumInstancesAllowed = Single) */
static double I  = 0.0;   /* integrador          I[k]   */
static double xf = 0.0;   /* filtro derivativo   xf[k]  */

void pid_aw_init(void)
{
    I  = 0.0;
    xf = 0.0;
}

double pid_aw_step(double r, double y, const double *p, double *v)
{
    const double Kp = p[0], Ki = p[1], Kd = p[2], Tf = p[3];
    const double Kb = p[4], Ts = p[5], umin = p[6], umax = p[7];

    /* 1) error */
    const double e = r - y;

    /* 2) salida con los estados en k:  D[k] = Kd/Tf (e[k] - xf[k]),  u[k] = Kp e[k] + I[k] + D[k] */
    const double D = (Kd / Tf) * (e - xf);
    const double u = Kp * e + I + D;

    /* 3) saturacion */
    double us = u;
    if (us > umax) us = umax;
    if (us < umin) us = umin;

    /* 4) actualizacion Forward Euler de los estados a k+1 */
    I  = I  + Ts * (Ki * e + Kb * (us - u));     /* back-calculation      */
    xf = xf + (Ts / Tf) * (e - xf);              /* filtro de la derivada */

    *v = u;
    return us;
}
