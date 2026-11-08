#ifndef MPC_MOTOR_H
#define MPC_MOTOR_H
/* MPC de velocidad para un motor DC, resuelto con OSQP embebido (osqp/, workspace.c). */
void mpc_motor_init(void);   /* llamar una vez al inicio */
int  mpc_motor_step(const double *x, double r, double Imax, double TL, double *V);
     /* x = [i; w] medido, r = referencia de w, Imax = limite de corriente, TL = torque de carga.
        Devuelve la tension V a aplicar (en *V) y el numero de iteraciones de OSQP. */
#endif
