/* MPC de velocidad para un motor DC con OSQP embebido.
 * En cada muestra: (1) estado estacionario para la referencia, (2) limites l, u del QP,
 * (3) OSQP resuelve, (4) se aplica la primera tension.
 * El QP (matrices P y A) esta en workspace.c, generado con OSQP codegen (generar_mpc_motor.py).
 * Los valores de abajo deben coincidir con los usados para generar workspace.c. */

#include "mpc_motor.h"
#include "osqp.h"
#include "workspace.h"

#define N      10      /* horizonte                                    */
#define N_VAR  32      /* variables: x~0..x~N (22) + u~0..u~(N-1) (10)  */
#define M_RES  42      /* restricciones: 22 dinamica + 10 i + 10 V      */
#define R_M    1.0     /* resistencia de armadura [ohm]                 */
#define K_M    0.01    /* constante del motor [V s/rad] = [N m/A]       */
#define B_M    0.1     /* friccion viscosa [N m s]                      */
#define V_MAX  24.0    /* |V| <= V_MAX [V]                              */

void mpc_motor_init(void)
{
    c_float x0[N_VAR] = {0}, y0[M_RES] = {0};
    osqp_warm_start(&workspace, x0, y0);          /* parte desde cero */
}

int mpc_motor_step(const double *x, double r, double Imax, double TL, double *V)
{
    c_float l[M_RES], u[M_RES];
    int k;

    /* (1) Estado estacionario: con w = r,  K i = b w + TL  y  V = R i + K w */
    double w_inf = r;
    double i_inf = (B_M * r + TL) / K_M;
    double V_inf = R_M * i_inf + K_M * r;

    /* (2) Limites en variables de desviacion (x~ = x - x_inf, u~ = V - V_inf) */
    l[0] = u[0] = -(x[0] - i_inf);                 /* x~0 = x medido - x_inf     */
    l[1] = u[1] = -(x[1] - w_inf);

    for (k = 2; k < 2 * (N + 1); k++)
    {    
        l[k] = u[k] = 0.0; /* dinamica: igualdades */
    }
    
    for (k = 0; k < N; k++) 
    { 
        l[22 + k] = -Imax  - i_inf;
        u[22 + k] = Imax  - i_inf; /* |i| <= Imax */
    }  

    for (k = 0; k < N; k++) {
        l[32 + k] = -V_MAX - V_inf;
        u[32 + k] = V_MAX - V_inf; /* |V| <= Vmax */
    }  

    /* (3) Resolver el QP */
    osqp_update_bounds(&workspace, l, u);
    osqp_solve(&workspace);

    /* (4) Primera tension del horizonte, de vuelta a variable absoluta */
    *V = workspace.solution->x[2 * (N + 1)] + V_inf;
    return (int)workspace.info->iter;
}
