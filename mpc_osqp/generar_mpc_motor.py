"""
generar_mpc_motor.py  -  Diseña el MPC del motor DC y regenera osqp/src/workspace.c y osqp/include/workspace.h

Uso (en esta carpeta):
    python generar_mpc_motor.py               # regenera workspace.c/.h y muestra qué ajustar en mpc_motor.c
    python generar_mpc_motor.py --referencia  # además guarda referencia_python.csv (lazo cerrado en Python)

Requiere: pip install "osqp==0.6.7.post3" numpy scipy

QUÉ HACE
  1. Arma el QP del MPC (matrices P y A, vector q) a partir de los parámetros de abajo.
  2. Con el codegen de OSQP genera workspace.c/.h y los copia a osqp/src y osqp/include.
     El resto de osqp/ (copia exacta del repo TiempoMPC, OSQP 0.6.2) NO se toca.

QUÉ PARÁMETROS REQUIEREN REGENERAR workspace.c (y volver a compilar / correr PIL)
  planta (J, b, K, R, L), Ts, N, Q, Rw, Vmax (si lo tocas, también V_MAX en mpc_motor.c),
  y los ajustes de OSQP (max_iter, eps_abs, eps_rel, ...).
QUÉ NO LO REQUIERE (cambian en línea): referencia, Imax, TL (entradas del bloque C Function).

DESPUÉS DE REGENERAR, EDITA A MANO mpc_motor.c si cambiaste algo de esto
  (el script imprime los valores al final):
     N, N_VAR, M_RES, R_M, K_M, B_M, V_MAX   y los índices 22 y 32 de los límites.
"""
import numpy as np, scipy.linalg as sl, scipy.signal as ss, scipy.sparse as sp, osqp, os, shutil, sys

# ======================= PARÁMETROS (edita aquí) =======================
# Planta: motor DC (modelo de velocidad CTMS, Univ. Michigan).  x = [i; w],  u = V,  p = TL
J, b, K, R, L = 0.01, 0.1, 0.01, 1.0, 0.5   # kg m^2, N m s, V/(rad/s)=N m/A, ohm, H
Ts   = 0.01                                  # periodo de muestreo [s]  (debe coincidir con el Sample time en Simulink)

# MPC
N    = 10                                    # horizonte
Q    = np.diag([0.01, 1000.0])               # peso en [i, w]
Rw   = np.array([[0.01]])                    # peso en V
Vmax = 24.0                                  # |V| <= Vmax [V]

# Ajustes de OSQP (quedan fijos dentro de workspace.c)
ajustes = dict(verbose=False, warm_start=True, adaptive_rho=False, polish=False,
               max_iter=100, check_termination=5, eps_abs=1e-3, eps_rel=1e-3)

# Escenario solo para el lazo de referencia en Python (igual al de lazo.slx)
def escenario(t):
    r  = 0.0 if t < 0.1 else 1.0
    TL = 0.0 if t < 0.8 else 0.005
    return r, TL, 11.0                       # referencia, torque de carga, Imax
# =======================================================================

A  = np.array([[-R/L, -K/L], [K/J, -b/J]])
B  = np.array([[1/L], [0.0]])
Bp = np.array([[0.0], [-1/J]])
Ad, BB, _, _, _ = ss.cont2discrete((A, np.hstack([B, Bp]), np.eye(2), np.zeros((2, 2))), Ts, 'zoh')
Bd, Bpd = BB[:, [0]], BB[:, [1]]
nx, nu = 2, 1
Qp = sl.solve_discrete_are(Ad, Bd, Q, Rw)    # costo terminal (DARE, igual a dlqr en MATLAB)

# z = [x0..xN, u0..u(N-1)] en variables de desviación (x~ = x - x_inf, u~ = V - V_inf)
P   = sp.block_diag([sp.kron(sp.eye(N), Q), Qp, sp.kron(sp.eye(N), Rw)], format='csc')
q   = np.zeros(P.shape[0])
Ax  = sp.kron(sp.eye(N + 1), -sp.eye(nx)) + sp.kron(sp.eye(N + 1, k=-1), Ad)
Bu  = sp.kron(sp.vstack([sp.csc_matrix((1, N)), sp.eye(N)]), Bd)
Aeq = sp.hstack([Ax, Bu])                                                    # 2(N+1) filas: x0 y dinámica
Ai  = sp.hstack([sp.csc_matrix((N, nx)), sp.kron(sp.eye(N), np.array([[1, 0]])),
                 sp.csc_matrix((N, N * nu))])                                # N filas: |i| en x1..xN
Au  = sp.hstack([sp.csc_matrix((N, (N + 1) * nx)), sp.eye(N * nu)])          # N filas: |V| en u0..u(N-1)
Aqp = sp.vstack([Aeq, Ai, Au], format='csc')

def estado_estacionario(r, TL):
    M = np.block([[np.eye(2) - Ad, -Bd], [np.array([[0, 1]]), np.zeros((1, 1))]])
    s = np.linalg.solve(M, np.concatenate([(Bpd * TL).ravel(), [r]]))
    return s[:2], s[2:]

def limites(x, xinf, uinf, Imax):
    leq = np.concatenate([-(x - xinf), np.zeros(N * nx)])
    return (np.concatenate([leq, np.full(N, -Imax - xinf[0]), np.full(N, -Vmax - uinf[0])]),
            np.concatenate([leq, np.full(N,  Imax - xinf[0]), np.full(N,  Vmax - uinf[0])]))

if __name__ == '__main__':
    aqui = os.path.dirname(os.path.abspath(__file__))
    l0, u0 = limites(np.zeros(2), np.zeros(2), np.zeros(1), 11.0)

    # ---- codegen de OSQP: solo se usan workspace.c y workspace.h ----
    prob = osqp.OSQP(); prob.setup(P, q, Aqp, l0, u0, **ajustes)
    gen = os.path.join(aqui, '_codegen_tmp'); shutil.rmtree(gen, ignore_errors=True)
    prob.codegen(gen, project_type='', parameters='vectors', force_rewrite=True, compile_python_ext=False)
    os.makedirs(os.path.join(aqui, 'osqp', 'src'), exist_ok=True)
    os.makedirs(os.path.join(aqui, 'osqp', 'include'), exist_ok=True)
    shutil.copyfile(os.path.join(gen, 'src', 'osqp', 'workspace.c'), os.path.join(aqui, 'osqp', 'src', 'workspace.c'))
    shutil.copyfile(os.path.join(gen, 'include', 'workspace.h'),     os.path.join(aqui, 'osqp', 'include', 'workspace.h'))
    shutil.rmtree(gen)

    n, m = P.shape[0], Aqp.shape[0]
    print('workspace.c / workspace.h regenerados en osqp/')
    print('\nValores que deben coincidir en mpc_motor.c:')
    print(f'  #define N      {N}')
    print(f'  #define N_VAR  {n}      /* 2(N+1) + N */')
    print(f'  #define M_RES  {m}      /* 2(N+1) dinamica + N corriente + N tension */')
    print(f'  #define R_M    {R}\n  #define K_M    {K}\n  #define B_M    {b}\n  #define V_MAX  {Vmax}')
    print(f'  limites de corriente: filas {2*(N+1)} .. {2*(N+1)+N-1}   (l[22 + k] en el codigo actual)')
    print(f'  limites de tension  : filas {2*(N+1)+N} .. {2*(N+1)+2*N-1}   (l[32 + k] en el codigo actual)')
    print(f'  primera tension     : solution->x[{2*(N+1)}]')
    print(f'  Ts = {Ts} s  (Sample time del bloque y de mpc_ctrl.slx)')

    # ---- lazo cerrado de referencia en Python (opcional) ----
    if '--referencia' in sys.argv:
        ref = osqp.OSQP(); ref.setup(P, q, Aqp, l0, u0, **ajustes)
        x = np.zeros(2); filas = []
        for k in range(150):
            r, TL, Im = escenario(k * Ts)
            xi, ui = estado_estacionario(r, TL); l, u = limites(x, xi, ui, Im)
            ref.update(l=l, u=u); res = ref.solve()
            V = res.x[(N + 1) * nx] + ui[0]
            filas.append([k * Ts, x[0], x[1], V, res.info.iter, r, TL])
            x = Ad @ x + (Bd * V).ravel() + (Bpd * TL).ravel()
        np.savetxt(os.path.join(aqui, 'referencia_python.csv'), np.array(filas), delimiter=',',
                   header='t,i,w,V,iter,ref,TL', comments='')
        print('referencia_python.csv guardado')
