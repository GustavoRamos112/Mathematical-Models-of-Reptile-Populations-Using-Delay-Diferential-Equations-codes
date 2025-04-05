import numpy as np
import matplotlib.pyplot as plt
from ddeint import ddeint

# 1. VALORES DE PARÁMETROS 
tau_2 = 8 * b2   # Retraso 1 (8*b2)
tau_3 = 11 * b2  # Retraso 2 (11*b2)
tau_5 = 31 * b2  # Retraso 3 (31*b2)
retrasos = [tau_2, tau_3, tau_5]  # Lista de retrasos

# 2. INICIALIZACIÓN DE VARIABLES
t0 = 0       # Tiempo inicial
tF = 200     # Tiempo final
b2 = 0.844   # Tasa de natalidad

# 3. SOLUCIÓN NUMÉRICA DE EDDR
def historia(t):
    """Función de historia para valores iniciales constantes"""
    return np.array([
        10, 10, 10, 10,    # Af1, Af2, Am2, Am3
        10, 7, 4, 4,       # Bf1, Bf2, Bm2, Bm3
        2.5, 1.5, 1, 1,    # Cf1, Cf2, Cm2, Cm3
        2.5, 1.5, 1, 1     # Df1, Df2, Dm2, Dm3
    ])

def modelo(Y, t):
    """Sistema de ecuaciones diferenciales retardadas"""
    # Obtener estados retardados
    Y_tau2 = Y(t - retrasos[0])
    Y_tau3 = Y(t - retrasos[1])
    Y_tau5 = Y(t - retrasos[2])
    y = Y(t)  # Estado actual
    
    # Parámetros del sistema
    d1, d2, d_new, d3 = 0.245, 0.151, 0.001, 0.139
    b1, k1 = 0.286, 0.79*1.5
    k2, k3 = 0.15*1.5, 0.06*1.5
    s8 = np.exp(-d1*8)
    s11 = np.exp(-d2*3)
    s31 = np.exp(-d_new*20)
    
    # Cálculo de parámetros derivados
    c1, c2 = k2/k1, k3/k1
    r1, a1 = b1/b2, d1/b2
    a2, a_new = d2/b2, d_new/b2
    a3 = d3/b2
    
    # Variables auxiliares (Q, R, S)
    Q1 = (y[4] + y[8] + y[12]) / (1 + y[4] + y[8] + y[12])
    Q2 = c1 / (c1 + y[4] + y[5] + y[8] + y[9] + y[12] + y[13])
    Q3 = (Y_tau2[4] + Y_tau2[8] + Y_tau2[12]) / (1 + Y_tau2[4] + Y_tau2[8] + Y_tau2[12])
    Q4 = c1 / (c1 + Y_tau2[4] + Y_tau2[5] + Y_tau2[8] + Y_tau2[9] + Y_tau2[12] + Y_tau2[13])
    Q5 = c2 / (c2 + y[4] + y[5] + y[8] + y[9] + y[12] + y[13])
    Q6 = c2 / (c2 + Y_tau2[4] + Y_tau2[5] + Y_tau2[8] + Y_tau2[9] + Y_tau2[12] + Y_tau2[13])
    
    R3 = (Y_tau3[4] + Y_tau3[8] + Y_tau3[12]) / (1 + Y_tau3[4] + Y_tau3[8] + Y_tau3[12])
    R4 = c1 / (c1 + Y_tau3[4] + Y_tau3[5] + Y_tau3[8] + Y_tau3[9] + Y_tau3[12] + Y_tau3[13])
    R6 = c2 / (c2 + Y_tau3[4] + Y_tau3[5] + Y_tau3[8] + Y_tau3[9] + Y_tau3[12] + Y_tau3[13])
    
    S3 = (Y_tau5[4] + Y_tau5[8] + Y_tau5[12]) / (1 + Y_tau5[4] + Y_tau5[8] + Y_tau5[12])
    S4 = c1 / (c1 + Y_tau5[4] + Y_tau5[5] + Y_tau5[8] + Y_tau5[9] + Y_tau5[12] + Y_tau5[13])
    S6 = c2 / (c2 + Y_tau5[4] + Y_tau5[5] + Y_tau5[8] + Y_tau5[9] + Y_tau5[12] + Y_tau5[13])

    # Sistema de ecuaciones diferenciales
    dydt = np.zeros(16)
    
    # Población 0-8 años 
    # R1 hembras 0-8: y[0]
    dydt[0] = (r1*(y[4]*(1 - Q1)) 
        - s8*r1*Y_tau2[4]*(1 - Q3) 
        - a1*y[0] 
        + (y[8] + y[12])*(1 - Q1) 
        - s8*(Y_tau2[8] + Y_tau2[12])*(1 - Q3))
    # R2 hembras 0-8: y[1]
    dydt[1] = (0.5*r1*((y[4] + y[5])*Q1*Q2 
        - s8*(Y_tau2[4] + Y_tau2[5])*Q3*Q4) 
        - a1*y[1] 
        + 0.5*((y[8] + y[9] + y[12] + y[13])*Q1*Q2 
        - s8*(Y_tau2[8] + Y_tau2[9] + Y_tau2[12] + Y_tau2[13])*Q3*Q4))
    # R2 machos 0-8: y[2]
    dydt[2] = (0.5*r1*((y[4] + y[5])*Q1*Q2 
        - s8*(Y_tau2[4] + Y_tau2[5])*Q3*Q4) 
        - a1*y[2] 
        + 0.5*((y[8] + y[9] + y[12] + y[13])*Q1*Q2 
        - s8*(Y_tau2[8] + Y_tau2[9] + Y_tau2[12] + Y_tau2[13])*Q3*Q4))
    # R3 machos 0-8: y[3]
    dydt[3] = (r1*(Q5*(y[4]*Q1 + y[5])*(1 - Q2) 
        - s8*Q6*(Y_tau2[4]*Q3 + Y_tau2[5])*(1 - Q4)) 
        - a1*y[3] 
        + Q5*((y[8] + y[12])*Q1 + y[9] + y[13])*(1 - Q2) 
        - s8*Q6*((Y_tau2[8] + Y_tau2[12])*Q3 + Y_tau2[9] + Y_tau2[13])*(1 - Q4))

    # Población 8-11 años
    # R1 hembras 8-11: y[4]
    dydt[4] = (s8*r1*Y_tau2[4]*(1 - Q3) 
        + s8*(Y_tau2[8] + Y_tau2[12])*(1 - Q3) 
        - a2*y[4] 
        - s8*s11*r1*Y_tau3[4]*(1 - R3) 
        - s8*s11*(Y_tau3[8] + Y_tau3[12])*(1 - R3))
    # R2 hembras 8-11: y[5]
    dydt[5] = (s8*0.5*r1*(Y_tau2[4] + Y_tau2[5])*Q3*Q4 
        + s8*(Y_tau2[8] + Y_tau2[9] + Y_tau2[12] + Y_tau2[13])*Q3*Q4 
        - a2*y[5] 
        - s8*s11*0.5*r1*(Y_tau3[4] + Y_tau3[5])*R3*R4 
        - s8*s11*(Y_tau3[8] + Y_tau3[9] + Y_tau3[12] + Y_tau3[13])*R3*R4)
    # R2 machos 8-11: y[6]
    dydt[6] = (s8*0.5*r1*(Y_tau2[4] + Y_tau2[5])*Q3*Q4 
        + s8*(Y_tau2[8] + Y_tau2[9] + Y_tau2[12] + Y_tau2[13])*Q3*Q4 
        - a2*y[6] 
        - s8*s11*0.5*r1*(Y_tau3[4] + Y_tau3[5])*R3*R4 
        - s8*s11*(Y_tau3[8] + Y_tau3[9] + Y_tau3[12] + Y_tau3[13])*R3*R4)
    # R3 machos 8-11: y[7]
    dydt[7] = (s8*Q6*r1*(Y_tau2[4]*Q3 + Y_tau2[5])*(1 - Q4) 
        + s8*Q6*((Y_tau2[8] + Y_tau2[12])*Q3 + Y_tau2[9] + Y_tau2[13])*(1 - Q4) 
        - a2*y[7] 
        - s8*s11*R6*r1*(Y_tau3[4]*R3 + Y_tau3[5])*(1 - R4) 
        - s8*s11*R6*((Y_tau3[8] + Y_tau3[12])*R3 + Y_tau3[9] + Y_tau3[13])*(1 - R4))

    # Población 11-31 años #
    # R1 hembras 11-31: y[8]
    dydt[8] = (s8*s11*r1*Y_tau3[4]*(1 - R3) 
        + s8*s11*(Y_tau3[8] + Y_tau3[12])*(1 - R3) 
        - a_new*y[8] 
        - s8*s11*s31*r1*Y_tau5[4]*(1 - S3) 
        - s8*s11*s31*(Y_tau5[8] + Y_tau5[12])*(1 - S3))
    # R2 hembras 11-31: y[9]
    dydt[9] = (s8*s11*0.5*r1*(Y_tau3[4] + Y_tau3[5])*R3*R4 
        + s8*s11*(Y_tau3[8] + Y_tau3[9] + Y_tau3[12] + Y_tau3[13])*R3*R4 
        - a_new*y[9] 
        - s8*s11*s31*0.5*r1*(Y_tau5[4] + Y_tau5[5])*S3*S4 
        - s8*s11*s31*(Y_tau5[8] + Y_tau5[9] + Y_tau5[12] + Y_tau5[13])*S3*S4)
    # R2 machos 11-31: y[10]
    dydt[10] = (s8*s11*0.5*r1*(Y_tau3[4] + Y_tau3[5])*R3*R4 
        + s8*s11*(Y_tau3[8] + Y_tau3[9] + Y_tau3[12] + Y_tau3[13])*R3*R4 
        - a_new*y[10] 
        - s8*s11*s31*0.5*r1*(Y_tau5[4] + Y_tau5[5])*S3*S4 
        - s8*s11*s31*(Y_tau5[8] + Y_tau5[9] + Y_tau5[12] + Y_tau5[13])*S3*S4)
    # R3 machos 11-31: y[11]
    dydt[11] = (s8*s11*R6*r1*(Y_tau3[4]*R3 + Y_tau3[5])*(1 - R4) 
        + s8*s11*R6*((Y_tau3[8] + Y_tau3[12])*R3 + Y_tau3[9] + Y_tau3[13])*(1 - R4) 
        - a_new*y[11] 
        - s8*s11*s31*S6*r1*(Y_tau5[4]*S3 + Y_tau5[5])*(1 - S4) 
        - s8*s11*s31*S6*((Y_tau5[8] + Y_tau5[12])*S3 + Y_tau5[9] + Y_tau5[13])*(1 - S4))

    # Población 31+ años
    # R1 hembras 31+: y[12]
    dydt[12] = (s8*s11*s31*r1*Y_tau5[4]*(1 - S3) 
        + s8*s11*s31*(Y_tau5[8] + Y_tau5[12])*(1 - S3) 
        - a3*y[12])
    # R2 hembras 31+: y[13]
    dydt[13] = (s8*s11*s31*0.5*r1*(Y_tau5[4] + Y_tau5[5])*S3*S4 
        + s8*s11*s31*(Y_tau5[8] + Y_tau5[9] + Y_tau5[12] + Y_tau5[13])*S3*S4 
        - a3*y[13])
    # R2 machos 31+: y[14]
    dydt[14] = (s8*s11*s31*0.5*r1*(Y_tau5[4] + Y_tau5[5])*S3*S4 
        + s8*s11*s31*(Y_tau5[8] + Y_tau5[9] + Y_tau5[12] + Y_tau5[13])*S3*S4 
        - a3*y[14])
    # R3 machos 31+: y[15]
    dydt[15] = (s8*s11*s31*S6*r1*(Y_tau5[4]*S3 + Y_tau5[5])*(1 - S4) 
        + s8*s11*s31*S6*((Y_tau5[8] + Y_tau5[12])*S3 + Y_tau5[9] + Y_tau5[13])*(1 - S4) 
        - a3*y[15])

    return dydt

# 4. Resolver el sistema
t = np.linspace(t0, tF, 2000)
sol = ddeint(modelo, historia, t)

# Extraer variables
Af1 = sol[:, 0]    # Hembras 0-8 R1
Af2 = sol[:, 1]    # Hembras 0-8 R2
Am2 = sol[:, 2]    # Machos 0-8 R2
Am3 = sol[:, 3]    # Machos 0-8 R3
Bf1 = sol[:, 4]    # Hembras 8-11 R1
Bf2 = sol[:, 5]    # Hembras 8-11 R2
Bm2 = sol[:, 6]    # Machos 8-11 R2
Bm3 = sol[:, 7]    # Machos 8-11 R3
Cf1 = sol[:, 8]    # Hembras 11-31 R1
Cf2 = sol[:, 9]    # Hembras 11-31 R2
Cm2 = sol[:, 10]   # Machos 11-31 R2
Cm3 = sol[:, 11]   # Machos 11-31 R3
Df1 = sol[:, 12]   # Hembras 31+ R1
Df2 = sol[:, 13]   # Hembras 31+ R2
Dm2 = sol[:, 14]   # Machos 31+ R2
Dm3 = sol[:, 15]   # Machos 31+ R3
# Nota importante sobre la diferencia con MATLAB:
# En el código original de MATLAB había un error en las asignaciones de Df2, Dm2 y Dm3
# donde todas usaban el índice 13. Esto se corrigió en la versión Python usando los índices
# correctos 13, 14 y 15 respectivamente para mantener coherencia con las 16 variables.

# Cálculos adicionales
u = Af1 + Bf1 + Cf1 + Df1  # Población femenina total R1
v = Af2 + Bf2 + Cf2 + Df2  # Población femenina total R2
w = Am2 + Bm2 + Cm2 + Dm2  # Población masculina total R2
q = Am3 + Bm3 + Cm3 + Dm3  # Población masculina total R3
r = (w + q) / (u + v + w + q)  # Proporción de machos

# 5. GRÁFICOS 
# Primer gráfico: Distribución por edades
plt.figure(figsize=(12,6))
# Trazar todas las series con estilos correspondientes
plt.plot(t, Af1, 'g-', lw=3, label='Af1 (0-8 R1)')
plt.plot(t, Af2, 'r-', lw=3, label='Af2 (0-8 R2)')
plt.plot(t, Am2, 'y-', lw=3, label='Am2 (0-8 R2)')
plt.plot(t, Am3, 'b-', lw=3, label='Am3 (0-8 R3)')
plt.plot(t, Bf1, 'g:', lw=3, label='Bf1 (8-11 R1)')
plt.plot(t, Bf2, 'r:', lw=3, label='Bf2 (8-11 R2)')
plt.plot(t, Bm2, 'y:', lw=3, label='Bm2 (8-11 R2)')
plt.plot(t, Bm3, 'b:', lw=3, label='Bm3 (8-11 R3)')
plt.plot(t, Cf1, 'g--', lw=3, label='Cf1 (11-31 R1)')
plt.plot(t, Cf2, 'r--', lw=3, label='Cf2 (11-31 R2)')
plt.plot(t, Cm2, 'y--', lw=3, label='Cm2 (11-31 R2)')
plt.plot(t, Cm3, 'b--', lw=3, label='Cm3 (11-31 R3)')
plt.plot(t, Df1, 'g-.', lw=3, label='Df1 (31+ R1)')
plt.plot(t, Df2, 'r-.', lw=3, label='Df2 (31+ R2)')  
plt.plot(t, Dm2, 'y-.', lw=3, label='Dm2 (31+ R2)')  
plt.plot(t, Dm3, 'b-.', lw=3, label='Dm3 (31+ R3)')  

# Ajustar leyenda y límites
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

plt.title('Distribución de población de caimanes')
plt.xlabel('Tiempo')
plt.ylabel('Individuos')
plt.legend()

# Segundo gráfico: Poblaciones totales
plt.figure(figsize=(12,6))
plt.plot(t, u, 'g-', lw=3, label='Hembras R1')
plt.plot(t, v, 'r-', lw=3, label='Hembras R2')
plt.plot(t, w, 'y-', lw=3, label='Machos R2')
plt.plot(t, q, 'b-', lw=3, label='Machos R3')
plt.title('Poblaciones totales por región')
plt.xlabel('Tiempo')
plt.ylabel('Individuos')
plt.legend()

# Tercer gráfico: Proporción de machos
plt.figure(figsize=(12,6))
plt.plot(t, r, 'r-', lw=3)
plt.title('Proporción de machos en la población')
plt.xlabel('Tiempo')
plt.ylabel('Proporción')

plt.show()