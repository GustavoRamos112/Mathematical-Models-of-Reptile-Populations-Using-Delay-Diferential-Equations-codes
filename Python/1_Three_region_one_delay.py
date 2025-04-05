# Importar bibliotecas necesarias
import numpy as np
import matplotlib.pyplot as plt
from ddeint import ddeint

# 1. VALORES DE LOS PARÁMETROS
# Parámetros demográficos
b = 0.826      # Tasa de natalidad
dJ = 0.2       # Tasa de mortalidad juvenil
dA = 0.0928    # Tasa de mortalidad adulta
k1 = 0.79      # Capacidad de carga región 1
k2 = 0.15      # Capacidad de carga región 2
k3 = 0.06      # Capacidad de carga región 3

# Cálculo de parámetros derivados
a1 = dJ / b    # Tasa mortalidad juvenil normalizada
a2 = dA / b    # Tasa mortalidad adulta normalizada
c1 = k2 / k1   # Ratio capacidad región 2 vs 1
c2 = k3 / k1   # Ratio capacidad región 3 vs 1
s = np.exp(-dJ * 10)  # Probabilidad de supervivencia hasta madurez

# 2. INICIALIZACIÓN DE VARIABLES
t0 = 0         # Tiempo inicial
tF = 200       # Tiempo final
tau = 10 * b   # Retardo temporal
lags = [tau]   # Lista de retardos

# 3. SOLUCIÓN NUMÉRICA DE LAS EDDR
def historia(t):
    """Función de historia (valores constantes para t <= 0)"""
    return np.array([10.0, 10.0, 10.0, 10.0, 15.0, 10.0, 5.0, 5.0])

def modelo(Y, t):
    """Sistema de ecuaciones diferenciales con retardo"""
    # Desempacar variables actuales
    f1, f2, m2, m3, F1, F2, M2, M3 = Y(t)

    # Obtener valores retardados
    Y_ret = Y(t - tau)
    F1_ret = Y_ret[4]  # F1 retardado
    F2_ret = Y_ret[5]  # F2 retardado

    # Inicializar vector de derivadas
    dydt = np.zeros(8)

    # 1. Ecuación para juveniles hembra R1
    dydt[0] = (F1/(1 + F1)) - a1*f1 - s*F1_ret/(1 + F1_ret)

    # 2. Ecuación para juveniles hembra R2
    term = (c1/(c1 + F1 + F2)) * (F1**2/(1 + F1) + F2)
    term_ret = (c1/(c1 + F1_ret + F2_ret)) * (F1_ret**2/(1 + F1_ret) + F2_ret)
    dydt[1] = 0.5*term - a1*f2 - 0.5*s*term_ret

    # 3. Ecuación para juveniles macho R2 (igual que hembra R2)
    dydt[2] = 0.5*term - a1*m2 - 0.5*s*term_ret

    # 4. Ecuación para juveniles macho R3
    term1 = (c2/(c2 + F1 + F2)) * (F1**2/(1 + F1) + F2)
    term2 = (F1 + F2)/(c1 + F1 + F2)
    term_ret1 = (c2/(c2 + F1_ret + F2_ret)) * (F1_ret**2/(1 + F1_ret) + F2_ret)
    term_ret2 = (F1_ret + F2_ret)/(c1 + F1_ret + F2_ret)
    dydt[3] = term1*term2 - a1*m3 - s*term_ret1*term_ret2

    # 5. Ecuación para adultos hembra R1
    dydt[4] = s*F1_ret/(1 + F1_ret) - a2*F1

    # 6. Ecuación para adultos hembra R2
    dydt[5] = 0.5*s*term_ret - a2*F2

    # 7. Ecuación para adultos macho R2 (igual que hembra R2)
    dydt[6] = 0.5*s*term_ret - a2*M2

    # 8. Ecuación para adultos macho R3
    dydt[7] = s*term_ret1*term_ret2 - a2*M3

    return dydt

# Resolver el sistema
t_span = np.linspace(t0, tF, 2000)
solucion = ddeint(modelo, historia, t_span)

# 4. GRÁFICOS
# Extraer componentes de la solución
f1 = solucion[:,0]
f2 = solucion[:,1]
m2 = solucion[:,2]
m3 = solucion[:,3]
F1 = solucion[:,4]
F2 = solucion[:,5]
M2 = solucion[:,6]
M3 = solucion[:,7]

# Calcular variables compuestas
total_m = m2 + m3 + M2 + M3
total_f = f1 + f2 + F1 + F2
r = total_m / (total_m + total_f)
q = f1 + F1  # Población total R1
o = f2 + F2  # Hembras R2
k = m2 + M2  # Machos R2
j = m3 + M3  # Población R3

# Primer gráfico: Juveniles
plt.figure(figsize=(10,6))
plt.subplot(2,1,1)
plt.title("Tres Regiones Un Retardo")
plt.plot(t_span, f1, 'b-', lw=2, label='FR1 Juvenil')
plt.plot(t_span, f2, 'g--', lw=2, label='FR2 Juvenil')
plt.plot(t_span, m2, 'c:', lw=2, label='MR2 Juvenil')
plt.plot(t_span, m3, 'r-.', lw=2, label='MR3 Juvenil')
plt.ylabel("Población")
plt.legend()

# Segundo gráfico: Adultos y proporciones
plt.figure(figsize=(10,6))
plt.subplot(2,1,1)
plt.title("Dinámica Poblacional")
plt.plot(t_span, q, 'b-', lw=2, label='R1 Total')
plt.plot(t_span, o, 'g--', lw=2, label='FR2 Total')
plt.plot(t_span, k, 'c:', lw=2, label='MR2 Total')
plt.plot(t_span, j, 'r-.', lw=2, label='MR3 Total')
plt.ylabel("Población")
plt.legend()

plt.subplot(2,1,2)
plt.plot(t_span, r, 'm-', lw=2)
plt.title("Proporción de Machos")
plt.xlabel("Tiempo")
plt.ylabel("Ratio")

plt.tight_layout()
plt.show()