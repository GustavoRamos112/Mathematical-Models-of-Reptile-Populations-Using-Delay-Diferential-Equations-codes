import numpy as np
import matplotlib.pyplot as plt
from ddeint import ddeint

# 1. Parámetros globales
# Tiempos inicial y final
t0 = 0
tF = 200

# Tasa de natalidad adulta
b2 = 0.844

# Cálculo de retardos
tau_2 = 8 * b2
tau_3 = 11 * b2
lags = [tau_2, tau_3]

# Parámetros del modelo
d1 = 0.14     # Tasa de mortalidad juveniles
d2 = 0.11     # Tasa de mortalidad adolescentes
d3 = 0.0695   # Tasa de mortalidad adultos
b1 = 0.286    # Tasa de natalidad adolescentes
k1 = 0.79     # Capacidad de carga R1
k2 = 0.15     # Capacidad de carga R2
k3 = 0.06     # Capacidad de carga R3

# Cálculos derivados
s8 = np.exp(-d1 * 8)     # Supervivencia a 8 años
s11 = np.exp(-d2 * 3)    # Supervivencia a 11 años
c1 = k2 / k1             # Ratio capacidad R2/R1
c2 = k3 / k1             # Ratio capacidad R3/R1
r1 = b1 / b2             # Ratio natalidad adolescentes/adultos
a1 = d1 / b2             # Mortalidad juvenil normalizada
a2 = d2 / b2             # Mortalidad adolescente normalizada
a3 = d3 / b2             # Mortalidad adulta normalizada

# 2. Función de historia (valores iniciales)
def historia(t):
    """Función que proporciona valores históricos constantes para t <= 0"""
    return [10, 10, 10, 10, 10, 7, 4, 4, 5, 3, 1, 1]

# 3. Sistema de ecuaciones diferenciales con retardos
def modelo(Y, t):
    # Función para acceder a estados pasados
    def get_delayed(var_index, delay):
        return Y(t - delay)[var_index] if t - delay > t0 else historia(t - delay)[var_index]

    # Variables actuales
    y = Y(t)

    # Variables retardadas
    y_tau2 = [get_delayed(i, tau_2) for i in range(12)]  # Retardo tau_2
    y_tau3 = [get_delayed(i, tau_3) for i in range(12)]  # Retardo tau_3

    # Calcular factores Q y R
    Q1 = (y[4] + y[8]) / (1 + y[4] + y[8])
    Q2 = c1 / (c1 + y[4] + y[5] + y[8] + y[9])
    Q3 = (y_tau2[4] + y_tau2[8]) / (1 + y_tau2[4] + y_tau2[8])
    Q4 = c1 / (c1 + y_tau2[4] + y_tau2[5] + y_tau2[8] + y_tau2[9])
    Q5 = c2 / (c2 + y[4] + y[5] + y[8] + y[9])
    Q6 = c2 / (c2 + y_tau2[4] + y_tau2[5] + y_tau2[8] + y_tau2[9])
    R3 = (y_tau3[4] + y_tau3[8]) / (1 + y_tau3[4] + y_tau3[8])
    R4 = c1 / (c1 + y_tau3[4] + y_tau3[5] + y_tau3[8] + y_tau3[9])
    R6 = c2 / (c2 + y_tau3[4] + y_tau3[5] + y_tau3[8] + y_tau3[9])

    # Inicializar vector de derivadas
    dydt = np.zeros(12)

    # Población juvenil
    # Af1: Juveniles hembras R1
    dydt[0] = r1*(((y[4])*(1-Q1)) - s8*r1*y_tau2[4]*(1-Q3)) - a1*y[0] + y[8]*(1-Q1) - s8*y_tau2[8]*(1-Q3)

    # Af2: Juveniles hembras R2
    dydt[1] = (r1/2)*(((y[4]+y[5])*Q1)*Q2 - s8*(y_tau2[4]+y_tau2[5])*Q3*Q4) - a1*y[1] + 0.5*((y[8]+y[9])*Q1)*Q2 - s8*(y_tau2[8]+y_tau2[9])*Q3*Q4

    # Am2: Juveniles machos R2
    dydt[2] = (r1/2)*((y[4]+y[5])*Q1*Q2 - s8*(y_tau2[4]+y_tau2[5])*Q3*Q4) - a1*y[2] + 0.5*((y[8]+y[9])*Q1)*Q2 - s8*(y_tau2[8]+y_tau2[9])*Q3*Q4 

    # Am3: Juveniles machos R3
    dydt[3] = r1*(Q5*(y[4]*Q1 + y[5])*(1-Q2) - s8*Q6*(y_tau2[4]*Q3 + y_tau2[5])*(1-Q4)) - a1*y[3] + Q5*(y[8]*Q1 + y[9])*(1-Q2) - s8*Q6*(y_tau2[8]*Q3 + y_tau2[9])*(1-Q4)

    # Población adolescente
    # Bf1: Adolescentes hembras R1
    dydt[4] = s8*r1*y_tau2[4]*(1-Q3) + s8*y_tau2[8]*(1-Q3) - a2*y[4] - s8*s11*r1*y_tau3[4]*(1-R3) - s8*s11*y_tau3[8]*(1-R3)

    # Bf2: Adolescentes hembras R2
    dydt[5] = s8*(r1/2)*(y_tau2[4]+y_tau2[5])*Q3*Q4 + s8*(y_tau2[8]+y_tau2[9])*Q3*Q4 - a2*y[5] - s8*s11*(r1/2)*(y_tau3[4]+y_tau3[5])*R3*R4 - s8*s11*(y_tau3[8]+y_tau3[9])*R3*R4

    # Bm2: Adolescentes machos R2
    dydt[6] = s8*(r1/2)*(y_tau2[4]+y_tau2[5])*Q3*Q4 + s8*(y_tau2[8]+y_tau2[9])*Q3*Q4 - a2*y[6] - s8*s11*(r1/2)*(y_tau3[4]+y_tau3[5])*R3*R4 - s8*s11*(y_tau3[8]+y_tau3[9])*R3*R4

    # Bm3: Adolescentes machos R3
    dydt[7] = s8*Q6*r1*(y_tau2[4]*Q3 + y_tau2[5])*(1-Q4) + s8*Q6*(y_tau2[8]*Q3 + y_tau2[9])*(1-Q4) - a2*y[7] - s8*R6*r1*(y_tau3[4]*R3 + y_tau3[5])*(1-R4) - s8*R6*(y_tau3[8]*R3 + y_tau3[9])*(1-R4)

    # Población adulta
    # Cf1: Adultas R1
    dydt[8] = s8*s11*r1*y_tau3[4]*(1-R3) + s8*s11*y_tau3[8]*(1-R3) - a3*y[8]

    # Cf2: Adultas R2
    dydt[9] = s8*s11*(r1/2)*(y_tau3[4]+y_tau3[5])*R3*R4 + s8*s11*(y_tau3[8]+y_tau3[9])*R3*R4 - a3*y[9]

    # Cm2: Adultos R2
    dydt[10] = s8*s11*(r1/2)*(y_tau3[4]+y_tau3[5])*R3*R4 + s8*s11*(y_tau3[8]+y_tau3[9])*R3*R4 - a3*y[10]

    # Cm3: Adultos R3
    dydt[11] = s8*R6*r1*(y_tau3[4]*R3 + y_tau3[5])*(1-R4) + s8*R6*(y_tau3[8]*R3 + y_tau3[9])*(1-R4) - a3*y[11]

    return dydt

# 4. Resolución del sistema

# Rango temporal
tt = np.linspace(t0, tF, 1000)

# Resolver ecuaciones diferenciales
sol = ddeint(modelo, historia, tt)

# Extraer variables de la solución
Af1 = sol[:, 0]
Af2 = sol[:, 1]
Am2 = sol[:, 2]
Am3 = sol[:, 3]
Bf1 = sol[:, 4]
Bf2 = sol[:, 5] 
Bm2 = sol[:, 6]
Bm3 = sol[:, 7]
Cf1 = sol[:, 8]
Cf2 = sol[:, 9]
Cm2 = sol[:, 10] 
Cm3 = sol[:, 11]

u = Af1 + Bf1 + Cf1  # Total hembras R1
v = Af2 + Bf2 + Cf2  # Total hembras R2
w = Am2 + Bm2 + Cm2  # Total machos R2
q = Am3 + Bm3 + Cm3  # Total machos R3
y_total = u + v      # Total hembras
z_total = w + q      # Total machos
r = z_total / (y_total + z_total)  # Proporción de machos

# 5. Graficación
# Figura 1: Poblaciones por categoría
plt.figure(figsize=(12, 6))
plt.plot(tt, Af1, 'g-', linewidth=3, label='Af1 (Juv. Hembras R1)')
plt.plot(tt, Af2, 'r-', linewidth=3, label='Af2 (Juv. Hembras R2)')
plt.plot(tt, Am2, 'y-', linewidth=3, label='Am2 (Juv. Machos R2)')
plt.plot(tt, Am3, 'b-', linewidth=3, label='Am3 (Juv. Machos R3)')
plt.plot(tt, Bf1, 'g:', linewidth=3, label='Bf1 (Adol. Hembras R1)')
plt.plot(tt, Bf2, 'r:', linewidth=3, label='Bf2 (Adol. Hembras R2)')
plt.plot(tt, Bm2, 'y:', linewidth=3, label='Bm2 (Adol. Machos R2)')
plt.plot(tt, Bm3, 'b:', linewidth=3, label='Bm3 (Adol. Machos R3)')
plt.plot(tt, Cf1, 'g--', linewidth=3, label='Cf1 (Adultas R1)')
plt.plot(tt, Cf2, 'r--', linewidth=3, label='Cf2 (Adultas R2)')
plt.plot(tt, Cm2, 'y--', linewidth=3, label='Cm2 (Adultos R2)')
plt.plot(tt, Cm3, 'b--', linewidth=3, label='Cm3 (Adultos R3)')
plt.xlabel('Tiempo')
plt.ylabel('Población')
plt.title('Dinámica poblacional de caimanes')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Figura 2: Poblaciones totales por región
plt.figure(figsize=(12, 6))
plt.plot(tt, u, 'g-', linewidth=3, label='Hembras R1')
plt.plot(tt, v, 'r-', linewidth=3, label='Hembras R2')
plt.plot(tt, w, 'y-', linewidth=3, label='Machos R2')
plt.plot(tt, q, 'b-', linewidth=3, label='Machos R3')
plt.xlabel('Tiempo')
plt.ylabel('Población total')
plt.title('Poblaciones totales por región y sexo')
plt.legend()
plt.tight_layout()

# Figura 3: Proporción de machos
plt.figure(figsize=(12, 6))
plt.plot(tt, r, 'r-', linewidth=3)
plt.xlabel('Tiempo')
plt.ylabel('Proporción de machos')
plt.title('Relación machos/población total ("Gráfico del amor")')
plt.tight_layout()

plt.show()