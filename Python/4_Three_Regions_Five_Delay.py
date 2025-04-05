import numpy as np
import matplotlib.pyplot as plt
from ddeint import ddeint

# 1. PARÁMETROS GLOBALES
d1 = 0.9
d2 = 0.151
d_new = 0.001
d3 = 0.139
b1 = 0.286
b2 = 0.844
k1 = 0.79*1.5
k2 = 0.15*1.5
k3 = 0.6*1.5
s1 = np.exp(-d1)
s8 = np.exp(-d2*7)
s11 = np.exp(-d2*3)
s12 = np.exp(-d_new)
s31 = np.exp(-d_new*19)
c1 = k2 / k1
c2 = k3 / k1
r1 = b1 / b2
a1 = d1 / b2
a2 = d2 / b2
a_new = d_new / b2
a3 = d3 / b2

taus = [1*b2, 8*b2, 11*b2, 12*b2, 31*b2]

# 2. FUNCIÓN DE HISTORIA
def history(t):
    return [
        # Valores iniciales (mismos que en MATLAB)
        5.0, 5.0, 5.0, 5.0,   # Af1, Af2, Am2, Am3
        5.0, 5.0, 5.0, 5.0,   # Bf1, Bf2, Bm2, Bm3
        5.0, 3.5, 2.0, 2.0,   # Cf1, Cf2, Cm2, Cm3
        5.0, 3.5, 2.0, 2.0,   # Df1, Df2, Dm2, Dm3
        2.5, 1.5, 1.0, 1.0,   # Ef1, Ef2, Em2, Em3
        2.5, 1.5, 1.0, 1.0    # Ff1, Ff2, Fm2, Fm3
    ]

# 3. ECUACIONES DDE
def model(Y, t):
    y = Y(t)
    # Iniciar el vector de derivadas
    dydt = np.zeros(24)
    
    # Obtener estados históricos
    Lag1 = Y(t - taus[0]) if t > taus[0] else history(t - taus[0])
    Lag2 = Y(t - taus[1]) if t > taus[1] else history(t - taus[1])
    Lag3 = Y(t - taus[2]) if t > taus[2] else history(t - taus[2])
    Lag4 = Y(t - taus[3]) if t > taus[3] else history(t - taus[3])
    Lag5 = Y(t - taus[4]) if t > taus[4] else history(t - taus[4])

    # Calcular variables auxiliares (Q1 a U6)
    Q1 = (y[8] + y[12] + y[16] + y[20]) / (1 + y[8] + y[12] + y[16] + y[20])
    Q2 = c1 / (c1 + y[8] + y[12] + y[16] + y[20] + y[9] + y[13] + y[17] + y[21])
    Q3 = (Lag1[8] + Lag1[12] + Lag1[16] + Lag1[20]) / (1 + Lag1[8] + Lag1[12] + Lag1[16] + Lag1[20])
    Q4 = c1 / (c1 + Lag1[8] + Lag1[12] + Lag1[16] + Lag1[20] + Lag1[9] + Lag1[13] + Lag1[17] + Lag1[21])
    Q5 = c2 / (c2 + y[8] + y[12] + y[16] + y[20] + y[9] + y[13] + y[17] + y[21])
    Q6 = c2 / (c2 + Lag1[8] + Lag1[12] + Lag1[16] + Lag1[20] + Lag1[9] + Lag1[13] + Lag1[17] + Lag1[21])
    
    R3 = (Lag2[8] + Lag2[12] + Lag2[16] + Lag2[20]) / (1 + Lag2[8] + Lag2[12] + Lag2[16] + Lag2[20])
    R4 = c1 / (c1 + Lag2[8] + Lag2[12] + Lag2[16] + Lag2[20] + Lag2[9] + Lag2[13] + Lag2[17] + Lag2[21])
    R6 = c2 / (c2 + Lag2[8] + Lag2[12] + Lag2[16] + Lag2[20] + Lag2[9] + Lag2[13] + Lag2[17] + Lag2[21])
    
    S3 = (Lag3[8] + Lag3[12] + Lag3[16] + Lag3[20]) / (1 + Lag3[8] + Lag3[12] + Lag3[16] + Lag3[20])
    S4 = c1 / (c1 + Lag3[8] + Lag3[12] + Lag3[16] + Lag3[20] + Lag3[9] + Lag3[13] + Lag3[17] + Lag3[21])
    S6 = c2 / (c2 + Lag3[8] + Lag3[12] + Lag3[16] + Lag3[20] + Lag3[9] + Lag3[13] + Lag3[17] + Lag3[21])
    
    T3 = (Lag4[8] + Lag4[12] + Lag4[16] + Lag4[20]) / (1 + Lag4[8] + Lag4[12] + Lag4[16] + Lag4[20])
    T4 = c1 / (c1 + Lag4[8] + Lag4[12] + Lag4[16] + Lag4[20] + Lag4[9] + Lag4[13] + Lag4[17] + Lag4[21])
    T6 = c2 / (c2 + Lag4[8] + Lag4[12] + Lag4[16] + Lag4[20] + Lag4[9] + Lag4[13] + Lag4[17] + Lag4[21])
    
    U3 = (Lag5[8] + Lag5[12] + Lag5[16] + Lag5[20]) / (1 + Lag5[8] + Lag5[12] + Lag5[16] + Lag5[20])
    U4 = c1 / (c1 + Lag5[8] + Lag5[12] + Lag5[16] + Lag5[20] + Lag5[9] + Lag5[13] + Lag5[17] + Lag5[21])
    U6 = c2 / (c2 + Lag5[8] + Lag5[12] + Lag5[16] + Lag5[20] + Lag5[9] + Lag5[13] + Lag5[17] + Lag5[21])
    
    # Población 0-1 años
    # Af1 Hembras 0-1 años R1
    dydt[0] = (r1*((y[8] + y[12])*(1 - Q1)) 
        - s1*r1*(Lag1[8] + Lag1[12])*(1 - Q3)
        - a1*y[0] 
        + (y[16] + y[20])*(1 - Q1) 
        - s1*(Lag1[16] + Lag1[20])*(1 - Q3))
    # Af2 Hembras 0-1 años R2
    dydt[1] = ((r1/2)*((y[9] + y[13] + (y[8] + y[12])*Q1)*Q2 
        - s1*(Lag1[9] + Lag1[13] + (Lag1[8] + Lag1[12])*Q3)*Q4)
        - a1*y[1] 
        + 0.5*((y[17] + y[21] + ((y[16] + y[20])*Q1)*Q2)
        - s1*(Lag1[17] + Lag1[21] + (Lag1[16] + Lag1[20])*Q3)*Q4))
    # Am2 Machos 0-1 años R2
    dydt[2] = ((r1/2)*((y[9] + y[13] + (y[8] + y[12])*Q1)*Q2 
        - s1*(Lag1[9] + Lag1[13] + (Lag1[8] + Lag1[12])*Q3)*Q4)
        - a1*y[2] 
        + 0.5*((y[17] + y[21] + ((y[16] + y[20])*Q1)*Q2)
        - s1*(Lag1[17] + Lag1[21] + (Lag1[16] + Lag1[20])*Q3)*Q4))
    # Am3 Machos 0-1 años R3
    dydt[3] = (r1*(Q5*((y[8] + y[12])*Q1 + y[9] + y[13])*(1 - Q2) 
        - s1*Q6*((Lag1[8] + Lag1[12])*Q3 + Lag1[9] + Lag1[13])*(1 - Q4)) 
        - a1*y[3] 
        + Q5*((y[16] + y[20])*Q1 + y[17] + y[21])*(1 - Q2) 
        - s1*Q6*((Lag1[16] + Lag1[20])*Q3 + Lag1[17] + Lag1[21])*(1 - Q4))
    
    # Poblacion 1-8 años
    # Bf1 Hembras 1-8 años R1
    dydt[4] = (s1*r1*(Lag1[8] + Lag1[12])*(1 - Q3) 
        + s1*(Lag1[16] + Lag1[20])*(1 - Q3) 
        - a2*y[4] 
        - s1*s8*r1*(Lag2[8] + Lag2[12])*(1 - R3) 
        - s1*s8*(Lag2[16] + Lag2[20])*(1 - R3))
    # Bf2 Hembras 1-8 años R2
    dydt[5] = (0.5*r1*s1*(Lag1[9] + Lag1[13] + (Lag1[8] + Lag1[12])*Q3)*Q4 
        + 0.5*s1*(Lag1[17] + Lag1[21] + (Lag1[16] + Lag1[20])*Q3)*Q4 
        - a2*y[5] 
        - 0.5*r1*s1*s8*(Lag2[9] + Lag2[13] + (Lag2[8] + Lag2[12])*R3)*R4 
        - 0.5*s1*s8*(Lag2[17] + Lag2[21] + (Lag2[16] + Lag2[20])*R3)*R4)
    # Bm2 Machos 1-8 años R2
    dydt[6] = (0.5*r1*s1*(Lag1[9] + Lag1[13] + (Lag1[8] + Lag1[12])*Q3)*Q4 
        + 0.5*s1*(Lag1[17] + Lag1[21] + (Lag1[16] + Lag1[20])*Q3)*Q4 
        - a2*y[6] 
        - 0.5*r1*s1*s8*(Lag2[9] + Lag2[13] + (Lag2[8] + Lag2[12])*R3)*R4 
        - 0.5*s1*s8*(Lag2[17] + Lag2[21] + (Lag2[16] + Lag2[20])*R3)*R4)
    # Bm3 Machos 1-8 años R3
    dydt[7] = (r1*s1*Q6*((Lag1[8] + Lag1[12])*Q3 + Lag1[9] + Lag1[13])*(1 - Q4) 
        + s1*Q6*((Lag1[16] + Lag1[20])*Q3 + Lag1[17] + Lag1[21])*(1 - Q4) 
        - a2*y[7] 
        - r1*s1*s8*R6*((Lag2[8] + Lag2[12])*R3 + Lag2[9] + Lag2[13])*(1 - R4) 
        - s1*s8*R6*((Lag2[16] + Lag2[20])*R3 + Lag2[17] + Lag2[21])*(1 - R4))
    
    # Poblacion 8-11 años
    # Cf1 Hembras 8-11 años R1  
    dydt[8] = (s1*s8*r1*(Lag2[8] + Lag2[12])*(1 - R3) 
        + s1*s8*(Lag2[16] + Lag2[20])*(1 - R3) 
        - a2*y[8] 
        - s1*s8*s11*r1*(Lag3[8] + Lag3[12])*(1 - S3) 
        - s1*s8*s11*(Lag3[16] + Lag3[20])*(1 - S3))
    # Cf2 Hembras 8-11 años R2
    dydt[9] = (0.5*r1*s1*s8*(Lag2[9] + Lag2[13] + (Lag2[8] + Lag2[12])*R3)*R4 
        + 0.5*s1*s8*(Lag2[17] + Lag2[21] + (Lag2[16] + Lag2[20])*R3)*R4 
        - a2*y[9] 
        - 0.5*r1*s1*s8*s11*(Lag3[9] + Lag3[13] + (Lag3[8] + Lag3[12])*S3)*S4 
        - 0.5*s1*s8*s11*(Lag3[17] + Lag3[21] + (Lag3[16] + Lag3[20])*S3)*S4)
    # Cm2 Machos 8-11 años R2
    dydt[10] = (0.5*r1*s1*s8*(Lag2[9] + Lag2[13] + (Lag2[8] + Lag2[12])*R3)*R4 
        + 0.5*s1*s8*(Lag2[17] + Lag2[21] + (Lag2[16] + Lag2[20])*R3)*R4 
        - a2*y[10] 
        - 0.5*r1*s1*s8*s11*(Lag3[9] + Lag3[13] + (Lag3[8] + Lag3[12])*S3)*S4 
        - 0.5*s1*s8*s11*(Lag3[17] + Lag3[21] + (Lag3[16] + Lag3[20])*S3)*S4)
    # Cm3 Machos 8-11 años R3
    dydt[11] = (r1*s1*s8*R6*((Lag2[8] + Lag2[12])*R3 + Lag2[9] + Lag2[13])*(1 - R4) 
        + s1*s8*R6*((Lag2[16] + Lag2[20])*R3 + Lag2[17] + Lag2[21])*(1 - R4) 
        - a2*y[11] 
        - r1*s1*s8*s11*S6*((Lag3[8] + Lag3[12])*S3 + Lag3[9] + Lag3[13])*(1 - S4) 
        - s1*s8*s11*S6*((Lag3[16] + Lag3[20])*S3 + Lag3[17] + Lag3[21])*(1 - S4))
    
    # Poblacion 11-12 años
    # Df1 Hembras 11-12 años R1
    dydt[12] = (s1*s8*s11*r1*(Lag3[8] + Lag3[12])*(1 - S3) 
        + s1*s8*s11*(Lag3[16] + Lag3[20])*(1 - S3) 
        - a_new*y[12] 
        - s1*s8*s11*s12*r1*(Lag4[8] + Lag4[12])*(1 - T3) 
        - s1*s8*s11*s12*(Lag4[16] + Lag4[20])*(1 - T3))
    # Df2 Hembras 11-12 años R2
    dydt[13] = (0.5*r1*s1*s8*s11*(Lag3[9] + Lag3[13] + (Lag3[8] + Lag3[12])*S3)*S4 
        + 0.5*s1*s8*s11*(Lag3[17] + Lag3[21] + (Lag3[16] + Lag3[20])*S3)*S4 
        - a_new*y[13] 
        - 0.5*r1*s1*s8*s11*s12*(Lag4[9] + Lag4[13] + (Lag4[8] + Lag4[12])*T3)*T4 
        - 0.5*s1*s8*s11*s12*(Lag4[17] + Lag4[21] + (Lag4[16] + Lag4[20])*T3)*T4)
    # Dm2 Machos 11-12 años R2
    dydt[14] = (0.5*r1*s1*s8*s11*(Lag3[9] + Lag3[13] + (Lag3[8] + Lag3[12])*S3)*S4 
        + 0.5*s1*s8*s11*(Lag3[17] + Lag3[21] + (Lag3[16] + Lag3[20])*S3)*S4 
        - a_new*y[14] 
        - 0.5*r1*s1*s8*s11*s12*(Lag4[9] + Lag4[13] + (Lag4[8] + Lag4[12])*T3)*T4 
        - 0.5*s1*s8*s11*s12*(Lag4[17] + Lag4[21] + (Lag4[16] + Lag4[20])*T3)*T4)
    # Dm3 Machos 11-12 años R3
    dydt[15] = (r1*s1*s8*s11*S6*((Lag3[8] + Lag3[12])*S3 + Lag3[9] + Lag3[13])*(1 - S4) 
        + s1*s8*s11*S6*((Lag3[16] + Lag3[20])*S3 + Lag3[17] + Lag3[21])*(1 - S4) 
        - a_new*y[15] 
        - r1*s1*s8*s11*s12*T6*((Lag4[8] + Lag4[12])*T3 + Lag4[9] + Lag4[13])*(1 - T4) 
        - s1*s8*s11*s12*T6*((Lag4[16] + Lag4[20])*T3 + Lag4[17] + Lag4[21])*(1 - T4))

    # Poblacion 12-31 años
    # Ef1 Hembras 12-31 años R1
    dydt[16] = (s1*s8*s11*s12*r1*(Lag4[8] + Lag4[12])*(1 - T3) 
        + s1*s8*s11*s12*(Lag4[16] + Lag4[20])*(1 - T3) 
        - a_new*y[16] 
        - s1*s8*s11*s12*s31*r1*(Lag5[8] + Lag5[12])*(1 - U3) 
        - s1*s8*s11*s12*s31*(Lag5[16] + Lag5[20])*(1 - U3))
    # Ef2 Hembras 12-31 años R2
    dydt[17] = (0.5*r1*s1*s8*s11*s12*(Lag4[9] + Lag4[13] + (Lag4[8] + Lag4[12])*T3)*T4 
        + 0.5*s1*s8*s11*s12*(Lag4[17] + Lag4[21] + (Lag4[16] + Lag4[20])*T3)*T4 
        - a_new*y[17] 
        - 0.5*r1*s1*s8*s11*s12*s31*(Lag5[9] + Lag5[13] + (Lag5[8] + Lag5[12])*U3)*U4 
        - 0.5*s1*s8*s11*s12*s31*(Lag5[17] + Lag5[21] + (Lag5[16] + Lag5[20])*U3)*U4)
    # Em2 Machos 12-31 años R2
    dydt[18] = (0.5*r1*s1*s8*s11*s12*(Lag4[9] + Lag4[13] + (Lag4[8] + Lag4[12])*T3)*T4 
        + 0.5*s1*s8*s11*s12*(Lag4[17] + Lag4[21] + (Lag4[16] + Lag4[20])*T3)*T4 
        - a_new*y[18] 
        - 0.5*r1*s1*s8*s11*s12*s31*(Lag5[9] + Lag5[13] + (Lag5[8] + Lag5[12])*U3)*U4 
        - 0.5*s1*s8*s11*s12*s31*(Lag5[17] + Lag5[21] + (Lag5[16] + Lag5[20])*U3)*U4)
    # Em3 Machos 12-31 años R3
    dydt[19] = (r1*s1*s8*s11*s12*T6*((Lag4[8] + Lag4[12])*T3 + Lag4[9] + Lag4[13])*(1 - T4) 
        + s1*s8*s11*s12*T6*((Lag4[16] + Lag4[20])*T3 + Lag4[17] + Lag4[21])*(1 - T4) 
        - a_new*y[19] 
        - r1*s1*s8*s11*s12*s31*U6*((Lag5[8] + Lag5[12])*U3 + Lag5[9] + Lag5[13])*(1 - U4) 
        - s1*s8*s11*s12*s31*U6*((Lag5[16] + Lag5[20])*U3 + Lag5[17] + Lag5[21])*(1 - U4))

    #Poblacion 31+ años
    # Ff1 Hembras 31+ años R1
    dydt[20] = (s1*s8*s11*s12*s31*r1*(Lag5[8] + Lag5[12])*(1/(1 + Lag5[8] + Lag5[12] + Lag5[16] + Lag5[20])) 
        + s1*s8*s11*s12*s31*(Lag5[16] + Lag5[20])*(1/(1 + Lag5[8] + Lag5[12] + Lag5[16] + Lag5[20])) 
        - a3*y[20])
    # Ff2 Hembras 31+ años R2
    dydt[21] = (0.5*r1*s1*s8*s11*s12*s31*(Lag5[9] + Lag5[13] + (Lag5[8] + Lag5[12])*U3)*U4 
        + 0.5*s1*s8*s11*s12*s31*(Lag5[17] + Lag5[21] + (Lag5[16] + Lag5[20])*U3)*U4 
        - a3*y[21])
    # Fm2 Machos 31+ años R2
    dydt[22] = (0.5*r1*s1*s8*s11*s12*s31*(Lag5[9] + Lag5[13] + (Lag5[8] + Lag5[12])*U3)*U4 
        + 0.5*s1*s8*s11*s12*s31*(Lag5[17] + Lag5[21] + (Lag5[16] + Lag5[20])*U3)*U4 
        - a3*y[22])
    # Fm3 Machos 31+ años R3
    dydt[23] = (r1*s1*s8*s11*s12*s31*U6*((Lag5[8] + Lag5[12])*U3 + Lag5[9] + Lag5[13])*(1 - U4) 
        + s1*s8*s11*s12*s31*U6*((Lag5[16] + Lag5[20])*U3 + Lag5[17] + Lag5[21])*(1 - U4) 
        - a3*y[23])

    return dydt

# 4. SOLUCIÓN NUMÉRICA
t = np.linspace(0, 200, 1000)
sol = ddeint(model, history, t)

# Extraer las variables de la solución

# Poblaciones por edad y región (24 variables)
Af1 = sol[:, 0]   # Hembras 0-1 años R1
Af2 = sol[:, 1]   # Hembras 0-1 años R2
Am2 = sol[:, 2]   # Machos 0-1 años R2
Am3 = sol[:, 3]   # Machos 0-1 años R3

Bf1 = sol[:, 4]   # Hembras 1-8 años R1
Bf2 = sol[:, 5]   # Hembras 1-8 años R2
Bm2 = sol[:, 6]   # Machos 1-8 años R2
Bm3 = sol[:, 7]   # Machos 1-8 años R3

Cf1 = sol[:, 8]   # Hembras 8-11 años R1
Cf2 = sol[:, 9]   # Hembras 8-11 años R2
Cm2 = sol[:, 10]  # Machos 8-11 años R2
Cm3 = sol[:, 11]  # Machos 8-11 años R3

Df1 = sol[:, 12]  # Hembras 11-12 años R1
Df2 = sol[:, 13]  # Hembras 11-12 años R2
Dm2 = sol[:, 14]  # Machos 11-12 años R2
Dm3 = sol[:, 15]  # Machos 11-12 años R3

Ef1 = sol[:, 16]  # Hembras 12-31 años R1
Ef2 = sol[:, 17]  # Hembras 12-31 años R2
Em2 = sol[:, 18]  # Machos 12-31 años R2
Em3 = sol[:, 19]  # Machos 12-31 años R3

Ff1 = sol[:, 20]  # Hembras 31+ años R1
Ff2 = sol[:, 21]  # Hembras 31+ años R2
Fm2 = sol[:, 22]  # Machos 31+ años R2
Fm3 = sol[:, 23]  # Machos 31+ años R3

# Variables compuestas
u = Af1 + Bf1 + Cf1 + Df1 + Ef1 + Ff1  # Total hembras R1
v = Af2 + Bf2 + Cf2 + Df2 + Ef2 + Ff2  # Total hembras R2
w = Am2 + Bm2 + Cm2 + Dm2 + Em2 + Fm2  # Total machos R2
x = Am3 + Bm3 + Cm3 + Dm3 + Em3 + Fm3  # Total machos R3

y_total = u + v  # Total población femenina
z_total = w + x  # Total población masculina

r = z_total / (y_total + z_total)  # Proporción masculina

# 5. GRÁFICOS

# Figura 1: Población 0-12 años
plt.figure(figsize=(12, 6))
plt.title("Población de 0-12 años")
plt.grid(True)

# Grupo 0-1 años
plt.plot(t, Af1, 'g--', linewidth=3, label='Af1 (0-1 R1)')
plt.plot(t, Af2, 'r--', linewidth=3, label='Af2 (0-1 R2)')
plt.plot(t, Am2, 'y--', linewidth=3, label='Am2 (0-1 R2)')
plt.plot(t, Am3, 'b--', linewidth=3, label='Am3 (0-1 R3)')

# Grupo 1-8 años
plt.plot(t, Bf1, 'g:', linewidth=3, label='Bf1 (1-8 R1)')
plt.plot(t, Bf2, 'r:', linewidth=3, label='Bf2 (1-8 R2)')
plt.plot(t, Bm2, 'y:', linewidth=3, label='Bm2 (1-8 R2)')
plt.plot(t, Bm3, 'b:', linewidth=3, label='Bm3 (1-8 R3)')

# Grupo 8-11 años
plt.plot(t, Cf1, 'g-.', linewidth=3, label='Cf1 (8-11 R1)')
plt.plot(t, Cf2, 'r-.', linewidth=3, label='Cf2 (8-11 R2)')
plt.plot(t, Cm2, 'y-.', linewidth=3, label='Cm2 (8-11 R2)')
plt.plot(t, Cm3, 'b-.', linewidth=3, label='Cm3 (8-11 R3)')

# Grupo 11-12 años
plt.plot(t, Df1, 'g-', linewidth=3, label='Df1 (11-12 R1)')
plt.plot(t, Df2, 'r-', linewidth=3, label='Df2 (11-12 R2)')
plt.plot(t, Dm2, 'y-', linewidth=3, label='Dm2 (11-12 R2)')
plt.plot(t, Dm3, 'b-', linewidth=3, label='Dm3 (11-12 R3)')

plt.xlabel('Tiempo (años)')
plt.ylabel('Población')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Figura 2: Población 12+ años
plt.figure(figsize=(12, 6))
plt.title("Población de 12+ años")
plt.grid(True)

# Grupo 12-31 años
plt.plot(t, Ef1, 'g--', linewidth=3, label='Ef1 (12-31 R1)')
plt.plot(t, Ef2, 'r--', linewidth=3, label='Ef2 (12-31 R2)')
plt.plot(t, Em2, 'y--', linewidth=3, label='Em2 (12-31 R2)')
plt.plot(t, Em3, 'b--', linewidth=3, label='Em3 (12-31 R3)')

# Grupo 31+ años
plt.plot(t, Ff1, 'g-', linewidth=3, label='Ff1 (31+ R1)')
plt.plot(t, Ff2, 'r-', linewidth=3, label='Ff2 (31+ R2)')
plt.plot(t, Fm2, 'y-', linewidth=3, label='Fm2 (31+ R2)')
plt.plot(t, Fm3, 'b-', linewidth=3, label='Fm3 (31+ R3)')

plt.xlabel('Tiempo (años)')
plt.ylabel('Población')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Figura 3: Ratio de machos
plt.figure(figsize=(12, 6))
plt.plot(t, r, 'm-', linewidth=3)
plt.title("Proporción de machos en la población total")
plt.xlabel('Tiempo (años)')
plt.ylabel('Proporción')
plt.grid(True)
plt.tight_layout()

# Mostrar todos los gráficos
plt.show()