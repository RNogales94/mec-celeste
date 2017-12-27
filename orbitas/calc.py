import numpy as np
import matplotlib.pyplot as plt

planets = np.array(["Mercurio","Venus","Tierra", "Marte","Júpiter","Saturno","Urano","Neptuno"])
epsilon = np.array([0.206,0.007,0.017,0.093,0.048,0.056,0.047,0.009])
a = np.array([0.387,0.723,1,1.524,5.203,9.546,19.2,30.09])
period = np.array([87.97,224.7,365.26,686.98,4332.6,10759,30687,60784])
masa = np.array([3.302*10**23,4.869*10**24,5.9742*10**24,6.4191*10**23,1.8987*10**27,5.6851*10**26,8.6849*10**25,1.0244*10**26])


# Método de Newton de Raphson
import math
import random
PI = math.pi

def newton_raphson(a, epsilon, period, t, t_ini, tolerance):
    ji = 2*PI*(t - t_ini)/period
    def phi(u):
        return (epsilon*(math.sin(u) - u*math.cos(u))+ji)/(1 - epsilon*math.cos(u))
    # Tomamos u_0 = PI
    u = PI
    phi_u = phi(u)

    while(abs(phi_u - u) > tolerance):
        u = phi_u
        phi_u = phi(u)
    return u

# Método de Newton de Bessel
import scipy.special as sp #import de las funciones especiales (como la de Bessel)

def bessel_method(a, epsilon, period, t, t_ini, tolerance):
    ji = 2*PI*(t-t_ini)/period
    def series_term(n):
        #sp.jv(n, z) es la funcion de bessel de orden n y argumento z
        return 2/n * sp.jv(n, n*epsilon)*math.sin(n*ji)

    u = ji + series_term(1)
    bessel_u = u + series_term(2)
    n = 3
    while(abs(bessel_u - u) > tolerance):
        u = bessel_u
        bessel_u = bessel_u + series_term(n)
        n = n+1

    return bessel_u

def position(a, e, u):
    x = a*(math.cos(u) - e)
    y = a*(math.sqrt(1 - e**2)*math.sin(u))
    return [x, y]

  # a*c(-sin(result.nr) * u.diff, sqrt(1-epsilon^2) * cos(result.nr) * u.diff, 0)

def velocidad(a, e, p, u):
    """
    Calculamos la velocidad como la derivada de la posicion
    Para lo cual vemos la posición como una función de u pos(u)
    Es por esto que necesitamos calcular la derivada de u: du
    """
    du =  2*PI/(p*(1-e * math.cos(u)))
    x = a*(-math.sin(u)*du)
    y = a*(math.sqrt(1 - e**2) * math.cos(u) * du)
    return [x,y]

def excent2real(u,e):
    return 2*math.atan(math.sqrt((1-e)/(1+e))*math.tan(u/2))

def energy(a, p):
    mu = a**3*4*PI**2 / p**2
    return -mu/(2*a)

def mu_2Cuerpos(masa_planeta):
    G = 6.67384*10**(-11)     #Constante de gravitacion universal
    K = 86400**2 / 149597870700**3 # Cambio de unidades: s --> dias , m --> ua
    G = K*G
    M = 1.98892*10**(30)    #Masa del Sol
    mu = G * M**3 / ((M + masa_planeta)**2)
    return mu

# calculo energia total
def energy2 (a, p, tol, e, t, ti, dos_cuerpos=False, planet=3):
    # (1/2)M vel^2
    vel = velocidad(a, e, p, newton_raphson(a, e, p, t, ti, tol))
    x = position(a, e, newton_raphson(a, e, p, t, ti, tol))

    if problema_2_cuerpos:
        masa_planeta = masa[planet]
        mu = mu_2Cuerpos(masa_planeta)
    else:
        mu = a**3*4*PI**2 / p**2

    #  energia = sum(velocidad^2)/2 - mu/sqrt( sum(posicion.nr^2) ),
    sum_vel2 = vel[0]**2 + vel[1]**2
    sum_pos2 = x[0]**2 + x[1]**2
    E_cin = sum_vel2/2
    E_pot = mu/math.sqrt(sum_pos2)
    E_total = E_cin - E_pot
    return E_total

"""
def momento(t, p, a, e, mu, u):
    #myFunc = lambda x: x * 5
    u = lambda t: newton_raphson(a, e, p, t, 0, 0.00000000001)
    du = lambda t: math.sqrt(mu) / a**(3/2)*(1-e*math.cos(u(t)))
    c = a**2 * math.sqrt(1-e**2) * du(t) * (1-e*cos(u(t)))
    return c
"""


if __name__ == "__main__":
    t_0 = 0
    tolerance = 0.00000000001
    problema_2_cuerpos = True
    while(True):
        planet = int(input("Elija un planeta: ")) -1
        t = int(input("Dias desde el pericentro: "))

        anomExc_nr = newton_raphson(a[planet], epsilon[planet], period[planet], t, t_0, tolerance)
        anomExc_be = bessel_method(a[planet], epsilon[planet], period[planet], t, t_0, tolerance)
        print("Anomalía excentrica (NR): ", anomExc_nr)
        print("Anomalía excentrica (Bessel): ", anomExc_be)
        anomReal = excent2real(anomExc_nr, epsilon[planet])
        print("Anomalía real: ", anomReal)
        if problema_2_cuerpos:
            energ2 = energy2(a[planet], period[planet], tolerance, epsilon[planet], t, t_0, True, planet)
            print("Energia2: ", energ2)
        else:
            energ = energy(a[planet], period[planet])
            print("Energia: ", energ)
            energ2 = energy2(a[planet], period[planet], tolerance, epsilon[planet], t, t_0)
            print("Energia2: ", energ2)
        pos = position(a[planet], epsilon[planet], anomExc_nr)
        print("Posicion NR: ", pos)
        pos = position(a[planet], epsilon[planet], anomExc_be)
        print("Posicion Be: ", pos)
        vel = velocidad(a[planet], epsilon[planet], period[planet], anomExc_be)
        print("Velocidad: ", vel)


        u = lambda t: newton_raphson(a[planet], epsilon[planet], period[planet], t, t_0, tolerance)
        pos = lambda t: position(a[planet], epsilon[planet], u(t))
        X = [pos(t)[0] for t in range(1000)]
        Y = [pos(t)[1] for t in range(1000)]
        #plt.scatter(X,Y)
        #plt.show()
