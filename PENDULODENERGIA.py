

from numpy import array, linspace
from math import sin, cos, pi
from pylab import plot, xlabel, ylabel, show
from scipy.integrate import odeint

from vpython import sphere, scene, vector, color, arrow, text, sleep,cylinder

arrow_size = 0.1

arrow_x = arrow(pos=vector(0,0,0), axis=vector(arrow_size,0,0), color=color.red)
arrow_y = arrow(pos=vector(0,0,0), axis=vector(0,arrow_size,0), color=color.green)
arrow_z = arrow(pos=vector(0,0,0), axis=vector(0,0,arrow_size))

R = 0.03 #Radio de la esfera

def func (r, t, g, l): 
    theta1 = r[0]
    omega1 = r[1]
    theta2 = r[2]
    omega2 = r[3]
    
    f_omega1 = - (omega1 ** 2 * sin(2 * theta1 - 2 * theta2) + 2 * omega2 ** 2 * sin(theta1 - theta2) + \
                  g / l * (sin(theta1 - 2 * theta2) + 3 * sin(theta1))) / (3 - cos(2 * theta1 - 2 * theta2))
    f_omega2 = (4 * omega1 ** 2 * sin(theta1 - theta2) + omega2 ** 2 * sin(2 * theta1 - 2 * theta2) + \
                 2 * g / l * (sin(2 * theta1 - theta2) - sin(theta2))) / (3 - cos(2 * theta1 - 2 * theta2))
    
    return array([omega1,f_omega1,omega2,f_omega2], float)


##############################################
def energia(r):
    theta1 = r[0]
    omega1 = r[1]
    theta2 = r[2]
    omega2 = r[3]
    ener=  - m * g * l * (2 * cos(theta1) + cos(theta2)) +  m * l ** 2 * (omega1 ** 2 + 0.5 * omega2 ** 2 + omega1 * omega2 * cos(theta1 - theta2))
    return array([ener],float)

##############################################

##################################
# Constants
g = 9.81  # m/s^2
m = 1. # kg
l = 0.4  # pendulum lengths in m
theta1_0 = pi/2.
theta2_0 = pi/2.
omega1_0 = 0.
omega2_0 = 0.

###################################
n_steps = 20000 
ti = 0.  
tf = 10.  
t_delta = (tf - ti) / n_steps 
t = linspace(ti, tf, n_steps) 

##################################

#initcond = array([phi1,a,phi2,b,da,db])
r = array([theta1_0, omega1_0, theta2_0, omega2_0], float)

solu, outodeint = odeint( func,r, t, args =(g,l),full_output=True) 
theta10, omega10, theta20, omega20 = solu.T 
#########################################################################################################################
# probemos esta energia para ver si sale 

ener=[]

for i in range(len())
    ener.append (- m * g * l * (2 * cos(theta10[i]) + cos(theta20[i])) +  m * l ** 2 * (omega10[i] ** 2 + 0.5 * omega20[i] ** 2 + omega10[i] * omega20[i] * cos(theta10[i] - theta20[i])) )



plot(t, ener)
xlabel('t (s)')
ylabel('energia (J)')
show()


# =====================

scene.range = 0.5 

xp = l*sin(theta1_0) #Pasa de coordenadas polares a cartesianas
yp = -l*cos(theta1_0)
zp = 0.

xs=l*(sin(theta1_0)+sin(theta2_0))
ys=-l*(cos(theta1_0)+cos(theta2_0))
zs=0.

sleeptime = 0.0001 #Tiempo con que se actualiza la posición de la partícula


esfera1=prtcl = sphere(pos=vector(xp,yp,zp), radius=l/10, color=color.yellow) #Define objeto con que se va a trabajar
esfera2=prtcls= sphere(pos=vector(xs,ys,zs), radius=l/10, color=color.blue)

cilindro1 = cylinder(pos=vector(0, 0, 0), axis=vector(xp, yp, 0), radius=l/40)
cilindro2 = cylinder(pos=vector(xp, yp, 0), axis=vector(xp, yp, 0), radius=l/40)
               
time_i = 0 #Contador que se mueve en el espacio temporal en el que se resolvió la ecuación diferencial


#for i in omega:
#    print(i)


while ti < tf: #ANIMACIÓN
    vector1 = vector( l*sin(theta10[time_i]), -l*cos(theta10[time_i]), zp )
    vector2= vector( l*(sin(theta10[time_i])+sin(theta20[time_i])), -l*(cos(theta10[time_i])+cos(theta20[time_i])), zs )
    
    cilindro1.axis = vector1
    esfera1.pos =  vector1
    cilindro2.pos =  vector1
    cilindro2.axis = vector2
    esfera2.pos = vector1 + vector2
    
    ti += t_delta
    sleep(sleeptime)
    time_i += 1




    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
