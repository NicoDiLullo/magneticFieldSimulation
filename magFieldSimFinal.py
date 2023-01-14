# Neutron star accretion basic simulation
import numpy as np
import matplotlib.pyplot as plt

# simulation variables; set to same as Dr. Mushtukov's code to test
time = 0.0 #time, seconds
step_time = 1.e-5
# can run sim by steps as opposed to by lambda
tot_time = 1.
steps = tot_time/step_time

# particle vars
initialVelocity = 1. #velocity m/s
l_deg = 0. #lamba; particle angle relative to star
r_m = 1.
velocity = initialVelocity

# ns info
mass = 1.4 # mass, kg, the neutron star
spinPeriod = .7 # spin period
omega = 2*np.pi/spinPeriod
RNS = 1.3 # radius, ^6
# old functions not used right now/combined into other methods.
'''
def field_lines(lambda_val, rsubm): #lambda in radians <-?
    lines = rsubm*np.power(np.cos(lambda_val), 2)
    return lines

def angle_grav_mag(lambda_val): #lambda still in radians <-?
    return np.arctan(0.5/np.tan(lambda_val))
'''
# acceleration due to gravitational force on a particle in a magnetic field line
def acceleration_gravity(l_deg, mass, r_m):
    x = np.arctan((0.5/np.arctan(l_deg))) # velo* grav force: chi = arctan[.5/tan(lambda)]
    toReturn = (1.328e5*mass/r_m**2/(np.cos(l_deg))**4 * np.cos(x))
    return toReturn

# acceleration due to centrifugal force on a particle in a magnetic field line
def centri_acceleration(l_deg, mass, r_m, omega):
    x = np.arctan((0.5/np.arctan(l_deg))) # same as above
    radiusFieldLine = r_m*(np.cos(l_deg)**2) # dipole magnetic field lines. r(lambda) = r_sub_m*cos^2(lambda)
    b = np.sin(x-l_deg)*np.sin(l_deg)*1.328e5*mass/radiusFieldLine**2
    toReturn = (np.cos(x-l_deg)*np.cos(l_deg)* ( 1.328e5*mass/radiusFieldLine**2 - 1.e3*omega**2*radiusFieldLine ) - b)
    return toReturn

det = 11
radiusFieldLine = r_m

#velocity negative search
#bissection search
#0-inner disc radius
#be sure code works for any initial coordinates and velocities
def bissection_search(mass):
    upper = 
    lower = 
    mid = (upper + lower)

#while time<5.e0: #better loop condition?
while radiusFieldLine > 1.e-2 and velocity>0: #stop code as soon as velocity < zero
#while time<tot_time: #better loop condition?
    # step
    #print("\n beta_1",velocity/3.e5)
    x_coord = step_time*velocity # distance = velocity*time
    deltaLambda = (x_coord/r_m*1.e-3)/np.cos(l_deg)/np.sqrt(1.+3.*(np.sin(l_deg))**2) # change in l_deg
    l_deg = l_deg + deltaLambda
    radiusFieldLine = r_m*(np.cos(l_deg)**2)
    gA = acceleration_gravity(l_deg, mass, r_m)
    cA = centri_acceleration(l_deg, mass, r_m, omega)
    velocity = velocity + cA*step_time
    #print("x_coord",x_coord)
    #print("deltaLambda",deltaLambda)
    #print("l_deg",l_deg)
    #print("radiusFieldLine",radiusFieldLine/1.e-2)
    #print("velocity_2",velocity)
    # account
    #step_time = 0.05/velocity
    time = time + step_time
    # plot
    #plt.new()
    #plt.figure(x_coord, lambda)
    #y_coord = x_coord*np.tan(l_deg)
    #plt.plot(r_m-x_coord, y_coord, 'o')
    plt.plot(velocity, time, 'o')
plt.show()
print(velocity)
    
#print("velocity")
#print(velocity)
#print("x_coord")
#print(x coord)
#print("lambda")
#ßprint(l_deg)


