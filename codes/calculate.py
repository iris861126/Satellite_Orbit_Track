#=====================================
# calculating code
# save data: lats, longs
#=====================================
import numpy as np
from datetime import datetime

while True:
    choice = input('Calculating:\nChoose the Satellite:  1) NOAA 11  2) Molniya\n')
    choice = int(choice)
    if choice == 1:
        M0          = 192.28166              #[deg]
        bomega0     = 29.31059               #[deg]
        somega0     = 167.74754              #[deg]
        epsilon     = 0.00119958 
        i           = 98.974460              #[deg]
        T           = 102.0764*60.           #[sec]
        a           = 7229.606               #[km] 
        break
    elif choice == 2:
        M0          = 40.                    #[deg]
        bomega0     = 0.                     #[deg]
        somega0     = 270.                   #[deg]
        epsilon     = 0.737 
        i           = 63.4                   #[deg]
        T           = 717.7*60.              #[sec]
        a           = 26553.                 #[km] 
        break
    else:
        print("Invalid value! Please do it again.")

origin = datetime(1990, 1,  1, 0,  0, 0, 0) 
star_t = datetime(1990, 3, 22, 2, 58, 0, 0)
end_t  = datetime(1990, 3, 22, 8,  4, 0, 0)
ree         = 6378.140               #[km]
j2          = 0.00198263
total_time  = 3*T
nT          = 10000
dT          = total_time/nT
deg2rad     = np.pi/180.
omega_GW0   = 100.38641              #[deg]; 1990/01/01 00:00 UTC
ro_rate     = 7.29*10**-5/deg2rad    #[deg/sec]
lats        = np.zeros(nT)
longs       = np.zeros(nT)
omega_GW    = np.zeros(nT)
x3          = np.zeros(nT)
y3          = np.zeros(nT)
z3          = np.zeros(nT)
omegas      = np.zeros(nT)

#---
n          = 2*np.pi/T  #[radians/sec]  
n_bar      = n*(1+3/2*j2*(ree/a)**2*(1-epsilon**2)**(-3/2)*(1-3/2*(np.sin(i*deg2rad))**2)) #[radians/sec] 
dM_dt      = n_bar/deg2rad #[deg/sec] 
dbomega_dt = (-n_bar*(3/2*j2*(ree/a)**2*(1-epsilon**2)**(-2)*np.cos(i*deg2rad)))/deg2rad #[deg/sec] 
dsomega_dt = (n_bar*(3/2*j2*(ree/a)**2*(1-epsilon**2)**(-2)*(2-5/2*(np.sin(i*deg2rad))**2)))/deg2rad #[deg/sec] 
 
#--------------
# initialize
#--------------
time = ((star_t-origin).days)*86400+((star_t-origin).seconds)
omega_GW[0] = (omega_GW0 + ro_rate*time)%360.  # omega(greenwich) at 1990/03/22 02:58 UTC
if omega_GW[0]> 180.:
    omega_GW[0] -= 360.
elif omega_GW[0] < -180:
    omega_GW[0] += 360.
    
M        = M0
bomega   = bomega0
somega   = somega0


iT = 0
while iT < nT:
    time = time + dT
    omega_GW[iT] = (omega_GW0 + ro_rate*time)%360.
    if omega_GW[iT] > 180.:
        omega_GW[iT] -= 360.
    elif omega_GW[iT] < -180:
        omega_GW[iT] += 360.
    
    M      += dM_dt*dT
    bomega += dbomega_dt*dT
    somega +=  dsomega_dt*dT
    
    # eccentric anomaly, e[deg]
    if choice == 1:
        e = M + ((2*epsilon-1/4*epsilon**3)*np.sin(M*deg2rad) + 5/4*epsilon**2*np.sin(2*M*deg2rad) + 13/12*epsilon**3*np.sin(3*M*deg2rad))/deg2rad
    elif choice == 2:
        e = 0. # [rad]
        e_increment = 2*np.pi/100. # [rad]
        while True:
            diff = np.radians(M) - (e-epsilon*np.sin(e))
            if diff > 1e-10:
                e = e + e_increment # [rad]
            else:
                e = e/deg2rad #[rad -> deg]
                #print("iteration done!")
                break
    else:
        print("something wrong!")
    # true anomaly, theta[deg]
    theta = 2*(np.arctan(((1+epsilon)/(1-epsilon))**(1/2)*np.tan(e*deg2rad/2)))/deg2rad
    # radius[km]
    r =  a*(1-epsilon**2)/(1+epsilon*np.cos(theta*deg2rad)) 
    # initial position in the plane of its orbit 
    x0 = r*np.cos(theta*deg2rad)
    y0 = r*np.sin(theta*deg2rad)
    z0 = 0.
    # rotate to the exact position
    # 1st rotation: small omega
    x1 = x0*np.cos(somega*deg2rad) - y0*np.sin(somega*deg2rad)
    y1 = x0*np.sin(somega*deg2rad) + y0*np.cos(somega*deg2rad)
    z1 = z0
    # 2nd rotation: inclination
    x2 = x1
    y2 = y1*np.cos(i*deg2rad) - z1*np.sin(i*deg2rad)
    z2 = y1*np.sin(i*deg2rad) + z1*np.cos(i*deg2rad)
    # 3rd rotation: big omega
    x3[iT] = x2*np.cos(bomega*deg2rad) - y2*np.sin(bomega*deg2rad)
    y3[iT] = x2*np.sin(bomega*deg2rad) + y2*np.cos(bomega*deg2rad)
    z3[iT] = z2
    
    # convert to spherical coordinate
    rs     = (x3[iT]**2+y3[iT]**2+z3[iT]**2)**(1/2)
    deltas = np.arcsin(z3[iT]/rs)/deg2rad
    omegas[iT] = np.arctan2(y3[iT],x3[iT])/deg2rad 
    
    lats[iT]   = deltas
    longs[iT]  = omegas[iT] - omega_GW[iT]
    if longs[iT] > 180.:
        longs[iT] -= 360.
    elif longs[iT] < -180:
        longs[iT] += 360.

    iT += 1

# save the data
if choice == 1:
    np.savetxt('lats.txt', lats)
    np.savetxt('longs.txt', longs)
    np.savetxt('x3.txt', x3[0:3333])
    np.savetxt('y3.txt', y3[0:3333])
    np.savetxt('z3.txt', z3[0:3333])
    print("NOAA 11 data saved!")
else:
    np.savetxt('lats_2.txt', lats)
    np.savetxt('longs_2.txt', longs)
    np.savetxt('x3_2.txt', x3[0:3333])
    np.savetxt('y3_2.txt', y3[0:3333])
    np.savetxt('z3_2.txt', z3[0:3333])
    print("Molniya data saved!")
