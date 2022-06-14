#=====================================
# plotting code II: the orbit
# load data: lats, longs
#=====================================
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.axes3d import Axes3D
import cartopy.crs as ccrs


fig = plt.figure()
ax = Axes3D(fig)
# draw earth 
img = plt.imread("earthmap5.png")

u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j] 
xe = 6360*np.cos(u)*np.sin(v) 
ye = 6360*np.sin(u)*np.sin(v) 
ze = 6360*np.cos(v) 
ax.plot_wireframe(xe, ye, ze, color="gray", linewidths=0.8, zorder=1) 
ax.scatter(0, 0, 0, c ='black', s=30, zorder=1)

while True:
    choice = input('Plot the 3d-orbit:\nChoose the Satellite:  1) NOAA 11  2) Molniya\n')
    choice = int(choice)
    if choice == 1:
        break
    elif choice == 2:
        break
    else:
        print("Invalid value! Please do it again.")
        
# load data
if choice == 1:
        x = np.genfromtxt('x3.txt')
        y = np.genfromtxt('y3.txt')
        z = np.genfromtxt('z3.txt')
        print("NOAA 11 data loaded!")
else:
        x = np.genfromtxt('x3_2.txt')
        y = np.genfromtxt('y3_2.txt')
        z = np.genfromtxt('z3_2.txt')
        print("Molniya data loaded!")
 
# draw the orbit
ccc = np.arange(3333)
ax.scatter(x, y, z, c = ccc,cmap='gist_rainbow', s=5, zorder=5)  # start from red
ax.scatter(0,0,6357, c ='black', s=30, zorder=1)
ax.text(0,0,6357, "North Pole", size=20, zorder=9, color='black',fontsize=20)

# plot the start pt, mid point & end pt
ax.scatter(x[0],y[0],z[0],zorder=7,color='salmon',s=80)
ax.text(x[0]+11,y[0]+11,z[0]+11, "Start", size=20, zorder=9, color='salmon',fontsize=20)
ax.scatter(x[1666],y[1666],z[1666],zorder=6, color='olive', s=80)
ax.text(x[1666]+11,y[1666]+11,z[1666]+11, "Mid", size=20, zorder=10, color='olive',fontsize=20)
ax.scatter(x[3332],y[3332],z[3332],zorder=8,color='steelblue',s=80)
ax.text(x[3332]+11,y[3332]-1100,z[3332]-1100, "End", size=20, zorder=13, color='steelblue',fontsize=20)

ax.set_xlabel('x [km]')
ax.set_ylabel('y [km]')
ax.set_zlabel('z [km]')
plt.show()
