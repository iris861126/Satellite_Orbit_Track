#=====================================
# plotting code I: the groundtrack
# load data: lats, longs
#=====================================
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

while True:
    choice = input('Plot the groundtrack:\nChoose the Satellite:  1) NOAA 11  2) Molniya\n')
    choice = int(choice)
    if choice == 1:
        break
    elif choice == 2:
        break
    else:
        print("Invalid value! Please do it again.")
        
# load data
if choice == 1:
    lats  = np.genfromtxt('lats.txt')
    longs = np.genfromtxt('longs.txt')
    print("NOAA 11 data loaded!")
else:
    lats  = np.genfromtxt('lats_2.txt')
    longs = np.genfromtxt('longs_2.txt')
    print("Molniya loaded!")
nT = np.size(lats)
#--------------
# draw the map
#--------------
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines(resolution='50m',color='gray')
ax.add_feature(cfeature.LAKES, edgecolor='gray', color='whitesmoke')
ax.add_feature(cfeature.OCEAN,color='whitesmoke')
ax.gridlines(linestyle='--')


# axes
plt.xlim(-180,180)
plt.ylim(-90,90)
tick_proj = ccrs.PlateCarree()
ax.set_xticks(np.arange(-180, 180+60, 60), crs=tick_proj)
ax.set_xticks(np.arange(-180, 180+30, 30), minor=True, crs=tick_proj)
ax.set_yticks(np.arange(-90, 90+30, 30), crs=tick_proj)
ax.set_yticks(np.arange(-90, 90+15, 15), minor=True, crs=tick_proj)
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())

#----------------
# draw the track
#----------------

if choice == 1:
    pts = np.zeros(nT)
    count = 0
    for k in range(1,nT):
         if abs(longs[k]-longs[k-1]) >= 180.:    
            pts[count] = k
            count += 1
    pts = pts[pts >0]
    ax.plot(longs[0:int(pts[0])],lats[0:int(pts[0])],color="black")
    ax.plot(longs[int(pts[count-1]):nT-1],lats[int(pts[count-1]):nT-1],color="black")
    for k in range(count-1):
        ax.plot(longs[int(pts[k]):int(pts[k+1])],lats[int(pts[k]):int(pts[k+1])],color="black")
    
    # plot the start pt, mid point & end pt
    plt.annotate('Start',fontsize=15,xy=(longs[0],lats[0]),xytext=(longs[0]+7,lats[0]+2),arrowprops={'headwidth':10,'facecolor':'salmon'})
    plt.annotate('End',fontsize=15,xy=(longs[nT-1],lats[nT-1]),xytext=(longs[nT-1]+7,lats[nT-1]+2),arrowprops={'headwidth':10,'facecolor':'steelblue'})
    plt.annotate('Mid',fontsize=15,xy=(longs[int(nT/2)],lats[int(nT/2)]),xytext=(longs[int(nT/2)]+7,lats[int(nT/2)]+2),arrowprops={'headwidth':10,'facecolor':'olive'})
    plt.scatter(longs[0],lats[0],color='salmon',s=50)
    plt.scatter(longs[nT-1],lats[nT-1],color='steelblue',s=50)
    plt.scatter(longs[int(nT/2)],lats[int(nT/2)],color='olive',s=50)
    # the title
    plt.title("NOAA II Groundtrack",fontsize=25)
    plt.title("1990.03.22",fontsize=10,loc='left')
    plt.title("Start time: 0258UTC   End time: 0804UTC",fontsize=10,loc='right')
else:
    plt.scatter(longs,lats,s=1.2,color= 'black')
    plt.arrow(longs[50],lats[50],longs[125]-longs[50],lats[125]-lats[50],width=2,overhang=0.5,color='olive')
    
    # plot the start pt & end pt
    plt.annotate('Start',fontsize=15,xy=(longs[0],lats[0]),xytext=(longs[0]+7,lats[0]+2),arrowprops={'headwidth':10,'facecolor':'salmon'})
    plt.annotate('End',fontsize=15,xy=(longs[nT-1],lats[nT-1]),xytext=(longs[nT-1]+7,lats[nT-1]+2),arrowprops={'headwidth':10,'facecolor':'steelblue'})
    plt.scatter(longs[0],lats[0],color='salmon',s=50)
    plt.scatter(longs[nT-1],lats[nT-1],color='steelblue',s=50)
    # the title
    plt.title("Molniya Groundtrack",fontsize=25)
    plt.title("Total Time: 3 periods",fontsize=10,loc='right')

plt.show()

