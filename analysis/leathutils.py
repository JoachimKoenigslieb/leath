import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.cm as cm

def get_end_site(x,y,index):
    if index==0: # clockwise up
        return (x,y+1)
    if index==1:
        return (x+1,y)
    if index==2:
        return (x,y-1)
    if index==3:
        return (x-1,y)

def bigplot(bonds):
    sites = [(site[0],site[1]) for burst in bonds for site in burst] #flatten
    maxx = abs(max(sites, key = lambda s: abs(s[0]))[0])
    maxy = abs(max(sites, key = lambda s: abs(s[1]))[1])
    M=np.zeros((2*maxx+20,2*maxy+20)) #center is at (0,0) and positions are allowed negative
    for (x,y) in sites:
        M[x+maxx][y+maxy]=1
    f,a = plt.subplots(1)
    a.imshow(M.transpose(),cmap='gray', origin='lower')
    return f,a

def plot3d(bonds):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for burst in bonds:
        x,y,z = zip(*[(i[0],i[1],i[2]) for i in burst])
        ax.scatter(x,y,z, s=1)
    return fig,ax

def plot(bonds, f = None, a = None):
    colors=[cm.viridis(i) for i in np.linspace(0,1,len(bonds))]
    for burst,color in zip(bonds,colors):
        end_points = [get_end_site(*i) for i in burst]
        start_points = [c[0:2] for c in burst]
        if f==None and a == None:
            f,a = plt.subplots(1,1)
        for s,e in zip(start_points,end_points):
            dx = [s[0],e[0]]
            dy = [s[1],e[1]]
            a.plot(dx,dy,color=color)
        try:
            a.plot(*burst[0][0:2],'o',color=color)
        except IndexError:
            pass
        a.axis('equal')

