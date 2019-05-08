import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from read_data import read

def get_end_site(x,y,index):
    if index==0: # clockwise up
        return (x,y+1)
    if index==1:
        return (x+1,y)
    if index==2:
        return (x,y-1)
    if index==3:
        return (x-1,y)

def plot(bonds, f = None, a = None, color=None):
    if color==None:
        colors=[cm.viridis(i) for i in np.linspace(0,1,len(bonds))]
    else:
        colors=[color for i in range(len(bonds))]
    for burst,color in zip(bonds,colors):
        end_points = [get_end_site(*i) for i in burst]
        start_points = [c[0:2] for c in burst]
        if f==None and a == None:
            f,a = plt.subplots(1,1, figsize = [10,10], dpi=100)
        for s,e in zip(start_points,end_points):
            dx = [s[0],e[0]]
            dy = [s[1],e[1]]
            a.plot(dx,dy,color=color)
        try:
            a.plot(*burst[0][0:2],'o', color='k', fillstyle='none')
        except IndexError:
            pass
        a.axis('equal')
    return f,a

def bigplot(bonds):
    sites = [(site[0],site[1]) for burst in bonds for site in burst] #flatten
    maxx = max(sites, key = lambda s: s[0])[0]
    maxy = max(sites, key = lambda s: s[1])[1]
    minx = min(sites, key = lambda s: s[0])[0]
    miny = min(sites, key = lambda s: s[1])[1]
    dx = maxx-minx+2
    dy = maxy-miny+2

    M=np.zeros((dx,dy)) #center is at (0,0) and positions are allowed negative
    
    for (x,y) in sites:
        M[(x-minx)][(y-miny)] = 1
    
    M=1-M
    
    dpi = 1000.0
    xpixels, ypixels = maxx, maxy
    
    fig = plt.figure(figsize=(ypixels/dpi, xpixels/dpi), dpi=dpi)
    fig.figimage(M, origin='lower', cmap='gray', resize=True)
    return fig

ns, bonds = read('../small')
f,a = plot(bonds)
f.savefig('small.png')
