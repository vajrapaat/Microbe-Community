import numpy as np
import math

bounds = 500
Dn = 250
bc = 10
Dc = 1000
an = 1
k = 0.333
gridspace = 301

points = int(2*bounds)
coord = np.linspace(-bounds, bounds, gridspace)
dx = coord[1]-coord[0]
bins = int(len(coord) - 1.0)

def nutrientfield():
    return np.ones((bins,bins))

def chemmfield():
    return np.zeros((bins,bins))

def density(cells,cell_area,specfilt = None):
    xs = []
    ys = []
    areas = []
    for c in cells:
        if specfilt is not None and c[5] != specfilt:
            continue
        xs.append(c[0])
        ys.append(c[1])
        areas.append(cell_area(c))
    if not xs:
        return np.zeros((bins, bins))
    hist, _, _ = np.histogram2d(xs,ys,bins=[coord, coord], weights = areas)
    gridcellar = dx**2 
    return hist/gridcellar

def laplacian(field):
    lap = np.zeros_like(field)
    lap[1:-1,1:-1] = (field[2:,1:-1] + field[:-2,1:-1] + field[1:-1,2:] + field[1:-1,:-2] - 4*field[1:-1,1:-1])/dx**2
    return lap

def substeps(D, dt):
    lim = dx**2/(4*D)
    if dt <= lim:
        return 1
    return int(math.ceil(dt/lim))
    
def nutrientstep(n_field, rho, dt):
    substep = substeps(Dn, dt)
    sub_dt = dt/substep
    for _ in range(substep):
        lap = laplacian(n_field)
        cons = an*rho*n_field/(k+n_field)
        n_field = n_field + (Dn*lap - cons)*sub_dt
        n_field[0,:] = 1.0
        n_field[-1,:] = 1.0
        n_field[:,0] = 1.0
        n_field[:,-1] = 1.0
        n_field = np.clip(n_field, 0.0, None) 
    return n_field

def chemmstep(c_field, rhoc, dt):
    substepc = substeps(Dc, dt)
    sub_dt = dt/substepc
    for _ in range(substepc):
        lap = laplacian(c_field)
        c_field = c_field + (Dc*lap - bc*c_field + rhoc)*sub_dt
        c_field[0,:] = 0.0
        c_field[-1,:] = 0.0
        c_field[:,0] = 0.0
        c_field[:,-1] = 0.0
        c_field = np.clip(c_field, 0.0, None)
    return c_field

def lookup(field, x, y):
    ix = int((x+bounds)/dx)
    iy = int((y+bounds)/dx)
    ix = min(max(ix, 0), bins-1)
    iy = min(max(iy, 0), bins-1)
    return field[ix,iy]
    
    
