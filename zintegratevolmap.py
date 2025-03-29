from gridData import Grid
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d

g = Grid("ca1_6ca.dx") #dx file name

mz= int(g.grid.shape[2]) - 1 #number of z slides

nz = 0 #integration of slides in Z
iz = g.grid[:,:,nz]
nz = 1
while nz<mz: 
    iz = iz + g.grid[:,:,nz]
    nz = nz + 1

f = np.array(iz)
cx= g.origin[0]
cy= g.origin[1]

for index, z in np.ndenumerate(f):  
    x = float(index[0] + cx) 
    y = float (index[1]+ cy)
    print(x, y, z) 






      
    
        

