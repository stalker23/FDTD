import numpy as np
from matplotlib import pyplot as plt 
import matplotlib.cm as cm 
import os
import timeit

startTime = timeit.default_timer()
###########################################################
###################### Problem Setup ######################
###########################################################
# Define Grid
# Define Domain
Nx = 15
Ny = 5
Lx = 4.0 # m
Ly = 8.0 # m
dx = Lx/(Nx-1) # m
dy = Ly/(Ny-1) # m
R_circle = np.array([2,4]) # m
radius_circle = 1 # m
dt = 0.01 # sec
numSteps = 100 # number of time steps
printSteps = 1;
figureDirectory = os.path.realpath('./')+'\\Figures\\step'

# Define Parameter Relations
c0 = 299792458 # m/s
# Permeability Water
mu_xx_H2O = 1;
mu_yy_H2O = 1;
mu_zz_H2O = 1;
# Permeability Gold
mu_xx_Gold = 1;
mu_yy_Gold = 1;
mu_zz_Gold = 1;
# Permittivity Water
epsilon_xx_H2O = 1
epsilon_yy_H2O = 1
epsilon_zz_H2O = 1
# Permittivity Gold
epsilon_xx_Gold = 1+1j
epsilon_yy_Gold = 1+1j
epsilon_zz_Gold = 1+1j

# Initialize Arrays
x = np.linspace(0,Lx,Nx+2)
y = np.linspace(0,Ly,Ny+2)
X, Y = np.meshgrid(x, y)
Hx_old = np.zeros((Ny+2,Nx+2))
Hy_old = np.zeros((Ny+2,Nx+2))
Hz_old = np.zeros((Ny+2,Nx+2))
Ex_old = np.zeros((Ny+2,Nx+2))
Ey_old = np.zeros((Ny+2,Nx+2))
Ez_old = np.zeros((Ny+2,Nx+2))
Hx = np.zeros((Ny+2,Nx+2))
Hy = np.zeros((Ny+2,Nx+2))
Hz = np.zeros((Ny+2,Nx+2))
Ex = np.zeros((Ny+2,Nx+2))
Ey = np.zeros((Ny+2,Nx+2))
Ez = np.zeros((Ny+2,Nx+2))
m_Hxx = np.zeros((Ny+2,Nx+2))
m_Hyy = np.zeros((Ny+2,Nx+2))
m_Hzz = np.zeros((Ny+2,Nx+2))
m_Exx = np.zeros((Ny+2,Nx+2))
m_Eyy = np.zeros((Ny+2,Nx+2))
m_Ezz = np.zeros((Ny+2,Nx+2))

# ICs and Define Domain Variables
for j in range(1,Nx-1):
    for i in range(1,Ny-1):
        # Calculate if point is in or on the circle
        R_point = np.array([X[i][j],Y[i][j]]) # Outside Circle
        if np.linalg.norm(np.subtract(R_point,R_circle)) > radius_circle:
            # ICs 
            # Everythings allready ICed to zero
            
            # Assign Parameters to points
            m_Hxx = -c0*dt/mu_xx_H2O
            m_Hyy = -c0*dt/mu_yy_H2O
            m_Hzz = -c0*dt/mu_zz_H2O
            m_Exx = c0*dt/epsilon_xx_H2O
            m_Eyy = c0*dt/epsilon_yy_H2O
            m_Ezz = c0*dt/epsilon_zz_H2O
            
        else: # Inside Circle
            # ICs 
            # Everythings allready ICed to zero
            
            # Assign Parameters to points
            m_Hxx = -c0*dt/mu_xx_Gold
            m_Hyy = -c0*dt/mu_yy_Gold
            m_Hzz = -c0*dt/mu_zz_Gold
            m_Exx = c0*dt/epsilon_xx_Gold
            m_Eyy = c0*dt/epsilon_yy_Gold
            m_Ezz = c0*dt/epsilon_zz_Gold
            
# Apply BCs
Hic = 5;
Eic = 5;
for i in range(1, Ny-1):
    j = 1
    Hx[i][j] = Hic
    Hy[i][j] = Hic
    Hz[i][j] = Hic
    Ex[i][j] = Eic
    Ey[i][j] = Eic
    Ez[i][j] = Eic
 
###########################################################           
###################### Iterate the Solution ###############
###########################################################
for k in range(0,numSteps):
    print(k)
    # Update Time Varying BCs
    
    # Apply Periodic BC on bottom edge 
    for j in range(1, Nx-1):
        i = 1
        Hx[i-1][j] = Hx[i][j]
        Hz[i-1][j] = Hz[i][j]
    # Apply Periodic BC on top edge 
    for j in range(1, Nx-1):
        i = Ny
        Ex[i+1][j] = Ex[i][j]
        Ez[i+1][j] = Ez[i][j]
    # Apply Periodic BC on left edge 
    for i in range(1, Ny-1):
        j = 1
        Hy[i][j-1] = Hy[i][j]
        Hz[i][j-1] = Hz[i][j]
     # Apply Periodic BC on right edge 
    for i in range(1, Ny-1):
        j = Nx
        Ey[i][j+1] = Ey[i][j]
        Ez[i][j+1] = Ez[i][j]

    for j in range(1,Nx-1):
        for i in range(1,Ny-1):
            # Calculate differences
            Ezy = (Ez[i+1][j]-Ez[i][j])/dy
            Ezx = (Ez[i][j+1]-Ez[i][j])/dx
            Eyx = (Ey[i][j+1]-Ey[i][j])/dx
            Exy = (Ex[i+1][j]-Ex[i][j])/dy
            Hzy = (Hz[i][j]-Hz[i-1][j])/dy
            Hzx = (Hz[i][j]-Hz[i][j-1])/dx
            Hyx = (Hy[i][j]-Hy[i][j-1])/dx
            Hxy = (Hx[i][j]-Hx[i-1][j])/dy                  
            # Update H from E
            Hx[i][j] = m_Hxx*(Ezy) + Hx_old[i][j]
            Hy[i][j] = m_Hyy*(-Ezx) + Hy_old[i][j]
            Hz[i][j] = m_Hzz*(Eyx-Exy) + Hz_old[i][j]
            # Update E from H
            Ex[i][j] = m_Exx*(Hzy) + Ex_old[i][j]
            Ey[i][j] = m_Eyy*(-Hzx) + Ey_old[i][j]
            Ez[i][j] = m_Ezz*(Hyx-Hxy) + Ez_old[i][j]
            
            # Print figure
            if (k+printSteps)%printSteps == 0:
                plt.pcolormesh(X, Y, Ex, cmap = cm.jet) 
                plt.savefig(figureDirectory + str(k) + '.png')
                
stopTime = timeit.default_timer()
runTime = stopTime-startTime
            
                    




