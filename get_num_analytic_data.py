#these functions are for the getting the data required for comparing the analytical
#and numerical results
from matplotlib import pyplot as plt
import numpy as np
from scipy.special import erf
import os


def get_phi_num_slab(phi_n_num,nx,ny):
    #we go for the middle of the slab
    mid_y = ny//2
    concent_data = np.zeros(nx-1)
    #we get the whole x coordinate w/ midplane
    concent_data = phi_n_num[mid_y,1:-1]
    return concent_data

def solve_inf_slab_analytical(t_current,dx,nx,r,x,C0,D):
    concent_data = np.zeros(nx-1)
    x_mid = dx*nx//2

    x0 = np.asanyarray((x-x_mid),dtype=float)
    if t_current <= 0.0:
        concent_data = C0 * (x0).astype(float)
        return concent_data
    denom = 2.0 * np.sqrt(D * t_current)
    concent_data = 2+(0.5*C0 * (erf((x0 + r)/denom) - erf((x0 - r)/denom)))
    return concent_data

def comparision(x,concent_data_num,concent_data_an,i):
    save_location = "C:\\Users\\mehme\\OneDrive\\Desktop\\projects\\2d\\adv_pth\\validation_cases"
    if not save_location:
        os.mkdir("C:\\Users\\mehme\\OneDrive\\Desktop\\projects\\2d\\adv_pth\\validation_cases")
    fig,(ax1,ax2) = plt.subplots(2,1)
    ax1.plot(x,concent_data_num)
    ax2.plot(x,concent_data_an)
    fig.align_ylabels()
    fig.savefig(save_location+"\\vtk_output\\"+f"data_{i}"+".png",format="png")
    