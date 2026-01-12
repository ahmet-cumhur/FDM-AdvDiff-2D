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
    r = r*dx
    concent_data = np.zeros(nx-1)
    x_in = x[1:-1]
    x_center = 0.5*(x_in[0]+x_in[-1])

    x0 = x - x_center
    if t_current <= 0.0:
        concent_data = C0 * (np.abs(x0)<=r).astype(float)
        return concent_data
    denom = 2.0 * np.sqrt(D * t_current)
    concent_data = (0.5*C0 * (erf((x0 + r)/denom) - erf((x0 - r)/denom)))
    return concent_data
#we need to store l2 error for each time than plot it 
def get_L2_error(concent_data_num,concent_data_an):
    error = concent_data_an-concent_data_num
    l2 = np.sqrt(np.mean(error**2))
    l2_rel = l2/(np.sqrt(np.mean(concent_data_an**2)) + 1e-14)
    linf = np.max(np.abs(error))
    return l2,l2_rel,linf

def comparision_concentration(x,concent_data_num,concent_data_an,i,C):
    save_location = "C:\\Users\\mehme\\OneDrive\\Desktop\\projects\\2d\\adv_pth\\validation_cases"
    if not save_location:
        os.mkdir("C:\\Users\\mehme\\OneDrive\\Desktop\\projects\\2d\\adv_pth\\validation_cases")
    
    fig,(ax1,ax2,ax3) = plt.subplots(3,1,figsize=(12,10),sharex=True)

    ax1.plot(x,concent_data_num,"b")
    ax1.set_title("Numerical Solution")
    ax1.set_xlabel("Domain length[m]")
    ax1.set_ylabel("Concentration[-]")
    ax1.set_ylim(-0.05,C*1.05)

    ax2.plot(x,concent_data_an,"r")
    ax2.set_title("Analytical Solution")
    ax2.set_xlabel("Domain length[m]")
    ax2.set_ylabel("Concentration[-]")
    ax2.set_ylim(-0.05,C*1.05)

    ax3.plot(x,concent_data_num,"b",label="Numerical Solution")
    ax3.plot(x,concent_data_an,"r",label="Analytical Solution")
    ax3.set_xlabel("Domain length[m]")
    ax3.set_ylabel("Concentration[-]")
    ax3.set_title("Numerical and Analytical Solution")
    ax3.set_ylim(-0.05,C*1.05)
    
    fig.legend()
    fig.tight_layout()
    fig.savefig(save_location+"\\vtk_output\\"+f"val_{i}"+".png",format="png")
    plt.close(fig)

def err_vis(t_current,l2,l2_rel,linf):
    save_location = "C:\\Users\\mehme\\OneDrive\\Desktop\\projects\\2d\\adv_pth\\validation_cases"
    fig,ax = plt.subplots()

    ax.plot(t_current,l2,"b",label="L2 Error")
    ax.plot(t_current,l2_rel,"r",label="L2 Relative Error")
    ax.plot(t_current,linf,"g",label="Linf Error")

    ax.set_xlabel("Time[s]")
    ax.set_ylabel("Error[%]")
    ax.set_title("1D Slab Diffusion error")

    fig.legend()
    fig.tight_layout()
    fig.savefig(save_location+"\\vtk_output\\"+"err_comp"+".png",format="png")
    plt.close(fig)


    