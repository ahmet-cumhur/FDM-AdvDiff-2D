from numba import njit, prange

@njit(parallel=True)
def get_phi(ny,nx,u,v,phi_n,phi_s,dt,dx,dy,nu):
        for j in prange(1,ny+1):
            for i in range(1,nx+1):
                adv_x = (u/dx)*(phi_n[j,i]-phi_n[j,i-1])
                adv_y = (v/dy)*(phi_n[j,i]-phi_n[j-1,i])
                diff_x = (nu/(dx*dx))*(phi_n[j,i-1]-2*phi_n[j,i]+phi_n[j,i+1])
                diff_y = (nu/(dy*dy))*(phi_n[j-1,i]-2*phi_n[j,i]+phi_n[j+1,i])
                phi_s[j,i] = phi_n[j,i]-dt*(adv_x+adv_y)+dt*(diff_x+diff_y)
                #dont forget to call apply_bc_s function after
        return phi_s
# needs a recheck!
def get_phi_vectorized(u,v,phi_n,phi_s,dt,dx,dy,nu,ny,nx,phi_inlet):
        adv_x_vector = (u/dx)*(phi_n[:,:]-phi_n[:,:-1])
        adv_y_vector = (v/dy)*(phi_n[:,:]-phi_n[:-1,:])
        diff_x_vector = (nu/(dx*dx))*(phi_n[:,:-1]-2*phi_n[:,:]+phi_n[:,:+1])
        diff_y_vector = (nu/(dy*dy))*(phi_n[:-1,:]-2*phi_n[:,:]+phi_n[:+1,:])
        phi_s[:,:] = phi_n[:,:]-dt*(adv_x_vector+adv_y_vector)+dt*(diff_x_vector+diff_y_vector)
        return phi_s

@njit(parallel=True)
def get_phi_central(ny,nx,u,v,phi_n,phi_s,dt,dx,dy,nu):
       for j in prange(1,ny+1):
            for i in range(1,nx+1):
                adv_x = (u/(2*dx))*(phi_n[j,i+1]-phi_n[j,i-1])
                adv_y = (v/(2*dy))*(phi_n[j+1,i]-phi_n[j-1,i])
                diff_x = (nu/(dx*dx))*(phi_n[j,i-1]-2*phi_n[j,i]+phi_n[j,i+1])
                diff_y = (nu/(dy*dy))*(phi_n[j-1,i]-2*phi_n[j,i]+phi_n[j+1,i])
                phi_s[j,i] = phi_n[j,i]-dt*(adv_x+adv_y)+dt*(diff_x+diff_y)
        
       return phi_s

def apply_bc_n(phi_n,phi_inlet,nx,ny):
        #add the bc
        #southern bc neuman
        phi_n[0, :] = phi_n[1, :]
        #northern bc neuman
        phi_n[ny+1, :] = phi_n[ny, :]   
        #eastern or outlet bc neuman 
        phi_n[:, nx+1] = phi_n[:, nx]
        #western bc neuman
        phi_n[:, 0] = 2.0*phi_inlet - phi_n[:, 1]
        
        return phi_n
def apply_bc_s(phi_s,phi_inlet,nx,ny):
        #add the bc
        #southern bc neuman
        phi_s[0, :] = phi_s[1, :]
        #northern bc neuman
        phi_s[ny+1, :] = phi_s[ny, :]   
        #eastern or outlet bc neuman 
        phi_s[:, nx+1] = phi_s[:, nx]
        #western bc neuman
        phi_s[:, 0] = 2.0*phi_inlet - phi_s[:, 1]
        
        return phi_s