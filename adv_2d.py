import numpy as np
import os

class Adv_Diff_Solver:
    def __init__(self,nx,ny,l,h,nu,u,v,phi_inlet):
        self.nx = nx
        self.ny = ny
        self.l = l
        self.h = h
        self.nu = nu
        self.u = u
        self.v = v
        self.phi_inlet = phi_inlet
        self.phi_n = np.zeros((ny+2,nx+2))
        self.phi_s = np.zeros((ny+2,nx+2))
        self.dx = l/nx
        self.dy = h/ny
        self.x = np.linspace(0,l,nx)
        self.y = np.linspace(0,h,ny)

    def time_roll(self,phi_n,phi_s):
        phi_n[:, :] = phi_s[:, :]
        return phi_n

    def apply_ic(self,phi_n,phi_s,nx,ny,r,validation_case,concentration):
        concentration=float(concentration)
        mid_x = nx//2
        mid_y = ny//2
        if validation_case == "source":
            phi_s[mid_y-r:mid_y,mid_x-r:mid_x+r] = concentration
            phi_n[mid_y-r:mid_y,mid_x-r:mid_x+r] = concentration
        elif validation_case == "slab":
            # +1 here comes for apply the ic in the middle otherwise
            #we see a 0.5 dx shift in the ic
            phi_s[:,mid_x-r+1:mid_x+r+1] = concentration
            phi_n[:,mid_x-r+1:mid_x+r+1] = concentration
        elif validation_case == "empty":
            phi_s[mid_y-r:mid_y+r,mid_x-r:mid_x+r] = concentration
            phi_n[mid_y-r:mid_y+r,mid_x-r:mid_x+r] = concentration
        return phi_s,phi_n
    
    def calculate_time_step(self,dx,dy,u,v,nu):
        eps = 1e-10
        coef_four = 0.2
        coef_cfl = 0.5
        dt_cfl_x = coef_cfl*dx/np.abs(u)
        if dt_cfl_x<eps:
            dt_cfl_x= max(1e10,dt_cfl_x)
        dt_cfl_y = coef_cfl*dy/np.abs(v)
        if dt_cfl_y<eps:
            dt_cfl_y= max(1e10,dt_cfl_y)

        dt_four_x = coef_four*(dx**2)/nu
        dt_four_y = coef_four*(dy**2)/nu
        dt = max(eps,min(dt_cfl_x,dt_cfl_y,dt_four_x,dt_four_y))
        return dt
    
    def calculate_time_step_diffusion(self,dx,dy,nu):
        eps = 1e-10
        coef_four = 0.05
        dt_four_x = coef_four*(dx**2)/nu
        dt_four_y = coef_four*(dy**2)/nu
        dt = max(eps,min(dt_four_x,dt_four_y))
        print("time step size x: ", dt_four_x,"time step size y: ",dt_four_y)
        return dt
    
    def calculate_time_step_2(self,dx,dy,u,v,nu):
        eps = 1e-10
        coef_four = 0.2
        coef_cfl = 0.5
        
        if np.abs(u)<eps:
            dt_cfl_x= 1e10
        else:
            dt_cfl_x = coef_cfl*dx/np.abs(u)
        if np.abs(v)<eps:
            dt_cfl_y= 1e10
        else:
            dt_cfl_y = coef_cfl*dy/np.abs(v)

        dt_four_x = coef_four*(dx**2)/nu
        dt_four_y = coef_four*(dy**2)/nu
        dt = max(eps,min(dt_cfl_x,dt_cfl_y,dt_four_x,dt_four_y))
        return dt
    
    def time_counter(self,t0,t1):
        return t1-t0
        
    def vtk_converter(self,nx,ny,phi_n,dx,dy,t_current,i):
        
        x0,y0 = 0.0,0.0
        npts = nx*ny
        work_dir = os.getcwd()
        if not work_dir:
            os.makedirs(work_dir+"\\vtk_output")
        filename = work_dir+f"\\vtk_output\\data_{i:07d}.vtk"
        temp_filename = work_dir+f"\\vtk_output\\data_{i:07d}.tmp"
        x0c = x0+0.5*dx; y0c = y0+0.5*dy
        
        # first we need to change the nx,ny to ny,nx and than,
        # we dont need the whole number so we can save the data as float32 which takes less space
        # than we choose how we compress the file like "C" or "F".
        phi_inner=phi_n[1:ny+1,1:nx+1]
        #phi_n has ny+2,nx+2 size because of it ghost cell application therefore,
        #we need to only print the inner cells
        phi_flat = phi_inner.ravel(order="C")
        
        f = open(temp_filename,"w")
        # Header
        f.write("# vtk DataFile Version 3.0\n")
        f.write("cfd-python 2D fields\n")
        f.write("ASCII\n")
        f.write("DATASET STRUCTURED_POINTS\n")
        f.write(f"DIMENSIONS {nx} {ny} 1\n")
        f.write(f"ORIGIN {x0c:.8g} {y0c:.8g} 0\n")
        f.write(f"SPACING {dx:.8g} {dy:.8g} 1\n")
        
        # ---- DATASET-LEVEL FIELD (time) ----
        f.write("FIELD FieldData 1\n")
        f.write("TIME 1 1 float\n")
        f.write(f"{t_current}\n")
        # ---- POINT DATA ----
        f.write(f"POINT_DATA {npts}\n")
        # Scalar omega (1 value per point)
        f.write("SCALARS C double 1\n")
        f.write("LOOKUP_TABLE default\n")
        np.savetxt(f, phi_flat, fmt="%.6g")

        f.flush()
        os.fsync(f.fileno())
        f.close()
        os.replace(temp_filename,filename)
        if f.closed == False:
            print("file is not closed properly... \n")
        
        

