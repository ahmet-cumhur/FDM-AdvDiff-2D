import time
from adv_2d import Adv_Diff_Solver
from get_phi_parr import get_phi,get_phi_vectorized,apply_bc_n,apply_bc_s
from numba import set_num_threads, get_num_threads

#set the number of threads here
set_num_threads(8)
#define the object here
adv_solver = Adv_Diff_Solver(400,200,2.0,1.0,0.01,0.5,0.1,0.0)
#initiate the variables here
nx = adv_solver.nx
ny = adv_solver.ny
dx = adv_solver.dx
dy = adv_solver.dy
l = adv_solver.l
h = adv_solver.h
nu = adv_solver.nu
u = adv_solver.u
v = adv_solver.v
phi_inlet = adv_solver.phi_inlet
phi_n = adv_solver.phi_n
phi_s = adv_solver.phi_s


i = 0
t_current = 0.0
t_final = 2.0
t0 = time.perf_counter()
#apply the inital condition
adv_solver.apply_ic(phi_n,phi_s,nx,ny)
#main time loop
while t_current<t_final:
    i+=1
    dt = adv_solver.calculate_time_step_2(dx,dy,u,v,nu)
    #dt = 1e-5
    t_current += dt
    #solver
    apply_bc_n(phi_n,phi_inlet,nx,ny)
    get_phi(ny,nx,u,v,phi_n,phi_s,dt,dx,dy,nu,phi_inlet)
    apply_bc_s(phi_s,phi_inlet,nx,ny)
    adv_solver.time_roll(phi_n,phi_s)
    #write the data
    if i%200 == 0:
        adv_solver.vtk_converter(nx,ny,phi_n,dx,dy,t_current,i)
    if i%200 == 0:
        print(f"Time step size: {dt}\n")
        print(f"Flow time: {t_current}\n")
        print("Numba threads:", get_num_threads())

t1 = time.perf_counter()
st = "Total simulation time: "+str(adv_solver.time_counter(t0,t1))
print(st)