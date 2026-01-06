# 2D Advection Diffusion Solver 
Finite difference method (FDM) solver for 2D advection diffusion problems on a structured grid.

## Demo (GIF)
![Demo GIF](media/demo.gif)

## Features
- Structured grid
- Boundary conditions handled with ghost cell application 
- Adaptive time step calculation 
- Numba parallelization of the stencil loop
- VTK output for ParaView

## Requirements
- Python 3.9+
- numpy
- numba
- ParaView (for output visualization) 