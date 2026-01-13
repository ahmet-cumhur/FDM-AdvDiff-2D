# Validation Of Numeric Solution Of 1D Slab Case With Analytical Solution 

## Purpose
valiadtion of the solvers accuracy + stability at highly diffusion dominated cases

## Governing Equation
$$
\frac{\partial C}{\partial t} + \mathbf{u}\cdot\nabla C = \nu \nabla^2 C
$$

<div align= "center" > 
in vectoral form
</div>

In 2D with derivative form:

$$
\frac{\partial C}{\partial t} + u*\frac{\partial C}{\partial x} = \nu*\frac{\partial^2 C}{\partial x^2}
$$

<div align= "center" > 
for x direction
</div>

$$
\frac{\partial C}{\partial t} + v*\frac{\partial C}{\partial y} = \nu*\frac{\partial^2 C}{\partial y^2}
$$ 

<div align= "center" > 
for y direction
</div>

we use the upwind scheme (US) and central differences scheme (CDS) for discretize these equations. 
For time discretation we apply euler explicit scheme (EES) because of its convenciece and easy apply. In future implicit schemes will be added too.

Time discretezation:
$$ 
\frac{\partical C}{\partial t} = \frac{C^{n+1}-C^{n}}{\Delta t}
$$

Spatial discretezation for x direction:

$$
u*\frac{\partial C}{\partial x} = \frac{C_{i,j}-C_{i-1,j}}{\Delta x}
$$

<div align="center">
US convective term 
</div>

$$
\nu*\frac{\partial^2 C}{\partial x^2} = \frac{C_{i-1,j}-2*C_{i,j}+C_{i+1,j}}{\Delta x^2}
$$

<div align="center">
CDS diffison term 
</div>

Spatial discretezation for y direction:

$$
v*\frac{\partial C}{\partial y} = \frac{C_{i,j}-C_{i,j-1}}{\Delta y}
$$

<div align="center">
US convective term 
</div>

$$
\nu*\frac{\partial^2 C}{\partial y^2} = \frac{C_{i,j-1}-2*C_{i,j}+C_{i,j+1}}{\Delta y^2}
$$

<div align="center">
CDS diffison term 
</div>

So total equation with time and spatial discretazation becomes:

$$
{C_{i,j}^{n+1}} = {C_{i,j}^{n}}-{\Delta t}*(\frac{C_{i,j}^{n}-C_{i-1,j}^{n}}{\Delta x} + \frac{C_{i,j}^{n}-C_{i,j-1}^{n}}{\Delta y})+{\Delta t}*(\frac{C_{i-1,j}^{n}-2*C_{i,j}^{n}+C_{i+1,j}^{n}}{\Delta x^2} + \frac{C_{i,j-1}^{n}-2*C_{i,j}^{n}+C_{i,j+1}^{n}}{\Delta y^2})
$$

TBC...