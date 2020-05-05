***************************************************************
	INFORMATION FOR MVP:PROBLEM SHEET 3 SUBMISSION
***************************************************************

Declan Mathews
s1610357
30/03/2020

This submission should contain two folders: programmes and data.


########################   Programmes   ########################


This folder should contain the following three programmes:

1. Cahn-Hilliard_solver.py
2. Jacobi_solver.py
3. Gauss-Seidel_solver.py

Each of these contains a class for the solver by the same name.
Each solver is tailored towards the desired Poisson's equation

They all follow a similar structure with the sections of methods separated with
various comment styles to follow it easier.

They all contain comments that should hopefully help enough, but more information
will be included here to be referred to if needs be.


******* Loop/Convolution Information *******


For the Cahn-Hilliard (CH) and Jacobi programmes there are both loop methods and 
convolution methods. The data calculation/gathering and animation is done using the convolution
methods for time-saving as there is no difference in the results. However, the loop
methods are included for completion and can be used by commenting them in and commenting out the 
convolution methods in the run_dynamics or run_convergence methods.

For Gauss-Seidel (GS) there is no convolution option.


# --- Cahn-Hilliard_solver.py --- #


Use: python Cahn-Hilliard_solver.py N nsweeps dx dt phi_0

where:
N is the system size NxN
nsweeps is the desired number of sweeps
dx is the position discretisation
dt is the time discretisation
phi_0 is the lattice initialisation value (there is standard noise of +/- 0.1 on top of this)

General Notes:

- The noise level can be controlled by the noise_dampening value
- M and kappa are also available to change if desired
- There is a Data Analysis section in main that was used to gather data for the free energy plots


# --- Jacobi_solver.py --- #


Use: python Jacobi_solver.py N c_lim phi_0 type view

where:
N is the system size NxNxN
c_lim is the convergence limit (done by magnitude of vector distance from previous state)
phi_0 is the lattice initialisation value (there is standard noise of +/- 0.1 on top of this)
type is monopole or wire (central positive charge or central positive wire in z direction)
view is anim or conv (animation of moving towards convergence or just displays converged image)

General Notes:

- The noise level can be controlled by the noise_dampening value


# --- Gauss-Seidel_solver.py --- #


Use: python Gauss-Seidel_solver.py N c_lim phi_0 omega type view

where:
N is the system size NxNxN
c_lim is the convergence limit (done by magnitude of vector distance from previous state)
phi_0 is the lattice initialisation value (there is standard noise of +/- 0.1 on top of this)
omega is the SOR over-relaxation variable (=1 for standard GS)
type is monopole or wire (central positive charge or central positive wire in z direction)
view is anim or conv (animation of moving towards convergence or just displays converged image)

General Notes:

- The noise level can be controlled by the noise_dampening value
- There is an extra view option of gather-omega used for gathering results for the number of 
  steps to converge for various omega values (1 to 2)
- There are methods for doing SOR-GS in 2D included for gathering the data of number of sweeps
  for convergence v omega in order to save time


#########################   Data   ##########################


All images/data are saved and should be named as to be easily identified.
All data is for N = 100, however the SOR data is from a 2D lattice.
The monopole and wire data was gathered using Jacobi convolution for c_lim = 0.001.
The SOR data was done in steps of 0.01 for omega between 1 and 2 with c_lim = 0.01.


# --- Cahn-Hilliard Results --- #


Gathered with N = 100, nsweeps = 100,000, dx = 1.0, dt = 1.0

Number of sweeps v free energy graphs:

Data format: Num_sweeps Free_energy

- phi_0 = 0.0

  data: 	CH-FE-phi0.0-data100.txt
  graph: 	CH-FE-phi0.0-100.pdf

- phi_0 = 0.5

  data: 	CH-FE-phi0.5-data100.txt
  graph:	CH-FE-phi0.5-100.pdf


# --- Monopole Results --- #


Gathered with Jacobi convolution methods.
Settings: N = 100, c_lim = 0.001, phi_0 = 0.5

- Potential contour: mid-plane

  data:		monopole_contour_data100.txt
  format:	x y z potential  # Note that there are gaps every new x value for pm3d use in gnuplot
  graph:	monopole-potential-100.pdf

- E-field vector plot: mid-plane (plot is normalised in gun-lot, data is not normalised)

  data:		monopole_vector_data100.txt
  format:	x y z Ex Ey Ez
  graph:	monopole-vector-100.pdf

- Distance v Potential: mid-plane

  format of file:   distance potential field_magnitude

  data: 	monopole_dist_data100.txt columns 1:2 
  graph: 	monopole-loglog-potential.pdf
  fit:		fit.log ; second result: gradient fit = -0.86 between 0 and 2

- Distance v E-field magnitude: mid-plane

  data: 	monopole_dist_data100.txt columns 1:3 
  graph: 	monopole-loglog-Emag.pdf
  fit:		fit.log ; first result: gradient fit = -2.07 between 0 and 2


# --- Wire Results --- #


Gathered with Jacobi convolution methods.
Settings: N = 100, c_lim = 0.001, phi_0 = 0.5

- Potential contour: mid-plane

  data:		wire_contour_data100.txt
  format:	x y z potential  # Note that there are gaps every new x value for pm3d use in gnuplot
  graph:	wire-potential-100.pdf

- B-field vector plot: mid-plane (plot is normalised in gun-lot, data is not normalised)

  data:		wire_vector_data100.txt
  format:	x y z Bx By Bz
  graph:	wire-vector-100.pdf

- Distance v Potential: mid-plane

  format of file:   distance potential field_magnitude

  data: 	wire_dist_data100.txt columns 1:2 
  graph-1: 	wire-loglog-potential.pdf
  fit-1:	fit.log ; third result: gradient fit = -0.40 between 0 and 2
  graph-2: 	wire-loglinear-potential.pdf
  fit-2:	fit.log ; fifth result: gradient fit = -0.16 between 0 and 2

- Distance v B-field magnitude: mid-plane

  data: 	wire_dist_data100.txt columns 1:3 
  graph: 	wire-loglog-Bmag.pdf
  fit:		fit.log ; fourth result: gradient fit = -1.02 between 0 and 2


# --- SOR Results --- #


Gathered with SOR-GS using a 2D lattice to save time.
Settings: N = 100, c_lim = 0.01, phi_0 = 0.5 monopole

- Omega v number of sweeps until convergence: Omega from 1 to 2 in steps of 0.01

  data: 	SOR_data100.txt
  format:	omega sweeps_to_converge
  graph:	SOR-omega-100.pdf











