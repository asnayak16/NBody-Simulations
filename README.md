# NBody-Simulations
Nbody Simulations done for UCSD's PHYS 241 Course
Last Update 22-Mar-2017

plummer_run.f90 is a FORTRAN90 Code written to perform N-Body simulations and is optimised to use threading using OpenMP directives
Compiled using Intel's Fortran Compiler : "ifort /Qopenmp" for best performance. 
May give errors if OMP_STACKSIZE env variable isn't set in a few MBs.

nbody_gpu.cu is a CUDA code to perform the same N-Body simulations. Written to run on the NVIDIA GTX 960M Graphics card. Compiled using nvcc and VC compiler on visual Studio. -gencode=arch=compute_50,code=\"sm_50,compute_50\ used.
FUTURE WORK: Needs update to double precision.

Try with 'pl0.dat' initial condition for stable plummer simulation.
Inform for any changes needed.
