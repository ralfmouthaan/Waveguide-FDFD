# Waveguide-FDFD

Finite difference frequency domain script for calculating progation-invariants supported by arbitrary geometry waveguides. Assumes a translationally-invariant geometry to express Maxwell's equations on a 2D Yee cell. The approach is fully described in my [paper](https://doi.org/10.1109/JLT.2021.3124469) and references therein. The user defines a problem geometry and passes this to ModeSolverFD.m for a full solution. 

Alternatively, the problem geometry can be passed to ModeSolverFD_LP, which calculates the linearly polarised modes supported by the waveguide. I have derived the maths for calculating the linearly polarised modes myself, and the details for these calculations are included in a pdf. I have not seen solutions that solve for only the linearly polarised modes before, but this more practically useful when exciting modes in waveguide with an SLM as I do. Furthermore, the problem size is reduced, and so the calculations require less memory and run faster. Even the vector mode solver runs faster than Lumerical - this is down to the careful use of sparse matrices and the strength of Matlab's eigensolver.

The scripts have been extensively benchmarked against well-defined problem geometries - see for example my [paper](https://doi.org/10.1109/JLT.2021.3124469) and especially my PhD thesis.

The outputs of these scripts can subsequently be used by Waveguide-Mode-Hologen to generate holograms for exciting these modes.
