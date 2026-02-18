# Evidence: palabos-getting-started

## Primary docs
- `examples/showCases/dsl2d/README.md`
- `examples/showCases/cavity2d/README.md`

## Primary source entry points
- `skills/palabos-getting-started/references/doc_map.md`
- `src/palabos3D.hh`
- `src/palabos3D.h`
- `src/palabos2D.hh`
- `src/palabos2D.h`
- `src/boundaryCondition/generalizedBoundaryDynamicsSolvers.hh`
- `src/boundaryCondition/generalizedBoundaryDynamicsSolvers.h`
- `src/particles/particleField3D.hh`
- `src/particles/particleField3D.h`
- `src/particles/particleField2D.hh`
- `src/particles/particleField2D.h`
- `src/particles/multiParticleField3D.hh`
- `src/particles/multiParticleField3D.h`
- `src/particles/multiParticleField2D.hh`
- `src/particles/multiParticleField2D.h`
- `src/parallelism/parallelMultiDataField3D.hh`
- `src/parallelism/parallelMultiDataField3D.h`
- `src/parallelism/parallelMultiDataField2D.hh`
- `src/parallelism/parallelMultiDataField2D.h`
- `src/offLattice/triangleBoundary3D.hh`

## Extracted headings
- Stability of collision models: Double shear layer
- Introduction
- Guide
- Exercice 1: Test the compilation and run the code in parallel
- Explanations: Code structure
- Exercice 2: Improving stability without modifying relaxation parameters
- Exercice 3: Improving stability and accuracy by increasing the bulk viscosity
- Exercice 4: Further improving stability and accuracy by regularizing high-order moments
- Exercice 5: Improving stability and accuracy without sacrifying acoustic waves
- Take home ideas
- Lid driven cavity in 2D
- Theory

## Executable command hints
- mpirun -np 2 ./dsl2d config.xml 10000 128 0.1
- $ cd build
- $ cmake ..
- $ make
- $ cd ..
- $ ./cavity2d

## Warnings and pitfalls
- # Stability of collision models: Double shear layer
- With this idea in mind, a 2D double shear layer will be simulated in under-resolved conditions since it allows to highlight dispersion and stability issues through very intuitive reasonings [1,2]. This test is carried out in a fully periodic domain without boundary conditions, and it consists in two shear layers evolving over time. By imposing a given (transverse velocity) perturbation at the initial time, these layers roll-up and form two counter-rotating vortices. Interestingly, when the Reynolds number is increased, while keeping the mesh resolution constant, one can observe the formation of secondary spurious vortices. The latter help us vizualizing dispersion issues that usually lead to the simulation blow-up, e.g., when the Reynolds or Mach numbers are increased, or if the mesh resolution if further decreased.
- ### Exercice 2: Improving stability without modifying relaxation parameters
- ### Exercice 3: Improving stability and accuracy by increasing the bulk viscosity
- ### Exercice 4: Further improving stability and accuracy by regularizing high-order moments
- ### Exercice 5: Improving stability and accuracy without sacrifying acoustic waves
- For pure aerodynamics, regularizing high-order moments and increasing the bulk viscosity are good ways to improve the stability of most collision models. Only the RR approach should be used in an SRT formalism. In case of computational aeroacoustic simulations, K-REG, CHM-REG and RR-SRT should be preferred.
- [5] Coreixas, C., Wissocq, G., Chopard, B., & Latt, J. (2020). Impact of collision models on the physical properties and the stability of lattice Boltzmann methods. Philosophical Transactions of the Royal Society A, 378(2175), 20190397. https://www.researchgate.net/publication/339136823_Impact_of_collision_models_on_the_physical_properties_and_the_stability_of_lattice_Boltzmann_methods
- The simulation domain/lattice is created with an instance of the `MultiBlockLattice2D` class. The arguments include the dynamic model, in this case, the BGK collision model (`BGKdynamics()`). It is also important to specify the model descriptor, in this case, the `D2Q9Descriptor` (line 125-127). Note that `Nx` and `Ny` are the number of cells in $`x`$ and $`y`$ direction and thus the size of the simulation domain.
