# Evidence: palabos-simulation-workflows

## Primary docs
- `examples/showCases/bloodFlowDefoBodies/README.md`
- `examples/showCases/aneurysm/README.md`
- `examples/codesByTopic/TRTdynamics/README.md`

## Primary source entry points
- `skills/palabos-simulation-workflows/references/doc_map.md`
- `src/core/runTimeDiagnostics.h`
- `src/core/runTimeDiagnostics.cpp`
- `src/complexDynamics/trtDynamics.hh`
- `src/complexDynamics/trtDynamics.h`
- `src/core/dynamics.hh`
- `src/core/dynamics.h`
- `src/boundaryCondition/generalizedBoundaryDynamicsSolvers.hh`
- `src/boundaryCondition/generalizedBoundaryDynamicsSolvers.h`
- `src/core/plbTimer.h`
- `src/core/plbTimer.cpp`
- `src/algorithm/timePeriodicSignal.hh`
- `src/algorithm/timePeriodicSignal.h`
- `src/latticeBoltzmann/externalFields.h`
- `src/complexDynamics/asinariModel.hh`
- `src/complexDynamics/asinariModel.h`
- `src/complexDynamics/advectionDiffusionBoundaryInstantiator3D.h`
- `src/complexDynamics/advectionDiffusionBoundaryInstantiator2D.h`
- `src/complexDynamics/advectionDiffusionBoundaryCondition3D.hh`
- `src/complexDynamics/advectionDiffusionBoundaryCondition3D.h`

## Extracted headings
- Instructions for Palabos-npFEM
- Compilation
- Important note on CUDA-capable Workstations
- Important note on Periodicity and MPI
- Perform Cell Packing
- Cellular Blood Flow Simulations
- Case study: Collision at an obstacle
- Flow inside a synthetic 3D aneurysm
- Concepts used in this example
- Before running
- Use TRT dynamics
- Lid driven cavity in 2D

## Executable command hints
- mpirun -n X ./bloodFlowDefoBodies cellPacking_params.xml
- mpirun -n X ./bloodFlowDefoBodies_gpu cellPacking_params.xml NumberOfNodes NumberOfGPUsPerNode
- mpirun -n X ./bloodFlowDefoBodies shear_params/poiseuille_params.xml
- mpirun -n X ./bloodFlowDefoBodies_gpu shear_params/poiseuille_params.xml NumberOfNodes NumberOfGPUsPerNode

## Warnings and pitfalls
- ## Important note on CUDA-capable Workstations
- ## Important note on Periodicity and MPI
- ***It is very important to understand that without the right CPs folder (generated through Cell Packing or provided) or initialPlacing.pos file, the Cellular Blood Flow Simuation is not going to start (it crashes).***
- by a factor of two at each intermediate convergence (the number
- of refinements and the convergence criterion are the
- * Use of `ValueTracer` to check for convergence.
