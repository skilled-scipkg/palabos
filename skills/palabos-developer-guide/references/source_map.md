# palabos source map: Developer Guide

Use this map only after `dco/README.md` and the topic docs in `references/doc_map.md`.

## Fast source navigation
- `rg -n "<symbol_or_keyword>" src coupledSimulators utility`
- `rg -n "class|struct|namespace|template" src coupledSimulators utility`
- Start with the symbol named in docs/issues, then inspect the owning header and implementation pair.

## Function-level source entry points

### Core API and extension surfaces
- `src/palabos2D.h` | top-level API include for 2D simulations and extension discovery.
- `src/palabos3D.h` | top-level API include for 3D simulations and extension discovery.
- `src/core/dynamics.h` | `Dynamics::collide`, `Dynamics::setParameter`, `Dynamics::computeRhoBarJ`.
- `src/core/dynamics.hh` | default dynamics behavior and boundary dynamics implementations.
- `src/core/cell.h` | cell-level data access path used by custom dynamics.
- `src/core/cell.hh` | inline cell behavior used by dynamics and processors.
- `src/core/runTimeDiagnostics.h` | runtime diagnostics API surfaces for checks and instrumentation.
- `src/core/runTimeDiagnostics.cpp` | diagnostics implementation and reporting behavior.

### Boundary behavior extension points
- `src/boundaryCondition/generalizedBoundaryDynamicsSolvers.h` | generalized boundary solver interfaces.
- `src/boundaryCondition/generalizedBoundaryDynamicsSolvers.hh` | boundary solver implementation details.

### Coupled solver integration points
- `coupledSimulators/npFEM/Solver.h` | CPU solver API (`initialize`, `solve`, `collisionDetection`).
- `coupledSimulators/npFEM/Solver.cpp` | CPU solver execution path and convergence internals.
- `coupledSimulators/npFEM/shapeOpWrapper.h` | coupling boundary between Palabos and ShapeOp-derived logic.
- `coupledSimulators/npFEM/shapeOpWrapper.cpp` | wrapper implementation used by coupled simulations.
