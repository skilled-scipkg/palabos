# palabos source map: Getting Started

Use this map after the starter docs in `references/doc_map.md`.

## Fast source navigation
- `rg -n "<symbol_or_keyword>" examples/showCases src`
- `rg -n "int main|cavitySetup|simulationSetup|getDynamics|writeVTK|writeGif" examples/showCases/cavity2d/cavity2d.cpp examples/showCases/dsl2d/dsl2d.cpp`

## Function-level source entry points

### Starter case drivers
- `examples/showCases/cavity2d/cavity2d.cpp` | `cavitySetup`, `writeGif`, `writeVTK`, `main`.
- `examples/showCases/dsl2d/dsl2d.cpp` | `getDynamics`, `simulationSetup`, `writeGif`, `writeVTK`, `main`.
- `examples/showCases/dsl2d/utility_dsl2d_param.h` | `Param` XML/CLI parsing, `writeLogFile`, collision-model parameter wiring.
- `examples/showCases/dsl2d/config.xml` | starter input for `dynName`, `hoOmega`, `bulk`, and output cadence.

### Core behavior checks (if starter files are insufficient)
- `src/palabos2D.h` | core 2D API access used by starter examples.
- `src/palabos2D.hh` | 2D inline implementation surfaces reached from starter code.
- `src/core/dynamics.h` | dynamics interface and runtime parameter hooks.
- `src/core/dynamics.hh` | default collide/compute behavior used by selected dynamics.
- `src/boundaryCondition/generalizedBoundaryDynamicsSolvers.h` | generalized boundary behavior references.
- `src/boundaryCondition/generalizedBoundaryDynamicsSolvers.hh` | boundary solver implementations.
