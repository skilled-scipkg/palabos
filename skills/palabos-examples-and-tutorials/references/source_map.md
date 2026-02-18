# palabos source map: Examples and Tutorials

Use this map after the example READMEs in `references/doc_map.md`.

## Fast source navigation
- `rg -n "<symbol_or_keyword>" examples/showCases src`
- `rg -n "int main|writeVTK|writeGif|writeGifs|setBoundaryVelocity|collideAndStream" examples/showCases/*/*.cpp`

## Function-level source entry points

### Showcase driver files
- `examples/showCases/cavity3d/cavity3d.cpp` | `cavitySetup`, `writeGifs`, `writeVTK`, `main`.
- `examples/showCases/rectangularChannelWithCylinder3d/rectangularChannelWithCylinder3d.cpp` | `poiseuilleVelocity`, `writeGifs`, `writeVTK`, `main`.
- `examples/showCases/movingWall/movingWall.cpp` | `setParam`, `computeLbParam`, `writeVTK`, immersed-boundary loop in `main`.
- `examples/showCases/settlingDrivenConvection/settlingDrivenConvection.cpp` | `ExpSetup`, `writeVTK`, `writeGif`, `main`.
- `examples/showCases/settlingDrivenConvection/SDCfunctions.h` | advection-diffusion helper functions used by settling-driven convection.
- `examples/showCases/partialBounceBack/fibrinolysis.cpp` | clot-loading flow setup and `main` run path.

### Shared core paths for behavior verification
- `src/palabos3D.h` | main 3D API include reached by showcase drivers.
- `src/palabos3D.hh` | inline 3D implementation details.
- `src/palabos2D.h` | fallback for 2D utility functions reused across examples.
- `src/core/dynamics.h` | dynamics API and collision parameter entry points.
- `src/boundaryCondition/generalizedBoundaryDynamicsSolvers.h` | generalized boundary behavior.
