# palabos source map: Simulation Workflows

Use this map after the workflow docs in `references/doc_map.md`.

## Fast source navigation
- `rg -n "<symbol_or_keyword>" examples/showCases examples/codesByTopic src`
- `rg -n "int main|readParameters|setOpenings|pointMeasures|selectDynamics|writeVTK|collideAndStream" examples/showCases/bloodFlowDefoBodies/bloodFlowDefoBodies.cpp examples/showCases/aneurysm/aneurysm.cpp examples/codesByTopic/TRTdynamics/TRTdynamics.cpp`

## Function-level source entry points

### Workflow drivers and inputs
- `examples/showCases/bloodFlowDefoBodies/bloodFlowDefoBodies.cpp` | `main`, cell-packing checks, `writeVTK`, flow loop.
- `examples/showCases/bloodFlowDefoBodies/cellPacking_params.xml` | starter cell-packing run configuration.
- `examples/showCases/bloodFlowDefoBodies/shear_params.xml` | shear-flow runtime configuration.
- `examples/showCases/bloodFlowDefoBodies/poiseuille_params.xml` | Poiseuille runtime configuration.
- `examples/showCases/bloodFlowDefoBodies/initialPlacing.pos` | fallback initial placement input.
- `examples/showCases/aneurysm/aneurysm.cpp` | `readParameters`, `setOpenings`, `pointMeasures`, multigrid run loop.
- `examples/showCases/aneurysm/aneurysm_bounceback.cpp` | bounce-back boundary variant for comparison.
- `examples/showCases/aneurysm/param.xml` | aneurysm geometry/opening/convergence settings.
- `examples/codesByTopic/TRTdynamics/TRTdynamics.cpp` | `selectDynamics`, `cavitySetup`, `writeVTK`, TRT variant switching.

### Core dynamics and diagnostics used by workflows
- `src/complexDynamics/trtDynamics.h` | TRT class interfaces and magic-parameter behavior.
- `src/complexDynamics/trtDynamics.hh` | TRT implementation details.
- `src/core/runTimeDiagnostics.h` | diagnostics API used for runtime checks.
- `src/core/runTimeDiagnostics.cpp` | diagnostics implementation.
- `src/core/dynamics.h` | common collision/dynamics base interfaces.
- `src/core/dynamics.hh` | common dynamics implementations.
- `src/core/plbTimer.h` | timer utilities for runtime/performance checkpoints.
