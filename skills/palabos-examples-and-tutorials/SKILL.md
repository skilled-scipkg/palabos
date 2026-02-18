---
name: palabos-examples-and-tutorials
description: This skill should be used when users ask about examples and tutorials in palabos; it prioritizes documentation references and then source inspection only for unresolved details.
---

# palabos: Examples and Tutorials

## High-Signal Playbook

### Route conditions
- Stay in this skill for walkthrough-style usage of showcase/tutorial examples.
- Route first-install verification or starter collision exercises to `palabos-getting-started`.
- Route multi-step blood-flow/aneurysm/TRT run pipelines to `palabos-simulation-workflows`.
- Route MPI/GPU scaling and cluster execution details to `palabos-parallel-hpc`.

### Triage questions
- Which example is targeted: `cavity3d`, `rectangularChannelWithCylinder3d`, `movingWall`, `settlingDrivenConvection`, or `partialBounceBack`?
- Is the user asking for run instructions, parameter interpretation, or code-structure mapping?
- Is the run single-process or MPI (`mpirun -np X ...`)?
- Are geometry/physics arguments expected in physical units or lattice units?
- Is the focus on post-processing outputs (GIF/VTK) or on numerical behavior?
- Are immersed-boundary motion parameters (`ampl`, `freq`) being modified?

### Canonical workflow
1. Open the specific README under `examples/showCases/*` and identify compile/run requirements.
2. Build in local `build/` with CMake and `make`.
3. Run the executable with documented command form (for `settlingDrivenConvection`, use MPI launch).
4. Confirm required runtime inputs before launch (for `partialBounceBack`, verify `clot_0.0003.txt` is present).
5. Inspect `tmp/` outputs (`.gif`, `.vtk`) to verify run success.
6. Map README concepts to example source (`*.cpp`, helper headers) before escalating to `src/`.
7. Escalate to core API entry points (`src/palabos2D.h`, `src/palabos3D.h`) only when example-level files are insufficient.

### Minimal working example
```bash
cd examples/showCases/settlingDrivenConvection
cd build && cmake .. && make -j
cd ..
mpirun -np 4 ./settlingDrivenConvection
```

```bash
cd examples/showCases/cavity3d
cd build && cmake .. && make -j
cd ..
./cavity3d
```

```bash
cd examples/showCases/partialBounceBack
mkdir -p build && cd build
cmake .. && make -j
cd ..
mpirun -n 4 ./fibrinolysis
```

### Pitfalls and fixes
- `rectangularChannelWithCylinder3d` geometry CLI arguments are physical units and internally non-dimensionalized; wrong assumptions cause scale errors (`examples/showCases/rectangularChannelWithCylinder3d/README.md`).
- `movingWall` behavior looks wrong when `ampl`/`freq` are changed inconsistently (`examples/showCases/movingWall/Readme.md`).
- `settlingDrivenConvection` at higher resolution becomes resource-heavy; start from documented moderate resolution (`examples/showCases/settlingDrivenConvection/README.md`).
- Users expect dedicated compile instructions in every showcase; some examples rely on the standard CMake pattern.
- `clot_0.0003.txt` is input data (solid-fraction map), not an executable script (`examples/showCases/partialBounceBack/README.md`).

### Convergence/validation checks
- Confirm each example produces expected artifact types in `tmp/` (GIF/VTK presence and update cadence).
- For cavity/channel examples, check stable average fields and physically plausible vortex structure.
- For settling-driven convection, compare qualitative instability evolution across resolutions.
- For moving-wall, verify periodic plate kinematics over time against configured oscillation parameters.
- For partial bounce-back, verify clot map ingestion and stable pressure-driven flow without startup abort.
- When comparing examples, keep discretization and boundary assumptions explicit to avoid false cross-case conclusions.

## Scope
- Handle questions about worked examples, tutorials, and cookbook usage.
- Keep responses abstract and architectural for large codebases; avoid exhaustive per-function documentation unless requested.

## Primary documentation references
- `examples/showCases/cavity3d/README.md`
- `examples/showCases/rectangularChannelWithCylinder3d/README.md`
- `examples/showCases/movingWall/Readme.md`
- `examples/showCases/settlingDrivenConvection/README.md`
- `examples/showCases/partialBounceBack/README.md`

## Required runtime inputs
- `examples/showCases/partialBounceBack/clot_0.0003.txt`

## Workflow
- Start with the primary references above.
- If details are missing, inspect `references/doc_map.md` for the complete topic document list.
- Use tutorials/examples as executable usage patterns when available.
- Use tests as behavior or regression references when available.
- If ambiguity remains after docs, inspect `references/source_map.md` and start with the ranked source entry points.
- Cite exact documentation file paths in responses.

## Tutorials and examples
- `examples/benchmarks`
- `examples/codesByTopic`
- `examples/gpuExamples`
- `examples/showCases`
- `examples/tutorial`

## Test references
- `examples/codesByTopic`
- `examples/showCases`

## Optional deeper inspection
- `coupledSimulators`
- `src`
- `utility`

## Source entry points for unresolved issues
- `examples/showCases/cavity3d/cavity3d.cpp`
- `examples/showCases/rectangularChannelWithCylinder3d/rectangularChannelWithCylinder3d.cpp`
- `examples/showCases/movingWall/movingWall.cpp`
- `examples/showCases/settlingDrivenConvection/settlingDrivenConvection.cpp`
- `examples/showCases/settlingDrivenConvection/SDCfunctions.h`
- `examples/showCases/partialBounceBack/fibrinolysis.cpp`
- `src/palabos3D.hh`
- `src/palabos3D.h`
- `src/palabos2D.hh`
- `src/palabos2D.h`
- `src/boundaryCondition/generalizedBoundaryDynamicsSolvers.hh`
- `src/boundaryCondition/generalizedBoundaryDynamicsSolvers.h`
- `src/core/dynamics.hh`
- `src/core/dynamics.h`
- Prefer targeted source search (for example: `rg -n "<symbol_or_keyword>" examples/showCases src`).
