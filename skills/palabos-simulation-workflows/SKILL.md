---
name: palabos-simulation-workflows
description: This skill should be used when users ask about simulation workflows in palabos; it prioritizes documentation references and then source inspection only for unresolved details.
---

# palabos: Simulation Workflows

## High-Signal Playbook

### Route conditions
- Stay in this skill for full run pipelines (`bloodFlowDefoBodies`, `aneurysm`, `TRTdynamics`) and restart/convergence flow.
- Route first-time compile smoke tests to `palabos-getting-started`.
- Route GPU/MPI scaling bottlenecks and cluster tuning to `palabos-parallel-hpc`.
- Route generic showcase walkthroughs to `palabos-examples-and-tutorials`.

### Triage questions
- Which workflow: `bloodFlowDefoBodies`, `aneurysm`, or `TRTdynamics`?
- CPU or GPU execution path?
- For blood flow runs, is a valid `CPs` folder or `initialPlacing.pos` available?
- Are periodic directions split across multiple MPI subdomains when periodicity is enabled?
- For aneurysm, was `aneurysm.stl.tgz` extracted and `param.xml` checked (`maxLevel`, `epsilon`, inlet mode)?
- For TRT, which collision variant (`BGK`, `TRT`, `BGKma1`, `TRTma1`) and magic parameter are required?

### Canonical workflow
1. Pick the scenario and open its README + XML inputs (`examples/showCases/bloodFlowDefoBodies/README.md`, `examples/showCases/aneurysm/README.md`, `examples/codesByTopic/TRTdynamics/README.md`).
2. Build the example in its local `build/` directory (`mkdir -p build && cd build && cmake .. && make -j`).
3. For blood flow, run cell packing first to generate initial positions, or provide `initialPlacing.pos`.
4. Launch the main simulation command (CPU or GPU form) with consistent MPI/task settings.
5. For aneurysm, ensure `aneurysm.stl` is extracted from `aneurysm.stl.tgz` and run `./aneurysm param.xml`.
6. For TRT, set collision model/magic parameter in code, build, and run `./TRTdynamics`.
7. Monitor convergence/diagnostics (`ValueTracer`, runtime statistics, probes) and validate generated VTK/GIF outputs.
8. Escalate unresolved behavior to targeted source files: `examples/showCases/bloodFlowDefoBodies/bloodFlowDefoBodies.cpp`, `examples/showCases/aneurysm/aneurysm.cpp`, `examples/codesByTopic/TRTdynamics/TRTdynamics.cpp`, then `src/core/runTimeDiagnostics.h` and `src/complexDynamics/trtDynamics.h`.

### Minimal working example
```bash
cd examples/showCases/bloodFlowDefoBodies
cd build && cmake .. && make -j
cd ..
mpirun -n 4 ./bloodFlowDefoBodies cellPacking_params.xml
mpirun -n 4 ./bloodFlowDefoBodies shear_params.xml
```

```bash
cd examples/showCases/aneurysm
tar -xzvf aneurysm.stl.tgz
cd build && cmake .. && make -j
cd ..
mpirun -n 4 ./aneurysm param.xml
```

```bash
cd examples/codesByTopic/TRTdynamics
cd build && cmake .. && make -j
cd ..
./TRTdynamics
```

### Pitfalls and fixes
- Missing `CPs` or `initialPlacing.pos` causes startup failure/crash in blood-flow runs (`examples/showCases/bloodFlowDefoBodies/README.md`).
- Periodic directions without MPI subdivision can trigger erroneous stretched-body behavior; split periodic axes across MPI ranks (`examples/showCases/bloodFlowDefoBodies/README.md`).
- CUDA workstation mismatch (compute capability/architecture flags) breaks or degrades GPU runs; align CMake GPU arch and device constraints (`examples/showCases/bloodFlowDefoBodies/README.md`).
- Aneurysm run fails before start because geometry archive was not extracted (`examples/showCases/aneurysm/README.md`).
- Inlet/outlet assignment appears wrong when `param.xml` direction/type ordering is inconsistent with STL openings (`examples/showCases/aneurysm/README.md`).
- TRT behavior differs from expectation when `magicParameter` is not set consistently (`examples/codesByTopic/TRTdynamics/README.md`).

### Convergence/validation checks
- Use `ValueTracer` and multigrid criteria (`maxLevel`, `epsilon`) to confirm aneurysm convergence (`examples/showCases/aneurysm/README.md`).
- Confirm blood-cell initialization outputs exist after cell packing before launching long production runs.
- Validate probes/surface outputs in aneurysm (`pointMeasures()`, `writeSurfaceVTK()` paths from README concepts).
- Compare TRT variants on the same cavity setup to verify expected stability/accuracy changes.
- Check that CPU/GPU runs use equivalent physical inputs before comparing results.

## Scope
- Handle questions about simulation setup, execution flow, and runtime controls.
- Keep responses abstract and architectural for large codebases; avoid exhaustive per-function documentation unless requested.

## Primary documentation references
- `examples/showCases/bloodFlowDefoBodies/README.md`
- `examples/showCases/aneurysm/README.md`
- `examples/codesByTopic/TRTdynamics/README.md`

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
- `examples/showCases/bloodFlowDefoBodies/bloodFlowDefoBodies.cpp`
- `examples/showCases/bloodFlowDefoBodies/cellPacking_params.xml`
- `examples/showCases/bloodFlowDefoBodies/shear_params.xml`
- `examples/showCases/aneurysm/aneurysm.cpp`
- `examples/showCases/aneurysm/aneurysm_bounceback.cpp`
- `examples/showCases/aneurysm/param.xml`
- `examples/codesByTopic/TRTdynamics/TRTdynamics.cpp`
- `src/core/runTimeDiagnostics.h`
- `src/core/runTimeDiagnostics.cpp`
- `src/complexDynamics/trtDynamics.hh`
- `src/complexDynamics/trtDynamics.h`
- `src/core/dynamics.hh`
- `src/core/dynamics.h`
- `src/boundaryCondition/generalizedBoundaryDynamicsSolvers.hh`
- `src/boundaryCondition/generalizedBoundaryDynamicsSolvers.h`
- Prefer targeted source search (for example: `rg -n "<symbol_or_keyword>" coupledSimulators src utility`).
