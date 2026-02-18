---
name: palabos-getting-started
description: This skill should be used when users ask about getting started in palabos; it prioritizes documentation references and then source inspection only for unresolved details.
---

# palabos: Getting Started

## High-Signal Playbook

### Route conditions
- Stay in this skill for first compile/run checks and starter cases in `examples/showCases/cavity2d` and `examples/showCases/dsl2d`.
- Route deformable-bodies and case-study pipelines to `palabos-simulation-workflows`.
- Route MPI/GPU scaling and cluster performance to `palabos-parallel-hpc`.
- Route contribution/legal or architecture-extension requests to `palabos-developer-guide`.

### Triage questions
- Which starter case is needed: `cavity2d` or `dsl2d`?
- Is this an installation smoke test, a stability study, or a performance run?
- Will execution be serial (`./cavity2d`) or MPI (`mpirun -np ...`)?
- For `dsl2d`, what Reynolds number, Mach number, and grid size are intended?
- Do collision settings need edits in `examples/showCases/dsl2d/config.xml` (`dynName`, `hoOmega`, `bulk`)?
- Should output files be generated for visualization, or reduced for timing?

### Canonical workflow
1. Build the selected example inside its local `build` directory (`examples/showCases/cavity2d/README.md`, `examples/showCases/dsl2d/README.md`).
2. Run baseline executable (`./cavity2d` or `mpirun -np 2 ./dsl2d config.xml 10000 128 0.1`).
3. Confirm health via console output and files in the example `tmp/` directory.
4. For DSL stability sweeps, vary `dynName`, `hoOmega`, and `bulk` exactly as guided by exercises in `examples/showCases/dsl2d/README.md`.
5. For benchmarking, reduce output overhead per README guidance (disable frequent GIF/VTK writes).
6. Escalate to source only if behavior is unclear: `examples/showCases/cavity2d/cavity2d.cpp`, `examples/showCases/dsl2d/dsl2d.cpp`, `examples/showCases/dsl2d/utility_dsl2d_param.h`.

### Minimal working example
```bash
cd examples/showCases/cavity2d
cd build && cmake .. && make -j
cd ..
./cavity2d
```

```bash
cd examples/showCases/dsl2d
cd build && cmake .. && make -j
cd ..
mpirun -np 2 ./dsl2d config.xml 10000 128 0.1
```

### Pitfalls and fixes
- `dsl2d` becomes unstable under under-resolved/high-Re settings: apply README collision-model guidance (`bulk=true`, `hoOmega=REG`, or RR/K variants) before arbitrary parameter changes (`examples/showCases/dsl2d/README.md`).
- Performance comparisons are skewed by output cost: suppress frequent GIF/VTK output as described in `examples/showCases/dsl2d/README.md`.
- Boundary-condition expectations mismatch implementation: verify actual boundary setup in `examples/showCases/cavity2d/cavity2d.cpp`.
- Missing expected artifacts: outputs are written under each example `tmp/` directory.
- Build/run path confusion: compile in local `build/`, then execute from the parent example directory.

### Convergence/validation checks
- For `cavity2d`, verify stable average density/energy logging and expected output cadence (`examples/showCases/cavity2d/README.md`).
- For `dsl2d`, treat spurious-vortex growth/blow-up as a resolution/collision-model warning (`examples/showCases/dsl2d/README.md`).
- Re-run with identical `config.xml` and CLI values; outputs/log trends should be reproducible.
- Increase grid resolution at fixed physics before concluding a model-quality issue.
- Ensure run mode matches intent: visualization mode keeps outputs; benchmark mode minimizes them.

## Scope
- Handle questions about initial setup, quickstarts, and core concepts.
- Keep responses abstract and architectural for large codebases; avoid exhaustive per-function documentation unless requested.

## Primary documentation references
- `examples/showCases/dsl2d/README.md`
- `examples/showCases/cavity2d/README.md`

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
- `examples/showCases/cavity2d/cavity2d.cpp`
- `examples/showCases/dsl2d/dsl2d.cpp`
- `examples/showCases/dsl2d/utility_dsl2d_param.h`
- `src/palabos2D.hh`
- `src/palabos2D.h`
- `src/palabos3D.hh`
- `src/palabos3D.h`
- `src/boundaryCondition/generalizedBoundaryDynamicsSolvers.hh`
- `src/boundaryCondition/generalizedBoundaryDynamicsSolvers.h`
- `src/core/dynamics.hh`
- `src/core/dynamics.h`
- Prefer targeted source search (for example: `rg -n "<symbol_or_keyword>" examples/showCases src`).
