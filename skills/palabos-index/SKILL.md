---
name: palabos-index
description: This skill should be used when users ask how to use palabos and the correct generated documentation skill must be selected before going deeper into source code.
---

# palabos Skills Index

## Route the request
- Classify the request into one of the generated topic skills listed below.
- Prefer abstract, workflow-level guidance for large scientific packages; do not attempt full function-by-function coverage unless explicitly requested.
- For immediate runnable commands from repo root, read `references/simulation_bootstrap.md` first.

## Generated topic skills
- `palabos-parallel-hpc`: Parallel and HPC (MPI/OpenMP/GPU execution, scaling, and batch systems)
- `palabos-examples-and-tutorials`: Examples and Tutorials (worked examples, tutorials, and cookbook usage)
- `palabos-simulation-workflows`: Simulation Workflows (simulation setup, execution flow, and runtime controls)
- `palabos-getting-started`: Getting Started (initial setup, quickstarts, and core concepts)
- `palabos-developer-guide`: Developer Guide (developer architecture, extension points, and contribution workflow)

## Documentation-first inputs
- `dco`
- `examples`
- `coupledSimulators/npFEM/npFEM_StandAlone_RhinoGH`

## Tutorials and examples roots
- `examples/benchmarks`
- `examples/codesByTopic`
- `examples/gpuExamples`
- `examples/showCases`
- `examples/tutorial`

## Test roots for behavior checks
- `examples/codesByTopic`
- `examples/showCases`

## Escalate only when needed
- Start from topic skill primary references.
- If those references are insufficient, search the topic skill `references/doc_map.md`.
- If documentation still leaves ambiguity, open `references/source_map.md` inside the same topic skill and inspect the suggested source entry points.
- Use targeted symbol search while inspecting source (e.g., `rg -n "<symbol_or_keyword>" coupledSimulators src utility`).

## Source directories for deeper inspection
- `coupledSimulators`
- `src`
- `utility`
