---
name: palabos-parallel-hpc
description: This skill should be used when users ask about parallel and hpc in palabos; it prioritizes documentation references and then source inspection only for unresolved details.
---

# palabos: Parallel and HPC

## High-Signal Playbook

### Route conditions
- Stay in this skill for MPI/GPU execution, scaling, compiler/toolchain selection, and cluster-oriented run patterns.
- Route full biological-flow pipeline logic (cell packing, CPs dependencies, case setup) to `palabos-simulation-workflows`.
- Route first local compile/run checks to `palabos-getting-started`.
- Route showcase modeling questions (non-HPC) to `palabos-examples-and-tutorials`.

### Triage questions
- Which parallel target: `partialBounceBack`, `sandstone`, `multiComponentPorous`, or `npFEM`?
- CPU-only MPI, GPU single node, or multi-node GPU?
- For sandstone cases, is `Berea.ascii` present in the execution directory?
- Is run intent benchmark mode (`1`) or production mode (`0`, VTK output)?
- For npFEM GPU, does hardware meet compute-capability expectations and CMake arch settings?
- For blood/deformable workflows, are `NumberOfNodes` and `NumberOfGPUsPerNode` consistent with launch?

### Canonical workflow
1. Select scenario and read the corresponding README (`examples/showCases/partialBounceBack/README.md`, `examples/gpuExamples/sandstone/README.md`, `examples/gpuExamples/multiComponentPorous/README.md`, `examples/showCases/bloodFlowDefoBodies/README.md`, `coupledSimulators/npFEM/npFEM_StandAlone_RhinoGH/README.md`).
2. Build with the recommended toolchain (`mkdir -p build && cd build && cmake .. && make -j`; for sandstone and multi-component GPU flows, set `export CXX=nvc++` before CMake).
3. Validate mandatory input assets (`clot_0.0003.txt`, `Berea.ascii`, and case XML/CP inputs when applicable).
4. Launch MPI/GPU command using the documented signature.
5. For production runs, verify output generation; for benchmark runs, disable heavy output paths.
6. If performance is poor, inspect doc-identified hotspots in `src/atomicBlock/atomicAcceleratedLattice3D.hh` and communication loops in `src/parallelism/parallelBlockCommunicator3D.cpp`.
7. For npFEM GPU constraints, inspect `coupledSimulators/npFEM/src_GPU/device_utilities.cu` and solver kernels.

### Minimal working example
```bash
cd examples/showCases/partialBounceBack
mkdir -p build && cd build
cmake .. && make -j
cd ..
mpirun -n 4 ./fibrinolysis
```

```bash
cd examples/gpuExamples/sandstone
cd build
export CXX=nvc++
cmake .. && make -j
cd ..
mpirun -np 4 ./sandstone 1 400 400 400
```

```bash
cd examples/gpuExamples/multiComponentPorous
cd build
export CXX=nvc++
cmake .. && make -j
cd ..
mpirun -np 4 ./multiComponentSandstone 1 400 400 400
```

### Pitfalls and fixes
- Sandstone fails immediately because `Berea.ascii` is missing in run directory (`examples/gpuExamples/sandstone/README.md`).
- Wrong compiler path for GPU build; set `CXX=nvc++` before build (`examples/gpuExamples/sandstone/README.md`).
- Confusing benchmark vs production flag (`1` vs `0`) causes unexpected output/performance behavior (`examples/gpuExamples/sandstone/README.md`).
- Startup appears slow (around minutes) due to CPU-side domain setup before GPU kernels; this is documented behavior (`examples/gpuExamples/sandstone/README.md`).
- npFEM GPU incompatibility on lower compute capability hardware; adjust arch/atomic strategy as documented (`examples/showCases/bloodFlowDefoBodies/README.md`).
- Large-RBC cases can require manual `maxThreadsPerBlock` tuning in `device_utilities.cu` (`examples/showCases/bloodFlowDefoBodies/README.md`).

### Convergence/validation checks
- Run short benchmark-mode sweeps over MPI sizes; check monotonic runtime scaling before long production runs.
- Confirm production mode emits expected VTK artifacts; benchmark mode should focus on timing.
- For partial bounce-back, verify clot map loading and physically stable pressure-driven flow (`examples/showCases/partialBounceBack/README.md`).
- For GPU runs, watch for communication overhead around packing/unpacking paths and confirm no kernel launch warnings.
- For sandstone and multi-component porous runs, verify `Berea.ascii` is found before timestep 0.
- Keep physical parameters fixed while changing parallel layout to isolate scaling effects.

## Scope
- Handle questions about MPI/OpenMP/GPU execution, scaling, and batch systems.
- Keep responses abstract and architectural for large codebases; avoid exhaustive per-function documentation unless requested.

## Primary documentation references
- `examples/showCases/partialBounceBack/README.md`
- `examples/gpuExamples/sandstone/README.md`
- `examples/gpuExamples/multiComponentPorous/README.md`
- `examples/showCases/bloodFlowDefoBodies/README.md`
- `coupledSimulators/npFEM/npFEM_StandAlone_RhinoGH/README.md`

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
- `examples/showCases/partialBounceBack/fibrinolysis.cpp`
- `examples/gpuExamples/sandstone/sandstone.cpp`
- `examples/gpuExamples/sandstone/runMultiGPU`
- `examples/gpuExamples/multiComponentPorous/multiComponentSandstone.cpp`
- `examples/gpuExamples/multiComponentPorous/runMultiGPU`
- `src/atomicBlock/atomicAcceleratedLattice3D.hh`
- `src/parallelism/parallelBlockCommunicator3D.cpp`
- `coupledSimulators/npFEM/npFEM_StandAlone_RhinoGH/api/API.h`
- `coupledSimulators/npFEM/npFEM_StandAlone_RhinoGH/api/API.cpp`
- `coupledSimulators/npFEM/npFEM_StandAlone_RhinoGH/CMakeLists.txt`
- `coupledSimulators/npFEM/src_GPU/Solver_GPU.h`
- `coupledSimulators/npFEM/src_GPU/Solver_GPU.cpp`
- `coupledSimulators/npFEM/src_GPU/device_utilities.cu`
- `coupledSimulators/npFEM/src_GPU/svd3_cuda.h`
- `coupledSimulators/npFEM/src_GPU/sum_cuda.h`
- `coupledSimulators/npFEM/src_GPU/sparse_matrix.h`
- Prefer targeted source search (for example: `rg -n "<symbol_or_keyword>" examples/gpuExamples examples/showCases coupledSimulators src`).
