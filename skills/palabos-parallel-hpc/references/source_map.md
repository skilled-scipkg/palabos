# palabos source map: Parallel and HPC

Use this map after the HPC docs in `references/doc_map.md`.

## Fast source navigation
- `rg -n "<symbol_or_keyword>" examples/gpuExamples examples/showCases coupledSimulators src`
- `rg -n "int main|collideAndStream|writeVTK|CUDA_VISIBLE_DEVICES|Solver_GPU" examples/gpuExamples coupledSimulators/npFEM`

## Function-level source entry points

### MPI / GPU example drivers
- `examples/showCases/partialBounceBack/fibrinolysis.cpp` | `main`, partial-bounce-back setup, clot-map loading.
- `examples/gpuExamples/sandstone/sandstone.cpp` | `main`, accelerated `collideAndStream`, `writeVTK`.
- `examples/gpuExamples/sandstone/runMultiGPU` | MPI local-rank to `CUDA_VISIBLE_DEVICES` mapping.
- `examples/gpuExamples/multiComponentPorous/multiComponentSandstone.cpp` | `main`, dual-lattice GPU path, `writeVTK`.
- `examples/gpuExamples/multiComponentPorous/runMultiGPU` | multi-GPU rank pinning helper.

### Performance-critical engine code
- `src/atomicBlock/atomicAcceleratedLattice3D.hh` | collision and pack/unpack kernels referenced by sandstone README.
- `src/parallelism/parallelBlockCommunicator3D.cpp` | communication loop surrounding repeated pack/unpack operations.

### npFEM CPU/GPU coupling internals
- `coupledSimulators/npFEM/Solver.h` | CPU solver API (`initialize`, `solve`, `collisionDetection`).
- `coupledSimulators/npFEM/Solver.cpp` | CPU solver execution path and convergence code.
- `coupledSimulators/npFEM/src_GPU/Solver_GPU.h` | GPU solver interface and data movement API.
- `coupledSimulators/npFEM/src_GPU/Solver_GPU.cpp` | GPU solver step loop and host/device orchestration.
- `coupledSimulators/npFEM/src_GPU/device_utilities.cu` | CUDA memory setup, kernels, and launch-bound-sensitive utilities.
- `coupledSimulators/npFEM/npFEM_StandAlone_RhinoGH/api/API.h` | stand-alone API used for Rhino/GH integration context.
- `coupledSimulators/npFEM/npFEM_StandAlone_RhinoGH/api/API.cpp` | stand-alone API implementation.
