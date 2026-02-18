# Evidence: palabos-parallel-hpc

## Primary docs
- `examples/showCases/partialBounceBack/README.md`
- `examples/gpuExamples/sandstone/README.md`
- `examples/gpuExamples/multiComponentPorous/README.md`
- `coupledSimulators/npFEM/npFEM_StandAlone_RhinoGH/README.md`

## Primary source entry points
- `skills/palabos-parallel-hpc/references/doc_map.md`
- `coupledSimulators/npFEM/npFEM_StandAlone_RhinoGH/api/API.h`
- `coupledSimulators/npFEM/npFEM_StandAlone_RhinoGH/api/API.cpp`
- `coupledSimulators/npFEM/npFEM_StandAlone_RhinoGH/CMakeLists.txt`
- `coupledSimulators/npFEM/src_GPU/Solver_GPU.h`
- `coupledSimulators/npFEM/src_GPU/Solver_GPU.cpp`
- `coupledSimulators/npFEM/src_GPU/svd3_cuda.h`
- `coupledSimulators/npFEM/src_GPU/sum_cuda.h`
- `coupledSimulators/npFEM/src_GPU/sparse_matrix.h`
- `coupledSimulators/npFEM/src_GPU/sparse_matrix.cpp`
- `coupledSimulators/npFEM/src_GPU/quasy_newton3.h`
- `coupledSimulators/npFEM/src_GPU/quasy_newton.h`
- `coupledSimulators/npFEM/src_GPU/projections_GPU_soa.h`
- `coupledSimulators/npFEM/src_GPU/projections_GPU_MATH.h`
- `coupledSimulators/npFEM/src_GPU/projections_GPU_deprecated.h`
- `coupledSimulators/npFEM/src_GPU/GPU_data.h`
- `coupledSimulators/npFEM/src_GPU/device_utilities.h`
- `coupledSimulators/npFEM/src_GPU/Constraint_Flattening.h`
- `coupledSimulators/npFEM/src_GPU/common.h`
- `coupledSimulators/npFEM/Solver.h`

## Extracted headings
- The Partially Saturated Method
- Compilation and execution
- Pore-level flow through a Berea sandstone
- Usage
- Performance-critical GPU code
- Ongoing improvements
- Instructions for npFEM (ShapeOp) with Rhino3D-Grasshopper
- Compilation
- Rhino-GH

## Executable command hints
- ./fibrinolysis
- mpirun -n NCores ./fibrinolysis

## Warnings and pitfalls
- (none extracted)
