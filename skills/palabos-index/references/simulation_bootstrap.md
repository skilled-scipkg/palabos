# Palabos Simulation Bootstrap (from repo root)

Use this file when you need a fast, reproducible start command and a minimal validation check.

## Preflight
- Verify tools: `cmake --version`, `mpirun --version`.
- Start from repository root (`palabos-palabos`).

## Quick-start runs

### 1) Installation smoke test (2D cavity)
```bash
cd examples/showCases/cavity2d
cd build && cmake .. && make -j
cd ..
./cavity2d
```
Validation checkpoints:
- Console prints step/energy/density updates.
- `examples/showCases/cavity2d/tmp/` contains fresh GIF/VTK outputs.

### 2) Collision-model stress test (DSL2D)
```bash
cd examples/showCases/dsl2d
cd build && cmake .. && make -j
cd ..
mpirun -np 2 ./dsl2d config.xml 10000 128 0.1
```
Validation checkpoints:
- No immediate blow-up for baseline parameters.
- `tmp/` outputs are updated and log file is written.

### 3) Porous clot benchmark (partial bounce-back)
```bash
cd examples/showCases/partialBounceBack
mkdir -p build && cd build
cmake .. && make -j
cd ..
mpirun -n 4 ./fibrinolysis
```
Required input:
- `examples/showCases/partialBounceBack/clot_0.0003.txt`
Validation checkpoints:
- Run starts without missing-input errors.
- Output appears under `tmp/`.

### 4) GPU porous benchmark (sandstone)
```bash
cd examples/gpuExamples/sandstone
cd build
export CXX=nvc++
cmake .. && make -j
cd ..
mpirun -np 4 ./sandstone 1 400 400 400
```
Required input:
- `Berea.ascii` in `examples/gpuExamples/sandstone/`
Validation checkpoints:
- Startup passes geometry load.
- Benchmark mode (`1`) runs without VTK overhead.

### 5) Full workflow run (blood cells + aneurysm)
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
Validation checkpoints:
- Blood-flow cell packing creates CP outputs before production run.
- Aneurysm run accepts `param.xml` and progresses through refinement levels.

## If behavior is unclear
- Use topic skills and their `references/source_map.md` files for function-level drill-down before broad code search.
