---
name: palabos-developer-guide
description: This skill should be used when users ask about developer guide in palabos; it prioritizes documentation references and then source inspection only for unresolved details.
---

# palabos: Developer Guide

## High-Signal Playbook

### Route conditions
- Stay in this skill for contribution legality, DCO certification, and licensing expectations.
- Route end-user simulation setup questions to `palabos-getting-started`, `palabos-simulation-workflows`, or `palabos-examples-and-tutorials`.
- Route MPI/GPU deployment/scaling to `palabos-parallel-hpc`.
- For architecture internals not covered by docs, use this skill as the escalation gateway to concrete source entry points.

### Triage questions
- Is the request legal/compliance (DCO/license) or technical implementation behavior?
- For contributions, which DCO path applies: original work, compatible prior work, or third-party certified work?
- Is AGPLv3-or-later licensing of contributions understood and acceptable?
- Is a persistent public contribution record (including sign-off identity) acceptable?
- Is the user requesting repository-extension guidance that needs source-level inspection?

### Canonical workflow
1. Start from `dco/README.md` and classify whether the question is legal/compliance or technical.
2. For compliance, map the contribution to DCO clauses (a)-(d) and confirm AGPLv3+ licensing expectations.
3. Keep the response auditable by citing `dco/README.md` paths and exact relevant clauses.
4. If the request is technical and docs are insufficient, escalate to implementation hubs (`src/palabos2D.h`, `src/palabos3D.h`, `src/core/dynamics.h`).
5. For subsystem-specific behavior, route to the corresponding topic skill before deeper source inspection.

### Minimal working example
```text
Contribution checklist (from dco/README.md):
- Origin is valid under DCO clause (a), (b), or (c)
- Contribution is accepted under AGPLv3-or-later
- Contributor accepts permanent public record of sign-off/metadata
```

### Pitfalls and fixes
- Contributors cannot certify origin rights under DCO clauses; resolve provenance before code submission (`dco/README.md`).
- License mismatch assumptions (non-AGPL-compatible code) can block contribution acceptance (`dco/README.md`).
- Public-record permanence is overlooked; clarify that contribution metadata is maintained indefinitely (`dco/README.md`).
- Legal/compliance questions get mixed with runtime debugging; split handling, then route technical work to the right topic skill.
- Developers jump into random internals without source entry anchors; start from core headers and nearest example drivers.

### Convergence/validation checks
- Every contribution request maps cleanly to one DCO certification path.
- License expectations are explicit (AGPLv3 or newer) before implementation advice proceeds.
- Compliance response cites `dco/README.md` instead of memory-based summaries.
- Technical escalation identifies at least one concrete source file and one owning topic skill.
- Final guidance is reproducible: legal checks first, then source drill-down only if needed.

## Scope
- Handle questions about developer architecture, extension points, and contribution workflow.
- Keep responses abstract and architectural for large codebases; avoid exhaustive per-function documentation unless requested.

## Primary documentation references
- `dco/README.md`

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
- `src/palabos2D.h`
- `src/palabos2D.hh`
- `src/palabos3D.h`
- `src/palabos3D.hh`
- `src/core/dynamics.h`
- `src/core/dynamics.hh`
- `src/core/runTimeDiagnostics.h`
- `src/core/runTimeDiagnostics.cpp`
- `src/boundaryCondition/generalizedBoundaryDynamicsSolvers.h`
- `src/boundaryCondition/generalizedBoundaryDynamicsSolvers.hh`
- `coupledSimulators/npFEM/Solver.h`
- `coupledSimulators/npFEM/Solver.cpp`
- Prefer targeted source search (for example: `rg -n "<symbol_or_keyword>" coupledSimulators src utility`).
