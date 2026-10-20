# Sayram2D Redesign: Objectives and Design Principles

## Purpose

Redesign Sayram2D as a modular, reusable C++ solver for two-dimensional Fokker–Planck (FP) diffusion equations. The numerical core should serve both standalone simulations and as a component within higher-dimensional models, while exposing convenient Python and C++ command-line interfaces.

The initial redesign should **preserve the validated numerical method and results**. Architectural refactoring and changes to numerical algorithms should be separate, independently testable steps.

## Primary objectives

### 1. Modularize the code

Separate responsibilities behind clear, minimal interfaces. In particular, distinguish:

- **Mesh and coordinates:** structured rectangular 2D grid geometry and coordinate information, allowing nonuniform spacing along either coordinate.
- **FP equation and geometry factors:** coefficients and terms associated with the chosen coordinate system, including the Jacobian factor \(G\).
- **Boundary conditions:** independently configurable treatments for each relevant boundary.
- **Diffusion coefficients:** a coherent object containing the tensor components, rather than several easily confused positional arguments.
- **Numerical solver:** discretization, matrix assembly, linear solution, and reusable workspace.
- **Simulation driver and I/O:** time evolution, initialization, configuration, and output, outside the reusable solver core.

Avoid excessive abstraction: introduce new classes or templates only when they clarify ownership, enable a required variation, or remove meaningful coupling.

### 2. Make boundary conditions easy to configure, replace, and extend

Boundary-condition logic should not be entangled with the main solver algorithm. It should be possible to select or change the condition on each boundary without rewriting matrix assembly or unrelated code.

Support the boundary conditions already implemented by Sayram2D, and design interfaces that can accommodate additional types later. Sayram2D uses **structured, rectangular 2D grids**: grid lines follow the two coordinate directions and cells may have nonuniform widths in either direction. Unstructured meshes and arbitrary cell shapes are **out of scope**. Boundary-condition handling must work on these structured, potentially nonuniform grids and preserve the existing numerical discretization. Define explicitly how boundary values or fluxes enter the discrete equations.

### 3. Support alternative coordinates and their Jacobian \(G\)

The FP solver must not hard-code the Jacobian \(G\) for its current coordinate choice. Users may adapt Sayram2D to another pair of coordinates by providing the appropriate \(G\) in the FP formulation. **Coordinate switching during a simulation is not required.**

Keep the Jacobian \(G\) configurable and separate from the generic solution machinery wherever practical. Clearly specify its mathematical convention and placement in the discrete FP equation. The solver should use the supplied \(G\) without embedding assumptions about a particular physical coordinate system.

**Scope and responsibility:** Sayram2D is not responsible for transforming diffusion coefficients or boundary conditions between coordinate systems. Users choosing different coordinates must supply diffusion coefficients and boundary conditions consistent with those coordinates. Within the solver, the intended coordinate-dependent change is the supplied Jacobian \(G\); no general coordinate-transformation machinery is required.

### 4. Make the 2D solver callable from 3D and 4D codes

Treat Sayram2D as a library component, not as a program that owns the full simulation lifecycle. A higher-dimensional driver should be able to solve independent 2D FP problems at fixed values of the remaining coordinates (for example, at each L shell or at each \((L,\mathrm{MLT})\) location).

- The **caller owns the evolving distribution** and controls time stepping/operator splitting.
- The solver accepts input and output 2D fields and a time step; it must not require a standalone simulation loop.
- Diffusion coefficients can be updated independently of a solve, allowing their update frequency to differ from the time-step frequency.
- Prefer compatibility with existing **xtensor** arrays and lightweight 2D views of higher-dimensional arrays, avoiding unnecessary slice copies.
- Define ownership, lifetimes, and thread-safety of per-solve workspaces so independent 2D solves can run in parallel.

An illustrative target API (not a fixed implementation specification) is:

```cpp
Solver solver(mesh, equation);
solver.set_diffusion(D);
solver.solve(input, output, dt);
```

### 5. Provide both Python and C++ command-line interfaces

Build around a **single C++ numerical core**, exposed through two entry points:

- **Python interface:** configure a problem, provide arrays and coefficients, invoke the solver, and retrieve results without duplicating the numerical implementation. Aim for natural interoperability with NumPy arrays and document copy/view and lifetime semantics.
- **C++ command-line executable:** run standalone problems from configuration/input files, with reproducible settings and clearly defined output. Keep argument parsing, file handling, and simulation orchestration outside the core solver library.

The Python and CLI interfaces should share the same numerical behavior and be tested against the same reference cases. The choice of binding technology and configuration format can be made after reviewing the current code and build system.

### 6. Use the numerical algorithm document as the reference

Use **`sayram2d_algorithm.pdf`** (*Sayram 2D: Numerical Algorithm*) as the reference for the intended numerical formulation. Preserve the existing validated behavior during architectural refactoring; if the current implementation differs from the document, report the difference and propose a separately testable correction rather than silently changing the algorithm.

- Preserve the positivity-preserving finite-volume formulation based on the nonlinear two-point flux approximation (NTPFA), including local flux conservation and the nonnegative-solution property.
- Support **orthogonal rectangular cells with potentially nonuniform spacing** in both coordinates; arbitrary quadrilateral and unstructured meshes are not required.
- Use the documented bilinear interpolation of interior vertex values, and handle boundary vertices consistently with the boundary conditions.
- Assemble the discrete system consistently with the documented $M$-matrix argument, including the Jacobian $G$.
- Treat **losses implicitly** ($\tau_K>0$) and **sources explicitly** ($\tau_K<0$), as described in the document.
- Keep the solver's numerical core separate from simulation drivers and I/O, while retaining a simple `solve(input, output, dt)`-style interface suitable for 2D views into higher-dimensional arrays.
- Establish baseline regression tests **before modifying the implementation**. Test positivity, local conservation, spatial convergence, boundary conditions, nonuniform grids, and loss/source handling. Distinguish tests of existing behavior from tests of the intended documented formulation.

## Constraints and acceptance criteria

1. **Numerical fidelity:** Existing test cases and boundary-condition behavior remain reproducible within documented numerical tolerances.
2. **Clear separation:** Mesh, coordinate/FP formulation, diffusion data, boundary conditions, solver, and drivers have distinct responsibilities.
3. **Coordinate extensibility:** A second coordinate formulation can be introduced without duplicating the entire solver or modifying unrelated components.
4. **Boundary-condition extensibility:** A boundary condition can be replaced or added through a localized implementation and targeted tests.
5. **Embedding:** A 3D/4D driver can call the 2D solver on array slices, with externally managed simulation state.
6. **Dual interfaces:** Python and C++ CLI examples solve an equivalent problem and agree numerically.
7. **Pragmatic scope:** Retain xtensor and Eigen initially; avoid changing the discretization or linear algebra backend solely as part of modularization.
8. **Algorithm reference:** Document discrepancies between the existing implementation and `sayram2d_algorithm.pdf`; address them only in separate, validated changes.
9. **Numerical validation:** Regression and targeted tests cover positivity, conservation, convergence, nonuniform grids, boundary conditions, and source/loss terms.

## Recommended workflow for the local AI agent

**Phase 1 — Audit and design (no code changes).** Read `sayram2d_algorithm.pdf` and inspect the current codebase, dependency graph, equation/Jacobian conventions, boundary-condition implementations, public interfaces, and existing tests. Identify concrete coupling problems and propose a minimal architecture, including ownership and data-flow diagrams. Flag ambiguities in the FP formulation rather than guessing.

**Phase 2 — Incremental refactoring.** Establish baseline numerical regression tests first, then extract reusable components and introduce the library-facing API while keeping the current numerical scheme. Make small, reviewable commits, each with passing tests.

**Phase 3 — Extensibility tests.** Demonstrate swapping boundary conditions, introducing an alternative coordinate/Jacobian formulation, and solving 2D slices from representative 3D/4D xtensor arrays.

**Phase 4 — User interfaces.** Add the Python bindings and C++ CLI as thin adapters to the same C++ core; document installation, inputs, outputs, and example usage.

**Phase 5 — Validation and documentation.** Compare against baseline results, test nonuniform meshes and boundary cases, and document the stable C++/Python APIs and their data ownership contracts.

## Initial instruction to the agent

> Review the existing Sayram2D implementation against the objectives above and `sayram2d_algorithm.pdf`, which specifies the intended numerical formulation. **Do not rewrite code yet.** First provide a concrete architecture proposal, identify which existing classes/functions can be retained, explain how boundary conditions and coordinate-dependent Jacobian \(G\) will be handled, and propose an incremental implementation and regression-testing plan. Keep the design as simple as the requirements allow.
