"""Time-dependent D and Dirichlet boundary example for the Python binding."""

import numpy as np

import sayram2d


def main() -> None:
    n1 = n2 = 40
    dt = 0.01
    steps = 100
    tau = 10.0

    q1_edges = np.linspace(0.0, 1.0, n1 + 1)
    q2_edges = np.linspace(0.0, 1.0, n2 + 1)
    mesh = sayram2d.Mesh(q1_edges, q2_edges)

    G = np.asfortranarray(np.ones((n1, n2)))
    D11 = np.asfortranarray(np.ones((n1, n2)))
    D12 = np.asfortranarray(np.full((n1, n2), 0.8))
    D22 = np.asfortranarray(np.ones((n1, n2)))

    bc = sayram2d.BoundaryConditions()
    bc.set_zero_flux(sayram2d.Side.Q1Max)
    bc.set_zero_flux(sayram2d.Side.Q2Min)
    bc.set_zero_flux(sayram2d.Side.Q2Max)

    boundary_values = np.empty(n2 + 1)
    equation = sayram2d.Equation(G, D11, D12, D22, bc)
    solver = sayram2d.Solver(mesh)

    f = np.asfortranarray(np.ones((n1, n2)))
    f_next = np.empty((n1, n2), order="F")

    for step in range(steps + 1):
        time = step * dt
        decay = np.exp(-time / tau)

        # Mutate the arrays in place: Equation holds views of this memory.
        D11[:, :] = decay
        D12[:, :] = 0.8 * decay
        D22[:, :] = decay

        boundary_values[:] = 0.5 * decay
        bc.set_dirichlet(sayram2d.Side.Q1Min, boundary_values)

        if step == steps:
            break

        solver.solve(equation, f, f_next, dt)
        f, f_next = f_next, f

        if step in (0, steps // 2):
            print(
                f"t={time:.2f} D11={D11[0, 0]:.6f} "
                f"bc={boundary_values[0]:.6f} min={f.min():.6f} max={f.max():.6f}"
            )


if __name__ == "__main__":
    main()
