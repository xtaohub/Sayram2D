# Sayram-2D

**Sayram-2D** is a C++ program within the **Sayram numerical framework**, designed to solve the **two-dimensional radiation belt diffusion equation** using a **positivity-preserving finite-volume method (PPFV)**. The solver supports a full diffusion tensor, is formulated in conservative form, and guarantees non-negativity of the solution without imposing CFL-type time-step restrictions.

While Sayram-2D is specifically developed for modeling **radiation belt dynamics**, the underlying numerical method is general. With minor modifications, the code can also be applied to **other two-dimensional diffusion equations** written in conservative form.

Prepare to extend Sayram-2D further so it can be easily used in 3D and 4D codes.

> **Note**  
> This code was previously released under the repository name **`fvm2d`**.  
> It has been renamed **Sayram-2D** to better reflect its role within the Sayram framework.  
> The numerical method and scientific results remain unchanged.

---

## Overview

In its current public implementation, Sayram-2D solves the two-dimensional Fokker–Planck equation governing radiation belt electron dynamics. The equation is formulated in **($\alpha_0$, $\log E$)** coordinates, where $\alpha_0$ is the equatorial pitch angle and $E$ is the particle kinetic energy.

The code advances the phase-space density on a structured finite-volume mesh using a full diffusion tensor, including cross-diffusion terms. This allows stable and accurate treatment of pitch-angle diffusion, energy diffusion, and their coupling arising from local wave–particle interactions.

Although the present implementation targets radiation belt applications, the numerical infrastructure is sufficiently general that, with slight modifications, it can be adapted to solve **other two-dimensional diffusion problems** written in conservative form.

Key features:

- Positivity-preserving finite-volume discretization  
- Full diffusion tensor, including cross-diffusion terms  
- No CFL-type time-step restriction  
- Extensible to general Fokker–Planck equations in conservative form  

---

## Governing Equation

Sayram-2D solves the two-dimensional diffusion equation

$$
\frac{\partial f}{\partial t} = \sum_{i,j} \frac{1}{G}\frac{\partial}{\partial Q_i}(GD_{Q_iQ_j}\frac{\partial f}{\partial Q_j}),
$$

where  
- $f$ is the phase-space density,  
- $Q_i$ are generalized coordinates,  
- $D_{Q_i Q_j}$ is the diffusion tensor, and  
- $G$ is the Jacobian factor.

The numerical scheme is based on the **positivity-preserving finite-volume (PPFV) method** developed by Gao and Wu, which guarantees non-negativity of the solution without imposing CFL constraints on the time step.

---

## Code Structure

- **`source/`**  
  Core numerical implementation: mesh construction, equation assembly, solver, and time integration.

- **`D/`**  
  Input diffusion coefficients.

- **`output/`**  
  Simulation outputs (default format: HDF5).

- **`plot/`**  
  Python scripts for visualization and comparison with reference results.  
  In particular, `cmp_ay.py` compares the default output with the results of Albert and Young (2005).

- **`p.ini`**  
  Default configuration file reproducing the test case of Albert and Young (2005, *GRL*).

---

## Dependencies

Sayram-2D depends on the following open-source C++ libraries:

- [Eigen](https://eigen.tuxfamily.org) (header-only)  
- [xtensor](https://github.com/xtensor-stack/xtensor) (header-only)  
- [xtl](https://github.com/xtensor-stack/xtl) (header-only)  
- [HDF5](https://github.com/HDFGroup/hdf5)

The first three libraries can be placed directly in the `source/` directory if desired.  
HDF5 must be installed and linked properly, as it is used for input and output by default.

---

## Compilation

### Makefile Configuration

Edit the following section in the `Makefile` according to your system setup:

```makefile
# -----------------
CC = clang++
LOCAL_INCLUDE = /Users/xtao/local/include
HDF5_INCLUDE = /opt/local/include
HDF5_LIB = /opt/local/lib
OPENMP =
OPT = -O2
# -----------------

```

Adjust compiler options and library paths as needed.

### Compile

Generate the executable by running:

```bash
make
```

If you encounter compilation issues, clean previous build files and recompile:

```bash
make clean
make
```

### Running the Code

Run the code with the default input parameter file p.ini:

```bash
./sayram-2d
```
The default input parameter file is "p.ini", which deals with the case in the paper by Albert and Young, GRL, 2005. 

To use a custom configuration file (e.g., new.ini), run:

```bash
./sayram-2d new.ini
```

## Boundary Conditions and Extensions

- **Default configuration**  
  The default setup reproduces the test case of Albert and Young (2005, *GRL*), with the boundary condition
  $$
  f(\alpha_0 = \alpha_{0,\text{LC}}) = 0 .
  $$

- **Alternative boundary condition**  
  An alternative condition,
  $$
  \left.\frac{\partial f}{\partial \alpha_0}\right|_{\alpha_0 = 0} = 0 ,
  $$
  can be implemented by modifying the corresponding source files. An example implementation is provided in the case **Albert_Young_LC**.

  To simulate this case, make the following changes:

  1. In `main.cc`, replace
     ```cpp
     #include "Albert_Young.h"
     ...

     Albert_Young eq(paras, m);
     ```
     with
     ```cpp
     #include "Albert_Young_LC.h"
     ...

     Albert_Young_LC eq(paras, m);
     ```

  2. Run the code using
     ```bash
     ./sayram-2d p_AlbertYoungLC.ini
     ```

- **Diffusion coefficients**  
  The input diffusion coefficients follow those of Albert and Young (2005, *GRL*). All three coefficients have units of \([p^2]/[t]\).

- **Additional loss or source term**  
  An additional loss (or source) term of the form \(f/\tau\) is incorporated directly into the coefficient matrix (with thanks to Mr. Bernhard Haas, GFZ).

- **Extensions**  
  Time-dependent diffusion coefficients or boundary conditions can be implemented by extending the **Equation** class and adding new cases under the **Cases** directory.

## Contributing

If you have any suggestions, please contact Peng Peng at pp140594 "AT" mail.ustc.edu.cn or Xin Tao at xtao "AT" ustc.edu.cn. 

## License and Citation

This code is released under the MIT License.

If you use Sayram-2D in your research, please cite:

Peng, P., Tao, X., Peng, Z., Jiang, Y., Gao, Z., Yang, D., et al. (2024). Modeling radiation belt dynamics using a positivity-preserving finite volume method on general meshes. Journal of Geophysical Research: Space Physics, 129, e2024JA032919. https://doi.org/10.1029/2024JA032919
