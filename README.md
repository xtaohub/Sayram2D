# Sayram2D

Sayram2D is a C++20 code for two-dimensional Fokker–Planck diffusion on
structured, orthogonal, possibly nonuniform rectangular grids. It implements
the positivity-preserving finite-volume NTPFA discretization described in
[`docs/sayram2d_algorithm.md`](docs/sayram2d_algorithm.md).

The project is organized around a reusable core library and a set of examples.
The core library depends only on Eigen. The standalone Albert–Young
application is one example; HDF5 is used by that example, not by the library.

This version builds upon the original Sayram2D developed by Xin Tao and Peng
Peng. The redesign was led and reviewed by Xin Tao, with implementation carried
out by DeepSeek V4.1 Flash using DeepSeek Harness. Its purpose was to make the
code modular, easier to maintain, and easier to use in 3D or 4D applications.

For questions about using Sayram2D, bug reports, or inquiries regarding
potential collaborations, please contact Xin Tao at xtao@ustc.edu.cn.

If you use Sayram2D in your research, coauthorship with Xin Tao or any other
Sayram2D developer is not required. However, please cite the following
reference:

> Peng, P., Tao, X., Peng, Z., Jiang, Y., Gao, Z., Yang, D., et al. (2024).
> Modeling radiation belt dynamics using a positivity-preserving finite volume
> method on general meshes. *Journal of Geophysical Research: Space Physics*,
> 129, e2024JA032919. https://doi.org/10.1029/2024JA032919

## Documentation

1. [Core library](docs/core_library.md) — build and use instructions, API,
   numerical method, implementation notes, and repository layout.
2. [Examples](docs/examples.md) — library examples, the Albert–Young
   application, data, and scripts.

## License and Copyright

Copyright (c) 2026 Xin Tao.

This code is released under the [MIT License](LICENSE).
