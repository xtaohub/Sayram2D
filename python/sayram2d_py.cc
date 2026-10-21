#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "boundary_conditions.h"
#include "equation.h"
#include "mesh.h"
#include "solver.h"

namespace py = pybind11;
using namespace sayram2d;

namespace {

// Require Fortran/column-major contiguous 2D double arrays. This matches the
// in-memory Array2D/ArrayView2D convention and keeps the binding zero-copy.
py::array_t<double> as_fortran_2d(py::array_t<double> array,
                                  const char* name) {
  const auto info = array.request();
  if (info.ndim != 2) {
    throw std::invalid_argument(std::string(name) + " must be 2D");
  }
  const auto n1 = static_cast<std::size_t>(info.shape[0]);
  const auto n2 = static_cast<std::size_t>(info.shape[1]);
  if (info.strides[0] != static_cast<py::ssize_t>(sizeof(double)) ||
      info.strides[1] !=
          static_cast<py::ssize_t>(sizeof(double) * n1)) {
    throw std::invalid_argument(
        std::string(name) +
        " must be Fortran-contiguous; use np.asfortranarray()");
  }
  if (array.size() != static_cast<py::ssize_t>(n1 * n2)) {
    throw std::invalid_argument(std::string(name) + " has an invalid size");
  }
  return array;
}

ArrayView2D<const double> const_view(const py::array_t<double>& array) {
  const auto info = array.request();
  return ArrayView2D<const double>(
      static_cast<const double*>(info.ptr),
      static_cast<std::size_t>(info.shape[0]),
      static_cast<std::size_t>(info.shape[1]));
}

ArrayView2D<double> mutable_view(py::array_t<double> array) {
  const auto info = array.request();
  if (info.readonly) {
    throw std::invalid_argument("output array must be writable");
  }
  return ArrayView2D<double>(static_cast<double*>(info.ptr),
                             static_cast<std::size_t>(info.shape[0]),
                             static_cast<std::size_t>(info.shape[1]));
}

class PyBoundaryConditions {
 public:
  BoundaryConditions& get() noexcept { return bc_; }
  const BoundaryConditions& get() const noexcept { return bc_; }

  void set_dirichlet(Side side, double value) {
    bc_.set_dirichlet(side, value);
  }

  void set_dirichlet(Side side,
                     py::array_t<double, py::array::c_style |
                                             py::array::forcecast>
                         values) {
    const auto info = values.request();
    if (info.ndim != 1) {
      throw std::invalid_argument("Dirichlet values must be one-dimensional");
    }
    const std::size_t index = static_cast<std::size_t>(side);
    auto& storage = values_[index];
    const auto* data = static_cast<const double*>(info.ptr);
    storage.assign(data, data + info.shape[0]);
    bc_.set_dirichlet(side, storage);
  }

  void set_zero_flux(Side side) { bc_.set_zero_flux(side); }
  void validate(const Mesh& mesh) const { bc_.validate(mesh); }
  BCType type(Side side) const { return bc_.type(side); }
  bool is_dirichlet(Side side) const { return bc_.is_dirichlet(side); }

 private:
  BoundaryConditions bc_;
  std::array<std::vector<double>, 4> values_{};
};

class PyEquation {
 public:
  PyEquation(py::array_t<double> G, py::array_t<double> D11,
             py::array_t<double> D12, py::array_t<double> D22,
             py::object bc)
      : G_(as_fortran_2d(G, "G")),
        D11_(as_fortran_2d(D11, "D11")),
        D12_(as_fortran_2d(D12, "D12")),
        D22_(as_fortran_2d(D22, "D22")),
        bc_object_(std::move(bc)),
        bc_(bc_object_.cast<PyBoundaryConditions&>()),
        equation_(const_view(G_), const_view(D11_), const_view(D12_),
                  const_view(D22_), bc_.get()) {}

  Equation& get() noexcept { return equation_; }
  const Equation& get() const noexcept { return equation_; }

  void set_loss(double value) { equation_.set_loss(value); }

  void set_loss(py::array_t<double> values) {
    loss_ = as_fortran_2d(values, "loss");
    equation_.set_loss(const_view(loss_));
  }

 private:
  py::array_t<double> G_;
  py::array_t<double> D11_;
  py::array_t<double> D12_;
  py::array_t<double> D22_;
  py::array_t<double> loss_;
  py::object bc_object_;
  PyBoundaryConditions& bc_;
  Equation equation_;
};

}  // namespace

PYBIND11_MODULE(sayram2d, m) {
  m.doc() = "Python bindings for Sayram2D";
  m.attr("__version__") = "0.1.0";

  py::enum_<Side>(m, "Side")
      .value("Q1Min", Side::Q1Min)
      .value("Q1Max", Side::Q1Max)
      .value("Q2Min", Side::Q2Min)
      .value("Q2Max", Side::Q2Max);

  py::enum_<BCType>(m, "BCType")
      .value("Dirichlet", BCType::Dirichlet)
      .value("ZeroFlux", BCType::ZeroFlux);

  py::class_<Mesh>(m, "Mesh")
      .def(py::init([](std::vector<double> q1_edges,
                       std::vector<double> q2_edges) {
             return Mesh(std::move(q1_edges), std::move(q2_edges));
           }),
           py::arg("q1_edges"), py::arg("q2_edges"))
      .def_property_readonly("n1", &Mesh::n1)
      .def_property_readonly("n2", &Mesh::n2)
      .def("q1_center", &Mesh::q1_center)
      .def("q2_center", &Mesh::q2_center)
      .def("q1_edge", &Mesh::q1_edge)
      .def("q2_edge", &Mesh::q2_edge);

  py::class_<PyBoundaryConditions>(m, "BoundaryConditions")
      .def(py::init<>())
      .def("set_dirichlet",
           py::overload_cast<Side, double>(
               &PyBoundaryConditions::set_dirichlet),
           py::arg("side"), py::arg("value"))
      .def("set_dirichlet",
           py::overload_cast<
               Side, py::array_t<double, py::array::c_style |
                                              py::array::forcecast>>(
               &PyBoundaryConditions::set_dirichlet),
           py::arg("side"), py::arg("values"))
      .def("set_zero_flux", &PyBoundaryConditions::set_zero_flux)
      .def("validate", &PyBoundaryConditions::validate)
      .def("type", &PyBoundaryConditions::type)
      .def("is_dirichlet", &PyBoundaryConditions::is_dirichlet);

  py::class_<PyEquation>(m, "Equation")
      .def(py::init<py::array_t<double>, py::array_t<double>,
                    py::array_t<double>, py::array_t<double>, py::object>(),
           py::arg("G"), py::arg("D11"), py::arg("D12"), py::arg("D22"),
           py::arg("bc"))
      .def("set_loss",
           py::overload_cast<double>(&PyEquation::set_loss),
           py::arg("value"))
      .def("set_loss",
           py::overload_cast<py::array_t<double>>(&PyEquation::set_loss),
           py::arg("values"));

  py::class_<Solver>(m, "Solver")
      .def(py::init<const Mesh&>(), py::arg("mesh"),
           py::keep_alive<1, 2>())
      .def(
          "solve",
          [](Solver& solver, PyEquation& equation, py::array_t<double> f_in,
             py::array_t<double> f_out, double dt) {
            auto input = as_fortran_2d(f_in, "f_in");
            auto output = as_fortran_2d(f_out, "f_out");
            const auto input_view = const_view(input);
            const auto output_view = mutable_view(output);
            py::gil_scoped_release release;
            solver.solve(equation.get(), input_view, output_view, dt);
          },
          py::arg("equation"), py::arg("f_in"), py::arg("f_out"),
          py::arg("dt"));
}
