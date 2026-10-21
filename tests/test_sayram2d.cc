#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <Eigen/Dense>

#include "albert_young_case.h"
#include "reference_case.h"
#include "array_view2d.h"
#include "boundary_conditions.h"
#include "container_adapter.h"
#include "equation.h"
#include "mesh.h"
#include "physics.h"
#include "solver.h"

namespace {

using namespace sayram2d;
using reference::CaseConfig;
using reference::CaseData;

void check(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

bool near(double a, double b, double tolerance = 1e-14) {
  const double scale = std::max({1.0, std::abs(a), std::abs(b)});
  return std::abs(a - b) <= tolerance * scale;
}

template <typename Function>
void expect_throw(Function&& function, const std::string& message) {
  try {
    function();
  } catch (const std::exception&) {
    return;
  }
  throw std::runtime_error(message);
}

double max_abs_diff(const std::vector<double>& a,
                    const std::vector<double>& b) {
  check(a.size() == b.size(), "max_abs_diff: size mismatch");
  double result = 0.0;
  for (std::size_t k = 0; k < a.size(); ++k) {
    result = std::max(result, std::abs(a[k] - b[k]));
  }
  return result;
}

double max_abs_diff(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b) {
  check(a.rows() == b.rows() && a.cols() == b.cols(),
        "max_abs_diff: Eigen size mismatch");
  return (a - b).cwiseAbs().maxCoeff();
}

void test_mesh_geometry() {
  Mesh mesh({0.0, 1.0, 3.0, 6.0}, {0.0, 2.0, 5.0});
  check(mesh.n1() == 3 && mesh.n2() == 2, "mesh sizes");
  check(near(mesh.q1_center(0), 0.5) && near(mesh.q1_center(1), 2.0) &&
            near(mesh.q1_center(2), 4.5),
        "q1 centers");
  check(near(mesh.q2_center(0), 1.0) && near(mesh.q2_center(1), 3.5),
        "q2 centers");
  check(near(mesh.dq1(0), 1.0) && near(mesh.dq1(1), 2.0) &&
            near(mesh.dq1(2), 3.0),
        "q1 widths");
  check(near(mesh.dq2(0), 2.0) && near(mesh.dq2(1), 3.0), "q2 widths");
  check(near(mesh.cell_measure({1, 1}), 6.0), "cell measure");
  check(mesh.cell_index({2, 1}) == 5, "cell index");
  static_assert(flat_index(2, 1, 3) == 5);
  check(flat_index(1, 2, 3) == 7, "flat index utility");
  check(mesh.is_boundary_vertex({0, 1}), "boundary vertex");
  check(!mesh.is_boundary_vertex({1, 1}), "interior vertex");

  const Face west = mesh.face({1, 1}, Side::Q1Min);
  check(west.owner.i == 1 && west.owner.j == 1, "face owner");
  check(near(west.a.q1, 1.0) && near(west.a.q2, 5.0), "west A");
  check(near(west.b.q1, 1.0) && near(west.b.q2, 2.0), "west B");
  check(near(west.length, 3.0), "west length");
  check(near(west.normal.q1, -1.0) && near(west.normal.q2, 0.0), "west normal");
  check(!west.boundary && west.neighbor.i == 0 && west.neighbor.j == 1,
        "west neighbor");

  const Face east = mesh.face({0, 1}, Side::Q1Max);
  check(near(east.a.q1, 1.0) && near(east.a.q2, 2.0), "east A");
  check(near(east.b.q1, 1.0) && near(east.b.q2, 5.0), "east B");
  check(near(east.normal.q1, 1.0), "east normal");
  check(near(east.length, west.length), "shared face length");
  check(east.neighbor.i == 1 && east.neighbor.j == 1, "east neighbor");

  const Face south = mesh.face({1, 1}, Side::Q2Min);
  check(near(south.length, 2.0), "south length");
  check(near(south.normal.q2, -1.0), "south normal");
  const Face north = mesh.face({1, 1}, Side::Q2Max);
  check(near(north.length, 2.0), "north length");
  check(near(north.normal.q2, 1.0), "north normal");

  Ind neighbor{};
  check(!mesh.neighbor({0, 0}, Side::Q1Min, &neighbor), "west boundary");
  check(!mesh.neighbor({0, 0}, Side::Q2Min, &neighbor), "south boundary");
  check(mesh.is_boundary_face({0, 0}, Side::Q1Min), "boundary face");
  check(!mesh.is_boundary_face({1, 1}, Side::Q1Min), "interior face");

  expect_throw([] { Mesh m({0.0, 1.0, 1.0}, {0.0, 1.0}); },
               "non-increasing q1 must throw");
  expect_throw([] { Mesh m({0.0, 1.0}, {0.0, 2.0, 1.0}); },
               "non-increasing q2 must throw");
  expect_throw([] { Mesh m({0.0}, {0.0, 1.0}); },
               "too few q1 edges must throw");
  expect_throw([] {
    Mesh m({0.0, std::numeric_limits<double>::quiet_NaN()}, {0.0, 1.0});
  }, "non-finite edge must throw");
}

struct CustomMatrix2D {
  std::vector<double> values;
  std::size_t rows = 0;
  std::size_t cols = 0;
  std::ptrdiff_t s1 = 0;
  std::ptrdiff_t s2 = 0;

  double* data() { return values.data(); }
  const double* data() const { return values.data(); }
  std::size_t dimension() const { return 2; }
  std::size_t shape(std::size_t axis) const {
    return axis == 0 ? rows : cols;
  }
  std::array<std::ptrdiff_t, 2> strides() const { return {s1, s2}; }
};

void test_array_views() {
  std::vector<double> values(12);
  for (std::size_t k = 0; k < values.size(); ++k) values[k] = double(k);

  const auto row = make_view2d(values, 3, 4);
  check(row.extent1() == 3 && row.extent2() == 4, "view extents");
  check(row.stride1() == 4 && row.stride2() == 1, "row-major strides");
  check(row(1, 2) == 6.0, "row-major access");
  check(row.is_contiguous(), "row-major contiguous");

  ArrayView2D<const double> const_row = row;
  check(const_row(2, 3) == 11.0, "mutable-to-const conversion");
  const std::vector<double>& const_values = values;
  const auto const_vector_view = make_view2d(const_values, 3, 4);
  check(const_vector_view(2, 3) == 11.0, "const vector view");

  ArrayView2D<double> column(values.data(), 3, 4, 1, 3);
  check(column.stride1() == 1 && column.stride2() == 3, "column strides");
  check(column(1, 2) == values[1 + 2 * 3], "column-major access");
  check(column.is_contiguous(), "column-major contiguous");

  std::vector<double> big(2 * 3 * 4);
  for (std::size_t k = 0; k < big.size(); ++k) big[k] = double(k);
  const auto slice = make_view2d(big.data() + 12, 3, 4, 4, 1);
  check(slice(1, 2) == big[12 + 4 + 2], "contiguous slice access");
  check(slice.is_contiguous(), "fixed-first slice is contiguous");
  const auto noncontiguous = make_view2d(big.data() + 2, 2, 3, 12, 4);
  check(!noncontiguous.is_contiguous(), "fixed-last slice is noncontiguous");
  check(noncontiguous(1, 2) == big[2 + 12 + 8], "noncontiguous access");

  CustomMatrix2D custom;
  custom.values = values;
  custom.rows = 3;
  custom.cols = 4;
  custom.s1 = 4;
  custom.s2 = 1;
  const auto custom_view = make_view2d(custom);
  check(custom_view(2, 1) == 9.0, "container adapter access");
  const CustomMatrix2D& const_custom = custom;
  const auto const_custom_view = make_view2d(const_custom);
  check(const_custom_view(2, 1) == 9.0, "const container adapter access");

  expect_throw([] {
    std::vector<double> wrong(3);
    (void)make_view2d(wrong, 2, 2);
  }, "wrong vector size must throw");
}

void test_boundary_conditions() {
  Mesh mesh({0.0, 1.0, 3.0}, {0.0, 1.0, 2.0});
  BoundaryConditions bc;
  for (Side side : all_sides()) {
    check(bc.type(side) == BCType::ZeroFlux, "default is ZeroFlux");
  }
  bc.validate(mesh);

  bc.set_dirichlet(Side::Q1Min, 0.5);
  check(bc.type(Side::Q1Min) == BCType::Dirichlet, "Dirichlet type");
  check(bc.is_constant(Side::Q1Min), "constant Dirichlet");
  check(near(bc.value(Side::Q1Min, 0), 0.5), "constant value");
  bc.validate(mesh);

  std::vector<double> q1_values(mesh.n2() + 1, 0.5);
  std::vector<double> q2_values(mesh.n1() + 1, 0.5);
  bc.set_dirichlet(Side::Q1Min, q1_values);
  bc.set_dirichlet(Side::Q2Min, q2_values);
  check(!bc.is_constant(Side::Q1Min), "array-valued Dirichlet");
  check(near(bc.value(Side::Q1Min, 1), 0.5), "array value");
  bc.validate(mesh);

  BoundaryConditions conflict;
  conflict.set_dirichlet(Side::Q1Min, 1.0);
  conflict.set_dirichlet(Side::Q2Min, 2.0);
  expect_throw([&] { conflict.validate(mesh); }, "corner conflict must throw");

  BoundaryConditions wrong_length;
  std::vector<double> short_values(mesh.n2(), 0.5);
  wrong_length.set_dirichlet(Side::Q1Max, short_values);
  expect_throw([&] { wrong_length.validate(mesh); },
               "wrong Dirichlet length must throw");

  BoundaryConditions nonfinite;
  nonfinite.set_dirichlet(Side::Q1Min,
                          std::numeric_limits<double>::quiet_NaN());
  expect_throw([&] { nonfinite.validate(mesh); },
               "non-finite Dirichlet value must throw");

  BoundaryConditions mixed_corner;
  mixed_corner.set_dirichlet(Side::Q1Min, 1.0);
  mixed_corner.set_zero_flux(Side::Q2Min);
  mixed_corner.validate(mesh);
}

void test_solve_and_inplace() {
  Mesh mesh({0.0, 0.4, 1.0}, {0.0, 0.5, 1.1, 1.7});
  CaseConfig config;
  config.bc = reference::mixed_dirichlet_zeroflux();
  config.inv_tau = 0.0;
  config.off_diagonal = true;
  config.varying_coefficients = true;
  CaseData data(mesh, config);
  Solver solver(mesh);

  const std::size_t n1 = mesh.n1();
  const std::size_t n2 = mesh.n2();
  const std::size_t n = n1 * n2;
  const std::vector<double> f0 = data.f();
  std::vector<double> out_of_place(n, 0.0);

  solver.solve(data.equation(), data.f_const_view(),
               make_view2d(out_of_place, n1, n2), 0.01);
  check(max_abs_diff(f0, out_of_place) > 0.0, "solve must change the field");
  check(*std::min_element(out_of_place.begin(), out_of_place.end()) > 0.0,
        "solution must remain positive");
  check(max_abs_diff(f0, data.f()) == 0.0,
        "out-of-place solve must preserve the input storage");

  // In-place solve.
  std::vector<double> in_place = f0;
  solver.solve(data.equation(), make_view2d(in_place, n1, n2),
               make_view2d(in_place, n1, n2), 0.01);
  check(max_abs_diff(in_place, out_of_place) < 1e-14,
        "in-place solve must match out-of-place");

  // Column-major views.
  std::vector<double> column(n);
  for (std::size_t i = 0; i < n1; ++i) {
    for (std::size_t j = 0; j < n2; ++j) {
      column[i + j * n1] = f0[i * n2 + j];
    }
  }
  solver.solve(data.equation(), make_view2d(f0.data(), n1, n2, n2, 1),
               make_view2d(column.data(), n1, n2, 1, n1), 0.01);
  double column_diff = 0.0;
  for (std::size_t i = 0; i < n1; ++i) {
    for (std::size_t j = 0; j < n2; ++j) {
      column_diff = std::max(
          column_diff, std::abs(column[i + j * n1] - out_of_place[i * n2 + j]));
    }
  }
  check(column_diff < 1e-14, "column-major output must match");

  // Noncontiguous slice of a larger buffer: shape (n1, n2, 2), fixed k=1.
  std::vector<double> big(2 * n, -1.0);
  for (std::size_t i = 0; i < n1; ++i) {
    for (std::size_t j = 0; j < n2; ++j) {
      big[i * n2 * 2 + j * 2 + 1] = f0[i * n2 + j];
    }
  }
  auto big_slice = make_view2d(big.data() + 1, n1, n2, n2 * 2, 2);
  solver.solve(data.equation(), data.f_const_view(), big_slice, 0.01);
  for (std::size_t i = 0; i < n1; ++i) {
    for (std::size_t j = 0; j < n2; ++j) {
      check(near(big[i * n2 * 2 + j * 2 + 1], out_of_place[i * n2 + j], 1e-14),
            "strided slice result");
    }
  }

  // Repeated solve with the same solver and unchanged data.
  std::vector<double> repeated(n);
  solver.solve(data.equation(), make_view2d(f0.data(), n1, n2),
               make_view2d(repeated.data(), n1, n2), 0.01);
  check(max_abs_diff(repeated, out_of_place) < 1e-14,
        "repeated solve must be deterministic");

  // Update coefficients in place and solve again.
  data.G()[0] *= 1.1;
  data.D11()[0] *= 1.05;
  std::vector<double> updated(n);
  solver.solve(data.equation(), make_view2d(f0.data(), n1, n2),
               make_view2d(updated.data(), n1, n2), 0.01);
  check(max_abs_diff(updated, out_of_place) > 1e-13,
        "coefficient update must affect the solve");
  check(*std::min_element(updated.begin(), updated.end()) > 0.0,
        "updated solution must remain positive");
}

void test_conservation_and_positivity() {
  Mesh mesh({0.0, 0.3, 0.8, 1.0}, {0.0, 0.4, 0.9, 1.4, 1.8});
  CaseConfig config;
  config.bc = reference::all_zeroflux();
  config.inv_tau = 0.0;
  config.off_diagonal = true;
  config.varying_coefficients = true;
  CaseData data(mesh, config);
  Solver solver(mesh);
  const double mass0 = data.total_mass();
  solver.solve(data.equation(), data.f_const_view(), data.f_view(), 0.01);
  const double mass1 = data.total_mass();
  check(std::abs(mass1 - mass0) <= 1e-12 * std::abs(mass0),
        "ZeroFlux mass conservation");
  check(*std::min_element(data.f().begin(), data.f().end()) > -1e-13,
        "conservative solve positivity");

  // A two-cell mesh has exactly one interior face. With ZeroFlux boundaries
  // and no source/loss, total-mass conservation directly checks that the
  // paired interior-face flux cancels.
  Mesh one_face_mesh({0.0, 0.5, 1.0}, {0.0, 1.0});
  CaseConfig one_face_config;
  one_face_config.bc = reference::all_zeroflux();
  one_face_config.inv_tau = 0.0;
  one_face_config.off_diagonal = true;
  one_face_config.varying_coefficients = true;
  CaseData one_face(one_face_mesh, one_face_config);
  Solver one_face_solver(one_face_mesh);
  const double one_face_mass0 = one_face.total_mass();
  one_face_solver.solve(one_face.equation(), one_face.f_const_view(),
                        one_face.f_view(), 0.01);
  const double one_face_mass1 = one_face.total_mass();
  check(std::abs(one_face_mass1 - one_face_mass0) <=
            1e-12 * std::abs(one_face_mass0),
        "single interior face conservation");

  CaseConfig loss_config = config;
  loss_config.bc = reference::mixed_dirichlet_zeroflux();
  loss_config.inv_tau = 0.5;
  CaseData loss_data(mesh, loss_config);
  const double loss_mass0 = loss_data.total_mass();
  solver.solve(loss_data.equation(), loss_data.f_const_view(),
               loss_data.f_view(), 0.01);
  check(loss_data.total_mass() < loss_mass0, "positive loss must reduce mass");
  check(*std::min_element(loss_data.f().begin(), loss_data.f().end()) > 0.0,
        "loss solution positivity");
}

void test_source_branch() {
  Mesh mesh({0.0, 0.25, 0.5, 0.75, 1.0}, {0.0, 0.3, 0.7, 1.0});
  CaseConfig source_config;
  source_config.bc = reference::all_zeroflux();
  source_config.inv_tau = -0.1;
  source_config.off_diagonal = false;
  source_config.varying_coefficients = false;
  CaseData source(mesh, source_config);
  Solver solver(mesh);
  const double source_mass0 = source.total_mass();
  solver.solve(source.equation(), source.f_const_view(), source.f_view(),
               0.01);
  const double source_mass1 = source.total_mass();
  check(source_mass1 > source_mass0,
        "negative inverse loss must add mass (source)");
  check(*std::min_element(source.f().begin(), source.f().end()) > 0.0,
        "source solution positivity");

  CaseConfig zero_config = source_config;
  zero_config.inv_tau = 0.0;
  CaseData zero(mesh, zero_config);
  solver.solve(zero.equation(), zero.f_const_view(), zero.f_view(), 0.01);
  check(source_mass1 > zero.total_mass(),
        "source case must have more mass than the no-source case");
}

void test_repeated_boundary_updates() {
  Mesh mesh({0.0, 0.3, 0.7, 1.0}, {0.0, 0.4, 0.9, 1.4, 1.8});
  CaseConfig config;
  config.bc = reference::mixed_dirichlet_zeroflux();
  config.inv_tau = 0.0;
  config.off_diagonal = false;
  config.varying_coefficients = true;
  CaseData data(mesh, config);
  Solver solver(mesh);
  const std::size_t n1 = mesh.n1();
  const std::size_t n2 = mesh.n2();
  const std::size_t n = n1 * n2;
  const std::vector<double> f0 = data.f();

  std::vector<double> first(n);
  solver.solve(data.equation(), make_view2d(f0.data(), n1, n2),
               make_view2d(first.data(), n1, n2), 0.01);

  data.set_boundary_value(Side::Q1Min, 1, 0.75);
  std::vector<double> second(n);
  solver.solve(data.equation(), make_view2d(f0.data(), n1, n2),
               make_view2d(second.data(), n1, n2), 0.01);
  check(max_abs_diff(first, second) > 1e-13,
        "in-place boundary update must affect the solve");

  // Reuse the same solver with a different boundary type. Both supported
  // types produce the same matrix sparsity pattern, so symbolic analysis is
  // still valid.
  BoundaryConditions all_dirichlet;
  for (Side side : all_sides()) all_dirichlet.set_dirichlet(side, 0.5);
  Equation dirichlet_equation(
      make_view2d(data.G(), n1, n2), make_view2d(data.D11(), n1, n2),
      make_view2d(data.D12(), n1, n2), make_view2d(data.D22(), n1, n2),
      all_dirichlet);
  std::vector<double> third(n);
  solver.solve(dirichlet_equation, make_view2d(f0.data(), n1, n2),
               make_view2d(third.data(), n1, n2), 0.01);
  check(*std::min_element(third.begin(), third.end()) > 0.0,
        "solver reuse with a different BC type");
}

void test_boundary_reference_semantics() {
  // Equation stores a reference, so descriptor/type/span changes made after
  // construction must be visible on the next solve.
  Mesh mesh({0.0, 0.4, 1.0}, {0.0, 0.5, 1.1, 1.7});
  CaseConfig config;
  config.bc = reference::all_zeroflux();
  config.inv_tau = 0.0;
  config.off_diagonal = false;
  config.varying_coefficients = false;
  CaseData data(mesh, config);
  Solver solver(mesh);

  const std::size_t n1 = mesh.n1();
  const std::size_t n2 = mesh.n2();
  const std::size_t n = n1 * n2;
  const std::vector<double> f0 = data.f();

  auto solve_to = [&](std::vector<double>& out) {
    solver.solve(data.equation(), make_view2d(f0.data(), n1, n2),
                 make_view2d(out.data(), n1, n2), 0.01);
  };

  std::vector<double> all_zero_flux(n);
  solve_to(all_zero_flux);

  // ZeroFlux -> Dirichlet after Equation construction.
  data.bc().set_dirichlet(Side::Q1Min, 0.0);
  std::vector<double> dirichlet_zero(n);
  solve_to(dirichlet_zero);
  check(max_abs_diff(all_zero_flux, dirichlet_zero) > 1e-13,
        "descriptor type change must affect solve");

  // Modify a constant Dirichlet value after construction.
  data.bc().set_dirichlet(Side::Q1Min, 0.5);
  std::vector<double> dirichlet_half(n);
  solve_to(dirichlet_half);
  check(max_abs_diff(dirichlet_zero, dirichlet_half) > 1e-13,
        "Dirichlet value change must affect solve");

  // Replace the boundary data span after construction.
  std::vector<double> replacement(n2 + 1, 0.2);
  data.bc().set_dirichlet(Side::Q1Min, replacement);
  std::vector<double> span_replaced(n);
  solve_to(span_replaced);
  check(max_abs_diff(dirichlet_half, span_replaced) > 1e-13,
        "boundary span replacement must affect solve");

  // Invalid length introduced after construction is rejected next solve.
  std::vector<double> wrong_length(n2, 0.1);
  data.bc().set_dirichlet(Side::Q1Max, wrong_length);
  expect_throw([&] { solve_to(span_replaced); },
               "invalid boundary configuration after construction must throw");

  // A temporary BoundaryConditions must not be accepted. The deleted
  // rvalue overload makes this construction ill-formed.
  static_assert(!std::is_constructible_v<
                    Equation, ArrayView2D<const double>,
                    ArrayView2D<const double>, ArrayView2D<const double>,
                    ArrayView2D<const double>, BoundaryConditions&&>,
                "Equation must reject temporary BoundaryConditions");
}


void test_eigen_support() {
  Mesh mesh({0.0, 0.4, 1.0, 1.5}, {0.0, 0.5, 1.2});
  const std::size_t n1 = mesh.n1();
  const std::size_t n2 = mesh.n2();

  Eigen::MatrixXd G(n1, n2);
  Eigen::MatrixXd D11(n1, n2);
  Eigen::MatrixXd D12(n1, n2);
  Eigen::MatrixXd D22(n1, n2);
  Eigen::MatrixXd inv_tau(n1, n2);
  for (std::size_t i = 0; i < n1; ++i) {
    for (std::size_t j = 0; j < n2; ++j) {
      const auto ii = static_cast<Eigen::Index>(i);
      const auto jj = static_cast<Eigen::Index>(j);
      G(ii, jj) = 1.0 + 0.1 * mesh.q1_center(i);
      D11(ii, jj) = 1.0 + 0.2 * mesh.q2_center(j);
      D12(ii, jj) = 0.05;
      D22(ii, jj) = 0.8 + 0.1 * mesh.q1_center(i);
      inv_tau(ii, jj) = 0.02;
    }
  }

  BoundaryConditions bc;
  bc.set_dirichlet(Side::Q1Min, 0.0);
  bc.set_zero_flux(Side::Q1Max);
  bc.set_zero_flux(Side::Q2Min);
  bc.set_zero_flux(Side::Q2Max);

  // Direct Eigen coefficients. Inverse loss defaults to zero and is set
  // separately when needed.
  Equation eigen_equation(G, D11, D12, D22, bc);
  check(!eigen_equation.has_loss(), "Eigen inverse loss defaults to zero");
  check(eigen_equation.inv_tau(0, 0) == 0.0, "default inverse-loss value");
  eigen_equation.set_loss(inv_tau);
  check(eigen_equation.has_loss() && !eigen_equation.loss_is_constant(),
        "Eigen inverse-loss matrix set manually");
  check(eigen_equation.G_view().data() == G.data(),
        "Eigen coefficient view must be zero-copy");
  check(eigen_equation.G_view().extent1() == n1 &&
            eigen_equation.G_view().extent2() == n2,
        "Eigen coefficient extents");
  check(eigen_equation.G_view().stride1() == 1 &&
            eigen_equation.G_view().stride2() ==
                static_cast<std::ptrdiff_t>(n1),
        "Eigen column-major strides");

  const Eigen::MatrixXd G_const = G;
  Equation const_equation(G_const, D11, D12, D22, bc);
  check(const_equation.G_view().data() == G_const.data(),
        "const Eigen coefficient matrices are supported");
  // Same data through the existing explicit view interface.
  ArrayView2D<const double> Gv(G.data(), n1, n2, 1, n1);
  ArrayView2D<const double> D11v(D11.data(), n1, n2, 1, n1);
  ArrayView2D<const double> D12v(D12.data(), n1, n2, 1, n1);
  ArrayView2D<const double> D22v(D22.data(), n1, n2, 1, n1);
  ArrayView2D<const double> tau_v(inv_tau.data(), n1, n2, 1, n1);
  Equation view_equation(Gv, D11v, D12v, D22v, bc);
  view_equation.set_loss(tau_v);

  Eigen::MatrixXd f0(n1, n2);
  for (std::size_t i = 0; i < n1; ++i) {
    for (std::size_t j = 0; j < n2; ++j) {
      f0(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
          1.0 + 0.1 * static_cast<double>(i + j);
    }
  }

  Solver solver(mesh);
  Eigen::MatrixXd eigen_out(n1, n2);
  Eigen::MatrixXd view_out(n1, n2);
  solver.solve(eigen_equation, f0, eigen_out, 0.01);
  solver.solve(view_equation, make_view2d(f0.data(), n1, n2, 1, n1),
               make_view2d(view_out.data(), n1, n2, 1, n1), 0.01);
  check(max_abs_diff(eigen_out, view_out) < 1e-14,
        "Eigen and view solves must agree");

  // In-place solve using the same Eigen matrix as input and output.
  Eigen::MatrixXd inplace = f0;
  solver.solve(eigen_equation, inplace, inplace, 0.01);
  check(max_abs_diff(eigen_out, inplace) < 1e-14,
        "Eigen in-place solve must agree");

  // Coefficient updates after Equation construction are visible.
  G(1, 1) += 0.25;
  D11(1, 1) += 0.10;
  Eigen::MatrixXd updated(n1, n2);
  solver.solve(eigen_equation, f0, updated, 0.01);
  check(max_abs_diff(eigen_out, updated) > 1e-12,
        "Eigen coefficient updates must affect solve");

  // Wrong field sizes are rejected by the existing Solver validation.
  Eigen::MatrixXd wrong(n1 + 1, n2);
  expect_throw([&] { solver.solve(eigen_equation, f0, wrong, 0.01); },
               "wrong Eigen output size must throw");
  expect_throw([&] { solver.solve(eigen_equation, wrong, eigen_out, 0.01); },
               "wrong Eigen input size must throw");

  // Mismatched coefficient extents are rejected by Equation.
  Eigen::MatrixXd bad_D11(n1, n2 + 1);
  expect_throw(
      [&] { Equation bad(G, bad_D11, D12, D22, bc); },
      "mismatched Eigen coefficient extents must throw");

  // Temporary coefficient matrices must not be accepted. The lvalue-only
  // constructor makes the four-matrix form ill-formed, and there is no
  // five-matrix Equation constructor.
  static_assert(!std::is_constructible_v<
                    Equation, Eigen::MatrixXd, Eigen::MatrixXd,
                    Eigen::MatrixXd, Eigen::MatrixXd, BoundaryConditions&>,
                "Equation must reject temporary Eigen coefficients");
  static_assert(!std::is_constructible_v<
                    Equation, Eigen::MatrixXd, Eigen::MatrixXd,
                    Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd,
                    BoundaryConditions&>,
                "Equation must not have an inv_tau constructor");
}

void test_albert_young_loader() {
  const std::filesystem::path source =
      "data/albert_young/ALBERT_YOUNG_80x80.bin";
  const auto data = albert_young::load_data(source.string());
  check(data.n1 == 80 && data.n2 == 80, "loader dimensions");
  check(data.q1_edges.size() == 81 && data.q2_edges.size() == 81,
        "loader edge sizes");
  check(data.G.size() == 6400 && data.D11.size() == 6400 &&
            data.D12.size() == 6400 && data.D22.size() == 6400,
        "loader coefficient sizes");

  std::ifstream input(source, std::ios::binary);
  check(static_cast<bool>(input), "cannot read reference binary");
  const std::vector<char> bytes((std::istreambuf_iterator<char>(input)),
                                std::istreambuf_iterator<char>());
  check(!bytes.empty(), "reference binary is empty");

  const auto temporary =
      std::filesystem::temp_directory_path() /
      "sayram2d_albert_young_loader_test.bin";
  auto write_temporary = [&](const std::vector<char>& content) {
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    check(static_cast<bool>(output), "cannot open loader test file");
    output.write(content.data(), static_cast<std::streamsize>(content.size()));
    check(static_cast<bool>(output), "cannot write loader test file");
  };
  auto expect_bad_file = [&](const std::vector<char>& content,
                             const char* message) {
    write_temporary(content);
    expect_throw(
        [&] { (void)albert_young::load_data(temporary.string()); }, message);
  };

  write_temporary(bytes);
  const auto round_trip = albert_young::load_data(temporary.string());
  check(round_trip.n1 == data.n1 && round_trip.n2 == data.n2 &&
            round_trip.G == data.G && round_trip.D22 == data.D22,
        "loader round trip");

  auto malformed = bytes;
  malformed[0] = 'X';
  expect_bad_file(malformed, "bad magic must be rejected");

  malformed = bytes;
  malformed.resize(malformed.size() - 1);
  expect_bad_file(malformed, "truncated binary must be rejected");

  malformed = bytes;
  malformed.push_back('\0');
  expect_bad_file(malformed, "trailing binary data must be rejected");

  std::filesystem::remove(temporary);
  expect_throw([] { (void)albert_young::load_data("missing.bin"); },
               "missing binary must be rejected");
}

void test_strided_3d_slice_integration() {
  // Two interleaved layers of a 3-D buffer, accessed through explicit
  // arbitrary-stride ArrayView2D slices. This is the caller-controlled path
  // for slices; the core never infers layout from a container.
  Mesh mesh({0.0, 0.4, 1.0}, {0.0, 0.5, 1.2});
  const std::size_t n1 = mesh.n1();
  const std::size_t n2 = mesh.n2();
  const std::size_t n = n1 * n2;
  const std::size_t layers = 2;

  std::vector<double> G(layers * n);
  std::vector<double> D11(layers * n);
  std::vector<double> D12(layers * n);
  std::vector<double> D22(layers * n);
  std::vector<double> f(layers * n);
  std::vector<double> out(layers * n, 0.0);

  auto slice = [&](std::vector<double>& values, std::size_t layer) {
    return make_view2d(values.data() + layer, n1, n2,
                       static_cast<std::ptrdiff_t>(layers * n2),
                       static_cast<std::ptrdiff_t>(layers));
  };

  for (std::size_t i = 0; i < n1; ++i) {
    for (std::size_t j = 0; j < n2; ++j) {
      for (std::size_t layer = 0; layer < layers; ++layer) {
        const std::size_t index = (i * n2 + j) * layers + layer;
        G[index] = 1.0 + 0.1 * static_cast<double>(layer);
        D11[index] = 1.0;
        D12[index] = 0.0;
        D22[index] = 0.9 + 0.1 * static_cast<double>(layer);
        f[index] = 1.0 + static_cast<double>(layer) + 0.1 * (i + j);
      }
    }
  }

  BoundaryConditions bc;
  bc.set_dirichlet(Side::Q1Min, 0.0);
  bc.set_zero_flux(Side::Q1Max);
  bc.set_zero_flux(Side::Q2Min);
  bc.set_zero_flux(Side::Q2Max);

  Solver solver(mesh);
  Equation equation0(slice(G, 0), slice(D11, 0), slice(D12, 0),
                     slice(D22, 0), bc);
  Equation equation1(slice(G, 1), slice(D11, 1), slice(D12, 1),
                     slice(D22, 1), bc);

  solver.solve(equation0, slice(f, 0), slice(out, 0), 0.01);
  solver.solve(equation1, slice(f, 1), slice(out, 1), 0.01);

  // Equivalent standalone layer-0 problem.
  std::vector<double> G0(n), D110(n), D120(n), D220(n), f0(n), out0(n);
  for (std::size_t k = 0; k < n; ++k) {
    G0[k] = G[k * layers];
    D110[k] = D11[k * layers];
    D120[k] = D12[k * layers];
    D220[k] = D22[k * layers];
    f0[k] = f[k * layers];
  }
  Equation standalone(make_view2d(G0, n1, n2), make_view2d(D110, n1, n2),
                     make_view2d(D120, n1, n2), make_view2d(D220, n1, n2),
                     bc);
  Solver standalone_solver(mesh);
  standalone_solver.solve(standalone, make_view2d(f0, n1, n2),
                          make_view2d(out0, n1, n2), 0.01);

  double baseline_diff = 0.0;
  double layer_diff = 0.0;
  for (std::size_t k = 0; k < n; ++k) {
    baseline_diff =
        std::max(baseline_diff, std::abs(out[k * layers] - out0[k]));
    layer_diff = std::max(
        layer_diff, std::abs(out[k * layers] - out[k * layers + 1]));
  }
  check(baseline_diff < 1e-14,
        "strided slice solve must agree with standalone 2-D solve");
  check(layer_diff > 1e-12, "3-D slices must remain isolated");

  auto out0_view = slice(out, 0);
  auto out1_view = slice(out, 1);
  check(out0_view.data() == out.data() && out1_view.data() == out.data() + 1,
        "slice data pointers must alias the 3-D buffer");
  check(out0_view.stride1() == static_cast<std::ptrdiff_t>(layers * n2) &&
            out0_view.stride2() == static_cast<std::ptrdiff_t>(layers),
        "slice strides must be explicit");

  // Change a coefficient and dt after construction, without rebuilding.
  G[0] += 0.25;
  std::vector<double> changed(layers * n, 0.0);
  solver.solve(equation0, slice(f, 0), slice(changed, 0), 0.02);
  check(std::abs(changed[0] - out[0]) > 1e-12,
        "slice coefficient/dt update must affect solve");

  // Boundary descriptor changes remain visible through the reference.
  bc.set_dirichlet(Side::Q1Min, 0.5);
  std::vector<double> boundary_changed(layers * n, 0.0);
  solver.solve(equation0, slice(f, 0), slice(boundary_changed, 0), 0.01);
  check(std::abs(boundary_changed[0] - changed[0]) > 1e-12,
        "slice boundary update must affect solve");
}

void test_physics_helpers() {
  check(near(physics::degrees_to_radians(180.0), physics::pi, 1e-15),
        "degree-to-radian conversion");
  check(near(physics::radians_to_degrees(physics::pi), 180.0, 1e-13),
        "radian-to-degree conversion");
  check(near(physics::joule_to_mev(physics::mev_to_joule(2.0)), 2.0, 1e-14),
        "MeV-to-joule round trip");
  check(near(physics::joule_to_kev(physics::kev_to_joule(2.0)), 2.0, 1e-14),
        "keV-to-joule round trip");

  const double rest_j = physics::ev_to_joule(0.511e6);
  const double kinetic_j = physics::mev_to_joule(1.0);
  const double momentum =
      physics::momentum_si(kinetic_j, rest_j);
  const double kinetic_back =
      physics::kinetic_energy_si(momentum, rest_j);
  check(momentum > 0.0, "SI momentum is positive");
  check(near(kinetic_back, kinetic_j, 1e-14),
        "momentum/kinetic-energy round trip");
  check(near(physics::momentum_mev_c(0.0), 0.0, 1e-15),
        "zero kinetic energy has zero momentum");
}

void test_validation_errors() {
  Mesh mesh({0.0, 0.4, 1.0}, {0.0, 0.5, 1.1});
  CaseConfig config;
  config.bc = reference::mixed_dirichlet_zeroflux();
  config.inv_tau = 0.0;
  config.off_diagonal = false;
  config.varying_coefficients = false;
  CaseData data(mesh, config);
  Solver solver(mesh);
  const std::size_t n1 = mesh.n1();
  const std::size_t n2 = mesh.n2();
  const std::size_t n = n1 * n2;

  std::vector<double> wrong((n1 + 1) * n2, 1.0);
  std::vector<double> out(n);
  expect_throw([&] {
    solver.solve(data.equation(), make_view2d(wrong, n1 + 1, n2),
                 make_view2d(out, n1, n2), 0.01);
  }, "wrong input extents must throw");

  expect_throw([&] {
    solver.solve(data.equation(), data.f_const_view(), data.f_view(), 0.0);
  }, "zero dt must throw");
  expect_throw([&] {
    solver.solve(data.equation(), data.f_const_view(), data.f_view(), -0.01);
  }, "negative dt must throw");

  data.G()[0] = 0.0;
  expect_throw([&] {
    solver.solve(data.equation(), data.f_const_view(), data.f_view(), 0.01);
  }, "zero G must throw");
  data.G()[0] = 1.0;

  data.D11()[0] = 0.0;
  expect_throw([&] {
    solver.solve(data.equation(), data.f_const_view(), data.f_view(), 0.01);
  }, "non-positive-definite D must throw");
  data.D11()[0] = 1.0;

  data.f()[0] = std::numeric_limits<double>::quiet_NaN();
  expect_throw([&] {
    solver.solve(data.equation(), data.f_const_view(), data.f_view(), 0.01);
  }, "non-finite input must throw");
  data.f()[0] = 1.0;

  expect_throw([&] {
    Equation bad(make_view2d(data.G(), n1, n2),
                 make_view2d(data.D11(), n1 + 1, n2),
                 make_view2d(data.D12(), n1, n2),
                 make_view2d(data.D22(), n1, n2), data.equation().bc());
  }, "mismatched coefficient extents must throw");

  std::vector<double> loss((n1 + 1) * n2, 0.0);
  Equation bad(make_view2d(data.G(), n1, n2),
               make_view2d(data.D11(), n1, n2),
               make_view2d(data.D12(), n1, n2),
               make_view2d(data.D22(), n1, n2), data.equation().bc());
  expect_throw([&] {
    bad.set_loss(make_view2d(loss, n1 + 1, n2));
  }, "mismatched loss extents must throw");
}

int run(const char* name, void (*function)()) {
  try {
    function();
    std::cout << "[PASS] " << name << "\n";
    return 0;
  } catch (const std::exception& error) {
    std::cout << "[FAIL] " << name << ": " << error.what() << "\n";
    return 1;
  }
}

}  // namespace

int main() {
  int failures = 0;
  failures += run("mesh_geometry", test_mesh_geometry);
  failures += run("array_views", test_array_views);
  failures += run("boundary_conditions", test_boundary_conditions);
  failures += run("solve_and_inplace", test_solve_and_inplace);
  failures += run("conservation_and_positivity", test_conservation_and_positivity);
  failures += run("source_branch", test_source_branch);
  failures += run("repeated_boundary_updates", test_repeated_boundary_updates);
  failures += run("boundary_reference_semantics", test_boundary_reference_semantics);
  failures += run("eigen_support", test_eigen_support);
  failures += run("physics_helpers", test_physics_helpers);
  failures += run("albert_young_loader", test_albert_young_loader);
  failures += run("strided_3d_slice_integration",
                  test_strided_3d_slice_integration);
  failures += run("validation_errors", test_validation_errors);
  if (failures == 0) {
    std::cout << "All tests passed.\n";
  }
  return failures == 0 ? 0 : 1;
}
