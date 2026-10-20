/*
 * File:        Solver.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        01/13/2025
 *
 * Copyright (c) Xin Tao
 *
 */

#include "Solver.h"
#include "BCTypes.h"

Solver::Solver(const Mesh& m_in, Equation* eqp)
  : m(m_in), eq(*eqp){

    std::size_t nx = m.nx();
    std::size_t ny = m.ny();

    M_.resize(nx*ny, nx*ny);

    f_.resize({nx,ny});

    R_.resize(nx*ny);
    ftmp_.resize(nx*ny);

    Lambda_.resize({nx, ny});
    vertex_f_.resize({nx+1,ny+1});

    istep_ = 0;

    init();

  }

void Solver::init(){

  for (std::size_t i = 0; i < m.nx(); i++){
    for (std::size_t j = 0; j < m.ny(); j++){
      f_(i,j) = eq.init_f({i,j});
    }
  }

  update_Lambda();
  update_vertex_f();

  // One-time symbolic analysis for SparseLU (pattern only).
  // Assumes the sparsity pattern will not change over time steps.
  R_.setZero();
  M_coeffs_.clear();
  assemble();

  M_.makeCompressed();
  eigen_solver.analyzePattern(M_);
}

void Solver::update_Lambda(){

  for (std::size_t i = 0; i < m.nx(); i++){
    for (std::size_t j = 0; j < m.ny(); j++){
      Lambda_(i,j) << eq.Dxx({i,j}) * eq.G({i,j}), eq.Dxy({i,j}) * eq.G({i,j}),
                      eq.Dxy({i,j}) * eq.G({i,j}), eq.Dyy({i,j}) * eq.G({i,j});
    }
  }
}


void Solver::a_sigma_func(const Ind& ind, int inbr, Vector2* a_sigma_i_p, double* a_sigmap){

  Edge edge;
  m.get_nbr_edge(ind, inbr, &edge);

  Point K = {m.x(ind.i), m.y(ind.j)};

  Point A = edge.v(0);
  Point B = edge.v(1);

  Vector2 vkb = B - K;  // as Vector KB (from K point to B)
  Vector2 vka = A - K;

  Vector2 rvkb = {vkb(1), -vkb(0)};  // (x,y) rotated by 90 clockwise, becoming (y, -x)
  Vector2 rvka = {vka(1), -vka(0)};

  double a_sigma_a = edge.length() * (edge.n().transpose() * Lambda(ind) * rvkb)(0) / (vka.transpose() * rvkb)(0);
  double a_sigma_b = edge.length() * (edge.n().transpose() * Lambda(ind) * rvka)(0) / (vkb.transpose() * rvka)(0);

  (*a_sigma_i_p)(0) = a_sigma_a;
  (*a_sigma_i_p)(1) = a_sigma_b;

  const VtxInd a = edge.vind(0);
  const VtxInd b = edge.vind(1);

  const double fA = vertex_f(a.i, a.j);
  const double fB = vertex_f(b.i, b.j);

  *a_sigmap = a_sigma_a * fA + a_sigma_b * fB;
}

void Solver::apply_inner_face_pair(const Ind& ind, int inbr){ // add coefficient from a inner (no-a-boundary) face
  //
  Ind nbr_ind;

  int rinb;

  double a_sigma_K, a_sigma_L;
  double mu_K, mu_L;

  double B_sigma, B_sigma_abs, B_sigma_plus, B_sigma_minus;
  double A_K, A_L;

  Vector2 a_sigma_i_K;
  Vector2 a_sigma_i_L;

  // calculate A_K
  a_sigma_func(ind, inbr, &a_sigma_i_K, &a_sigma_K);

  // calculate A_L
  m.get_nbr_ind(ind, inbr, &nbr_ind);

  rinb = m.rinbr(inbr);
  a_sigma_func(nbr_ind, rinb, &a_sigma_i_L, &a_sigma_L);

  calculate_mu(a_sigma_K, a_sigma_L, &mu_K, &mu_L);

  B_sigma = mu_L*a_sigma_L - mu_K*a_sigma_K;

  B_sigma_abs = std::abs(B_sigma);
  B_sigma_plus = (B_sigma_abs + B_sigma) / 2.0;
  B_sigma_minus = (B_sigma_abs - B_sigma) / 2.0;

  A_K = mu_K * (a_sigma_i_K[0] + a_sigma_i_K[1]) + B_sigma_plus/(f(ind)+gEPS);
  A_L = mu_L * (a_sigma_i_L[0] + a_sigma_i_L[1]) + B_sigma_minus/(f(nbr_ind)+gEPS);

  std::size_t K = m.flatten_cell_index(ind), L = m.flatten_cell_index(nbr_ind);

  M_coeffs_.push_back(T(K, K, A_K));
  M_coeffs_.push_back(T(K, L, -A_L));
  M_coeffs_.push_back(T(L, L, A_L));
  M_coeffs_.push_back(T(L, K, -A_K));

}

void Solver::apply_dirichlet_face(const Ind& ind, int inbr){ // add coefficient from a Dirichlet face
  // j
  double a_sigma_K;

  double B_sigma, B_sigma_abs, B_sigma_plus, B_sigma_minus;
  double A_K;
  Vector2 a_sigma_i_K;

  // calculate A_K
  a_sigma_func(ind, inbr, &a_sigma_i_K, &a_sigma_K);

  B_sigma = -a_sigma_K;
  B_sigma_abs = std::abs(B_sigma);
  B_sigma_plus = (B_sigma_abs + B_sigma) / 2.0;
  B_sigma_minus = (B_sigma_abs - B_sigma) / 2.0;

  A_K = a_sigma_i_K[0] + a_sigma_i_K[1] + B_sigma_plus/(f(ind)+gEPS);

  std::size_t K = m.flatten_cell_index(ind);
  M_coeffs_.push_back(T(K, K, A_K));
  R_(K) += B_sigma_minus;
}


void Solver::assemble(){ // obtain M and R

  std::size_t i, j;

  // 1) Interior face pairs (branch-free)
  //    west faces for i>=1
  for (i=1; i<m.nx(); ++i)
    for (j=0; j<m.ny(); ++j)
      apply_inner_face_pair({i, j}, m.inbr_im());
  
  //    south faces for j>=1
  for (i = 0; i < m.nx(); ++i) {
    for (j = 1; j < m.ny(); ++j) {
      apply_inner_face_pair({i, j}, m.inbr_jm());
    }
  }

  // 2) Boundary faces 
  apply_boundary_faces();

  long K;
  double UKK;

  for (std::size_t i=0; i<m.nx(); ++i) {
    for (std::size_t j=0; j<m.ny(); ++j) {
      K = m.flatten_cell_index({i,j});
      UKK = eq.G({i,j}) * m.cell_area_dt({i,j});

      M_coeffs_.push_back(T(K, K, UKK * (1 + m.dt() * eq.inv_tau({i,j}))));

      R_(K) += UKK * f_(i,j);
    }
  }

  M_.setFromTriplets(M_coeffs_.begin(), M_coeffs_.end());
}

void Solver::apply_boundary_faces() {

  // -------------------------
  // XMIN boundary: Face::IM on cells i = 0
  // -------------------------
  {
    const auto side = BoundaryID::XMIN;
    const BCType t = eq.bc_type(side);

    if (t == BCType::Dirichlet) {
      for (std::size_t j = 0; j < m.ny(); ++j) {
        apply_dirichlet_face({0, j}, m.inbr_im());
      }
    }     

  }

  // -------------------------
  // XMAX boundary: Face::IP on cells i = nx - 1
  // -------------------------
  {

    const auto side = BoundaryID::XMAX;
    const BCType t = eq.bc_type(side);

    if (t == BCType::Dirichlet) {
      for (std::size_t j = 0; j < m.ny(); ++j) {
        apply_dirichlet_face({m.nx() - 1, j}, m.inbr_ip());
      }
    } 

  }

  // -------------------------
  // YMIN boundary: Face::JM on cells j = 0
  // -------------------------
  {
    const auto side = BoundaryID::YMIN;
    const BCType t = eq.bc_type(side);

    if (t == BCType::Dirichlet) {
      for (std::size_t i = 0; i < m.nx(); ++i) {
        apply_dirichlet_face({i, 0}, m.inbr_jm());
      }
    } 

  }

  // -------------------------
  // YMAX boundary: Face::JP on cells j = ny - 1
  // -------------------------
  {

    const auto side = BoundaryID::YMAX; 
    const BCType t = eq.bc_type(side);

    if (t == BCType::Dirichlet) {
      for (std::size_t i = 0; i < m.nx(); ++i) {
        apply_dirichlet_face({i, m.ny() - 1}, m.inbr_jp());
      }
    } 
  }

}


void Solver::update() {
  R_.setZero();
  M_coeffs_.clear();

  assemble();

  M_.makeCompressed();
  eigen_solver.factorize(M_);
  ftmp_ = eigen_solver.solve(R_);
 
  for (std::size_t i=0; i<m.nx(); ++i){
    for (std::size_t j=0; j<m.ny(); ++j) {
      f_(i,j) = ftmp_(m.flatten_cell_index({i, j}));
    }
  }

  istep_ += 1;
  eq.update(t());
  update_Lambda();
  update_vertex_f();
}

void Solver::update_vertex_f(){

  fill_vertex_from_cells();
  fill_vertex_from_bcs();
}


void Solver::fill_vertex_from_cells() {
  const std::size_t nx = m.nx();
  const std::size_t ny = m.ny();


  auto lerp = [](double x, double x0, double x1, double f0, double f1) {
    const double denom = (x1 - x0);
    if (denom == 0.0) return 0.5 * (f0 + f1); // defensive
    const double w0 = (x1 - x) / denom;
    const double w1 = (x - x0) / denom;
    return w0 * f0 + w1 * f1;
  };


  // -----------------------------
  // 1) Corners (single adjacent cell)
  // -----------------------------
  vertex_f_(0,   0)   = f_(0,     0);
  vertex_f_(nx,  0)   = f_(nx-1,  0);
  vertex_f_(0,   ny)  = f_(0,     ny-1);
  vertex_f_(nx,  ny)  = f_(nx-1,  ny-1);

  // -----------------------------
  // 2) Interior vertices (bilinear)
  // i = 1..nx-1, j = 1..ny-1
  // -----------------------------

  for (std::size_t i=1; i<m.nx(); ++i){
    for (std::size_t j=1; j<m.ny(); ++j){

      const double xv = m.x_edge(i);
      const double yv = m.y_edge(j);

      const double xL = m.x(i - 1);
      const double xR = m.x(i);
      const double yB = m.y(j - 1);
      const double yT = m.y(j);

      const double denom_x = (xR - xL);
      const double denom_y = (yT - yB);

      const double wxL = (xR - xv) / denom_x;
      const double wxR = (xv - xL) / denom_x;
      const double wyB = (yT - yv) / denom_y;
      const double wyT = (yv - yB) / denom_y;

      // just in case the grid is nonuniform
      vertex_f_(i,j) = wxL * wyB * f_(i-1,j-1) + wxR * wyB * f_(i,j - 1) + wxL * wyT * f_(i-1,j) + wxR * wyT * f_(i,j);

    }
  }

  // -----------------------------
  // 3) Vertical boundary edges: i = 0 and i = nx (exclude corners)
  // Use linear interpolation in y between the two adjacent cells.
  // -----------------------------
  for (std::size_t j = 1; j < ny; ++j) {
    const double yv = m.y_edge(j);
    const double y0 = m.y(j - 1);
    const double y1 = m.y(j);

    // Left edge (i=0): adjacent cells are (0,j-1) and (0,j)
    vertex_f_(0, j) = lerp(yv, y0, y1, f_(0, j - 1), f_(0, j));

    // Right edge (i=nx): adjacent cells are (nx-1,j-1) and (nx-1,j)
    vertex_f_(nx, j) = lerp(yv, y0, y1, f_(nx - 1, j - 1), f_(nx - 1, j));
  }

  // -----------------------------
  // 4) Horizontal boundary edges: j = 0 and j = ny (exclude corners)
  // Use linear interpolation in x between the two adjacent cells.
  // -----------------------------
  for (std::size_t i = 1; i < nx; ++i) {
    const double xv = m.x_edge(i);
    const double x0 = m.x(i - 1);
    const double x1 = m.x(i);

    // Bottom edge (j=0): adjacent cells are (i-1,0) and (i,0)
    vertex_f_(i, 0) = lerp(xv, x0, x1, f_(i - 1, 0), f_(i, 0));

    // Top edge (j=ny): adjacent cells are (i-1,ny-1) and (i,ny-1)
    vertex_f_(i, ny) = lerp(xv, x0, x1, f_(i - 1, ny - 1), f_(i, ny - 1));
  }

}

void Solver::fill_vertex_from_bcs() {
  const std::size_t nx = m.nx();
  const std::size_t ny = m.ny();
  const double tt = t();

  auto set_dirichlet_vertex = [&](BoundaryID side, std::size_t i, std::size_t j) {
    double u = 0.0;
    if (!eq.dirichlet_vertex_value(side, i, j, tt, &u)) {
      throw std::runtime_error("Dirichlet BC: missing value.");
    }
    vertex_f_(i, j) = u;
  };

  // XMIN (i = 0)
  if (eq.bc_type(BoundaryID::XMIN) == BCType::Dirichlet) {
    for (std::size_t j = 0; j <= ny; ++j)
      set_dirichlet_vertex(BoundaryID::XMIN, 0, j);
  }

  // XMAX (i = nx)
  if (eq.bc_type(BoundaryID::XMAX) == BCType::Dirichlet) {
    for (std::size_t j = 0; j <= ny; ++j)
      set_dirichlet_vertex(BoundaryID::XMAX, nx, j);
  }

  // YMIN (j = 0)
  if (eq.bc_type(BoundaryID::YMIN) == BCType::Dirichlet) {
    for (std::size_t i = 0; i <= nx; ++i)
      set_dirichlet_vertex(BoundaryID::YMIN, i, 0);
  }

  // YMAX (j = ny)
  if (eq.bc_type(BoundaryID::YMAX) == BCType::Dirichlet) {
    for (std::size_t i = 0; i <= nx; ++i)
      set_dirichlet_vertex(BoundaryID::YMAX, i, ny);
  }
  
}

