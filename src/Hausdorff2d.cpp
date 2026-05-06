// Hausdorff2d.cpp
//
// H_stat_2s_2d_cpp: Hausdorff distance H(F_m^c, G_n^c) under the Chebyshev
// metric, computed via a single hsearch() call (Appendix A, paper).
//
// The algorithm is one-directional by construction:
//   projection_x = V_loc(F_m | G_n)  [locally-farthest vertices of F_m]
//   projection_y = all vertices of G_n
//
// Symmetry H(x,y) = H(y,x) follows from correct V_loc, not from two calls:
//   - Concave vertices (where F_m > G_n) capture F-above-G contributions.
//   - Convex vertices  (where G_n > F_m) capture G-above-F contributions.
// Both are in V_loc of F_m, so a single pass suffices.
//
// The key fix relative to the original code: the convex vertex condition
// is simply  G_n(ax,ay) > F_m(ax-,ay-)  (i.e. zy > z4).
// The extra guard |z4-z3| > eps that appeared in the R code is REMOVED:
// for generic continuous data z3 = z4 at every interior omnidirectional jump
// (because z3-z4 counts observations with x1 < ax AND x2 = ay exactly, which
// is zero), so the guard silently suppresses ALL interior convex vertices.
//
// Column layouts expected by hsearch() (0-indexed):
//   projection_x : [lambda0, lambda1, loc1, loc2, z]        (5 cols)
//   projection_y : [lambda0, lambda1, loc1, loc2, z, type]  (6 cols)

#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <memory>
#include "Hausdorff2d.h"
#include "Hausdorffsearch.h"
#include "fastCDF.h"
#include "parameter.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;
using namespace hparameter;

static MatrixXd rows_to_matrix(const std::vector<RowVectorXd> &rows, int ncols)
{
    if (rows.empty()) return MatrixXd(0, ncols);
    MatrixXd m(rows.size(), ncols);
    for (int i = 0; i < (int)rows.size(); ++i) m.row(i) = rows[i];
    return m;
}

double H_stat_2s_2d_cpp(const MatrixXd &x, const MatrixXd &y, double tol)
{
    const int mx = x.rows(), my = y.rows();

    // Sorted marginals
    VectorXd x1 = x.col(0), x2 = x.col(1);
    std::sort(x1.data(), x1.data() + mx);
    std::sort(x2.data(), x2.data() + mx);
    VectorXd y1 = y.col(0), y2 = y.col(1);
    std::sort(y1.data(), y1.data() + my);
    std::sort(y2.data(), y2.data() + my);

    const double M = std::round(std::max(x.maxCoeff(), y.maxCoeff())) + 2.5;

    // fastCDF expects (dimension x nbSim)
    ArrayXXd x_t(2, mx), y_t(2, my);
    for (int j = 0; j < mx; ++j) { x_t(0,j)=x(j,0); x_t(1,j)=x(j,1); }
    for (int j = 0; j < my; ++j) { y_t(0,j)=y(j,0); y_t(1,j)=y(j,1); }
    ArrayXd ones_x = ArrayXd::Ones(mx);
    ArrayXd ones_y = ArrayXd::Ones(my);

    // ------------------------------------------------------------------
    // Grid 1: F_m on (x1-tol, x1, M) x (x2-tol, x2, M)  size ng = 2*mx+1
    //   z_x[i0 + i1*ng] = F_m(g1[i0], g2[i1])
    //   Accessors (0-indexed xi,yi in 0..mx-1):
    //     F_m(x1[xi],   x2[yi]  ) = z_x_at(2*xi+1, 2*yi+1)   z1
    //     F_m(x1[xi],   x2[yi]-t) = z_x_at(2*xi+1, 2*yi  )   z2
    //     F_m(x1[xi]-t, x2[yi]  ) = z_x_at(2*xi,   2*yi+1)   z3
    //     F_m(x1[xi]-t, x2[yi]-t) = z_x_at(2*xi,   2*yi  )   z4
    //   Boundary:
    //     F_m(x1[xi], M)    = z_x_at(2*xi+1, 2*mx)
    //     F_m(x1[xi]-t, M)  = z_x_at(2*xi,   2*mx)
    //     F_m(M, x2[yi])    = z_x_at(2*mx,   2*yi+1)
    //     F_m(M, x2[yi]-t)  = z_x_at(2*mx,   2*yi  )
    // ------------------------------------------------------------------
    const int ng = 2*mx + 1;
    VectorXd g1(ng), g2(ng);
    for (int i=0;i<mx;++i){g1(2*i)=x1(i)-tol; g1(2*i+1)=x1(i);}
    for (int i=0;i<mx;++i){g2(2*i)=x2(i)-tol; g2(2*i+1)=x2(i);}
    g1(2*mx)=M; g2(2*mx)=M;
    auto pg1 = std::make_shared<ArrayXd>(g1.array());
    auto pg2 = std::make_shared<ArrayXd>(g2.array());
    ArrayXd z_x = fastCDF(x_t, {pg1, pg2}, ones_x);
    auto z_x_at = [&](int i0,int i1){ return z_x(i0 + i1*ng); };

    // ------------------------------------------------------------------
    // Grid 2: G_n on (x1, M) x (x2, M)  size ngc = mx+1
    //   z_yx_at(xi, yi) = G_n(x1[xi], x2[yi])
    //   z_yx_at(xi, mx) = G_n(x1[xi], M)
    //   z_yx_at(mx, yi) = G_n(M, x2[yi])
    // ------------------------------------------------------------------
    const int ngc = mx+1;
    VectorXd gc1(ngc), gc2(ngc);
    for (int i=0;i<mx;++i){gc1(i)=x1(i); gc2(i)=x2(i);}
    gc1(mx)=M; gc2(mx)=M;
    auto pgc1 = std::make_shared<ArrayXd>(gc1.array());
    auto pgc2 = std::make_shared<ArrayXd>(gc2.array());
    ArrayXd z_yx = fastCDF(y_t, {pgc1, pgc2}, ones_y);
    auto z_yx_at = [&](int i0,int i1){ return z_yx(i0 + i1*ngc); };

    // ==================================================================
    // Step 4: Build projection_x = V_loc(F_m | G_n)
    //
    // At each omnidirectional jump location (ax, ay) of F_m:
    //   Concave vertex (ax, ay, z1): in V_loc iff z1 > G_n(ax,ay)
    //   Convex  vertex (ax, ay, z4): in V_loc iff G_n(ax,ay) > z4
    //
    // NOTE: NO |z4-z3| guard.  For generic continuous data, z3=z4 at every
    // interior jump (z3-z4 = #{x1<ax, x2=ay exactly}/m = 0), so adding the
    // guard would suppress ALL interior convex vertices -- breaking symmetry.
    // ==================================================================
    std::vector<RowVectorXd> px_rows;

    // Interior: (x1[xi], x2[yi]) for xi,yi = 0..mx-1
    for (int xi=0; xi<mx; ++xi) {
        const double ax = x1(xi);
        for (int yi=0; yi<mx; ++yi) {
            const double ay = x2(yi);
            const double z1 = z_x_at(2*xi+1, 2*yi+1);
            const double z2 = z_x_at(2*xi+1, 2*yi  );
            const double z3 = z_x_at(2*xi,   2*yi+1);
            const double z4 = z_x_at(2*xi,   2*yi  );
            if (z1==z2 || z1==z3) continue;  // not omnidirectional
            const double zy = z_yx_at(xi, yi);
            if (z1 > zy) {                    // concave vertex in V_loc
                RowVectorXd r(PROJECTION_X_ELEMENTS); r << ax+z1, ay+z1, ax, ay, z1;
                px_rows.push_back(r);
            }
            if (zy > z4 + 1e-14) {            // convex vertex in V_loc
                RowVectorXd r(PROJECTION_X_ELEMENTS); r << ax+z4, ay+z4, ax, ay, z4;
                px_rows.push_back(r);
            }
        }
        // Boundary at (x1[xi], M)
        {
            const double z1b = z_x_at(2*xi+1, 2*mx);
            const double z4b = z_x_at(2*xi,   2*mx);
            if (z1b != z4b) {
                const double zyb = z_yx_at(xi, mx);
                if (z1b > zyb) {
                    RowVectorXd r(PROJECTION_X_ELEMENTS); r << x1(xi)+z1b, M+z1b, x1(xi), M, z1b;
                    px_rows.push_back(r);
                }
                if (zyb > z4b + 1e-14) {
                    RowVectorXd r(PROJECTION_X_ELEMENTS); r << x1(xi)+z4b, M+z4b, x1(xi), M, z4b;
                    px_rows.push_back(r);
                }
            }
        }
    }
    // Boundary at (M, x2[yi])
    for (int yi=0; yi<mx; ++yi) {
        const double z1b = z_x_at(2*mx, 2*yi+1);
        const double z4b = z_x_at(2*mx, 2*yi  );
        if (z1b != z4b) {
            const double zyb = z_yx_at(mx, yi);
            if (z1b > zyb) {
                RowVectorXd r(PROJECTION_X_ELEMENTS); r << M+z1b, x2(yi)+z1b, M, x2(yi), z1b;
                px_rows.push_back(r);
            }
            if (zyb > z4b + 1e-14) {
                RowVectorXd r(PROJECTION_X_ELEMENTS); r << M+z4b, x2(yi)+z4b, M, x2(yi), z4b;
                px_rows.push_back(r);
            }
        }
    }
    MatrixXd projection_x = rows_to_matrix(px_rows, 5);

    // ==================================================================
    // Step 5: Build projection_y = all vertices of G_n (types 1,2,3,4)
    //   Grid: (y1-tol, y1) x (y2-tol, y2)  size ngy = 2*my
    // ==================================================================
    const int ngy = 2*my;
    VectorXd gy1(ngy), gy2(ngy);
    for (int i=0;i<my;++i){gy1(2*i)=y1(i)-tol; gy1(2*i+1)=y1(i);}
    for (int i=0;i<my;++i){gy2(2*i)=y2(i)-tol; gy2(2*i+1)=y2(i);}
    auto pgy1 = std::make_shared<ArrayXd>(gy1.array());
    auto pgy2 = std::make_shared<ArrayXd>(gy2.array());
    ArrayXd z_y = fastCDF(y_t, {pgy1, pgy2}, ones_y);
    auto z_y_at = [&](int i0,int i1){ return z_y(i0 + i1*ngy); };

    std::vector<RowVectorXd> py_rows;
    for (int xi=0; xi<my; ++xi) {
        const double b1 = y1(xi);
        for (int yi=0; yi<my; ++yi) {
            const double b2 = y2(yi);
            const double z1 = z_y_at(2*xi+1, 2*yi+1);
            const double z2 = z_y_at(2*xi+1, 2*yi  );
            const double z3 = z_y_at(2*xi,   2*yi+1);
            const double z4 = z_y_at(2*xi,   2*yi  );
            if (z1==z2 || z1==z3) continue;
            { RowVectorXd r(PROJECTION_Y_ELEMENTS); r<<b1+z2,b2+z2,b1,b2,z2,2.0; py_rows.push_back(r); }
            if (z4!=z2) {
                { RowVectorXd r(PROJECTION_Y_ELEMENTS); r<<b1+z4,b2+z4,b1,b2,z4,4.0; py_rows.push_back(r); }
                if (z2!=z3) {
                    { RowVectorXd r(PROJECTION_Y_ELEMENTS); r<<b1+z3,b2+z3,b1,b2,z3,3.0; py_rows.push_back(r); }
                }
            }
            { RowVectorXd r(PROJECTION_Y_ELEMENTS); r<<b1+z1,b2+z1,b1,b2,z1,1.0; py_rows.push_back(r); }
        }
    }
    MatrixXd projection_y = rows_to_matrix(py_rows, 6);

    return hsearch(projection_x, projection_y);
}

// [[Rcpp::export]]
double H_stat_2s_2d(NumericMatrix x, NumericMatrix y, double tol = 1e-6)
{
    return H_stat_2s_2d_cpp(as<MatrixXd>(x), as<MatrixXd>(y), tol);
}


