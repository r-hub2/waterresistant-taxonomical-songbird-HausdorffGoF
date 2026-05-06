#ifndef HAUSDORFF2D_H
#define HAUSDORFF2D_H

#include <RcppEigen.h>

/** \file Hausdorff2d.h
 *  \brief Two-sample two-dimensional Hausdorff statistic between bivariate
 *         empirical CDF surfaces under the Chebyshev metric.
 *
 *  Implements H(F^c_m, G^c_n) as described in Appendix A of
 *  Dimitrova, Jia & Kaishev (2025).  The projection matrices fed to
 *  hsearch() are constructed entirely from bivariate ECDF counting;
 *  no fastCDF grid-ordering dependency is required.
 *
 *  Column layouts (0-indexed) passed to hsearch():
 *    projection_x  [proj1, proj2, loc1, loc2, z]         (5 cols)
 *    projection_y  [proj1, proj2, loc1, loc2, z, type]   (6 cols)
 *  where proj_i = loc_i + z for both matrices.
 */



/// \brief Two-sample 2-D Hausdorff statistic - core Eigen implementation
/// \param x    first bivariate sample, size m x 2
/// \param y    second bivariate sample, size n x 2
/// \param tol  tolerance used to probe left/below limits of the ECDF
/// \return Hausdorff distance H(F^c_m, G^c_n) under the Chebyshev metric
double H_stat_2s_2d_cpp(const Eigen::MatrixXd &x,
                        const Eigen::MatrixXd &y,
                        double tol = 1e-6);


#endif /* HAUSDORFF2D_H */
