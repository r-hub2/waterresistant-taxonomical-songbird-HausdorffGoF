#ifndef HAUSDORFF1D_H
#define HAUSDORFF1D_H

#include <vector>
#include <RcppEigen.h>

/** \file Hausdorff1d.h
 *  \brief Two-sample one-dimensional Hausdorff statistic between empirical CDFs,
 *         computed via a transformation-based method (tr) and a projection-based method (p).
 */



/// \brief Compute empirical CDF values at a set of evaluation points
/// \param data    sorted sample data
/// \param points  sorted evaluation points
/// \return ECDF value at each point in \p points
std::vector<double> compute_ecdf1d(const std::vector<double> &data,
                                   const std::vector<double> &points);

/// \brief Two-sample 1-D Hausdorff statistic - transformation method (core)
/// \param a  first sample as a 1-D Eigen array
/// \param b  second sample as a 1-D Eigen array
/// \return Hausdorff statistic between the two empirical CDFs
double H_stat_2s_1d_tr_cpp(const Eigen::ArrayXd &a, const Eigen::ArrayXd &b);

/// \brief Two-sample 1-D Hausdorff statistic - projection method (core)
/// \param a  first sample as a 1-D Eigen array
/// \param b  second sample as a 1-D Eigen array
/// \return Hausdorff statistic between the two empirical CDFs
double H_stat_2s_1d_p_cpp(const Eigen::ArrayXd &a, const Eigen::ArrayXd &b);


#endif /* HAUSDORFF1D_H */
