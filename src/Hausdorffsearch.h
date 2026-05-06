#ifndef HAUSDORFFSEARCH_H
#define HAUSDORFFSEARCH_H

#include <Eigen/Dense>

Eigen::MatrixXd sort_matrix(const Eigen::MatrixXd& mat);
double hsearch(const Eigen::MatrixXd& xmat0, const Eigen::MatrixXd& ymat0);

  
#endif /* HAUSDORFFSEARCH_H */
