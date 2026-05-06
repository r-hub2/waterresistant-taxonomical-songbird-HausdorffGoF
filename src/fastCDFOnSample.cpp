// Copyright (C) 2020 EDF
// All Rights Reserved
// This code is published under the GNU Lesser General Public License (GNU LGPL)
#include <RcppEigen.h>
#include <Rcpp.h>
#include "nDDominanceAlone.h"
#include <vector>
#include <memory>

// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen ;
using namespace Rcpp;


Eigen::ArrayXd fastCDFOnSample(const Eigen::ArrayXXd &p_x, const Eigen::ArrayXd &p_y)
{

    ArrayXd cdfDivid(p_y.size());

    // dominance excluding current point
    nDDominanceAlone(p_x, p_y, cdfDivid);
    cdfDivid += p_y;
    return cdfDivid / p_y.size();
}

// Wrapper function for R
// [[Rcpp::export]]
NumericVector fastCDFOnSample_Rcpp(NumericMatrix p_x_r) {
  // Convert R inputs to Eigen types (no copy)
  // ------------------------------------------------------------
  Map<ArrayXXd> p_x(as<Map<ArrayXXd>>(p_x_r));
  
  ArrayXd p_y = ArrayXd::Ones(p_x.cols());
  
  // ------------------------------------------------------------
  ArrayXd result = fastCDFOnSample(p_x, p_y);
  
  // Convert Eigen result to Rcpp NumericVector (safe copy)
  // ------------------------------------------------------------
  return NumericVector(result.data(), result.data() + result.size());
}
