#include <vector>
#include <memory>
#include <algorithm>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include "parameter.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace std;
using namespace Eigen ;
using namespace Rcpp;
using namespace hparameter;

// Custom comparator for stable two-column descending sort
struct ColumnComparator {
  const MatrixXd& mat;
  ColumnComparator(const MatrixXd& matrix) : mat(matrix) {}
  
  bool operator()(int a, int b) const {
    if(abs(mat(a, 0)-mat(b, 0))>EPSILON) 
      return mat(a, 0) > mat(b, 0); // Descending on col 0
    return mat(a, 1) > mat(b, 1);   // Descending on col 1 as tiebreak
  }
};

MatrixXd sort_matrix(const MatrixXd& mat) {
  const int n = mat.rows();
  std::vector<int> indices(n);
  std::iota(indices.begin(), indices.end(), 0);
  std::stable_sort(indices.begin(), indices.end(), ColumnComparator(mat));
  MatrixXd sorted(mat.rows(), mat.cols());
  for(int i = 0; i < n; ++i)
    sorted.row(i) = mat.row(indices[i]);
  return sorted;
}
  

double hsearch(const MatrixXd& xmat0, const MatrixXd& ymat0) {
  // FIX H5: removed `upper` and `visit` entirely.
  // In 2D the correct dominated y-vertex for x[i+1] can have a *larger*
  // lambda[0] than ymat[upper] (i.e. it sits earlier in the sorted ymat),
  // so the upper-pointer optimisation skips it and returns a wrong vertex.
  // Always scan from j=0.
  bool type1 = false;
  double max_y = 0, h0 = 0, dist = 0, h = 0;
  int tmp = 0, tmp_l = 0;
  int xrow = xmat0.rows(), yrow = ymat0.rows();
  
  MatrixXd xmat = sort_matrix(xmat0);
  MatrixXd ymat = sort_matrix(ymat0);
  
  // frontier4: indices of all type-4 vertices on V_tilde for current x-vertex
  std::vector<int> frontier4;
  frontier4.reserve(32);

  for(int i = 0; i < xrow; i++)
  {
    tmp = -1; max_y = SMALL_init; tmp_l = 0;
    frontier4.clear();

    // ----------------------------------------------------------------
    // First pass: collect V_tilde for x-vertex A.
    // Scan dominated vertices in descending lambda[0] order, tracking
    // the running maximum lambda[1].  A vertex B is on the Pareto
    // frontier iff its lambda[1] exceeds max_y at the moment it is seen.
    //
    // Stop early on type-1 or type-2 (Proposition A.2: at most one such
    // vertex on V_tilde; Lemma A.3 applies directly -- Case 2).
    // Collect ALL type-4 frontier vertices for Case 3.
    // ----------------------------------------------------------------
    for(int j = 0; j < yrow; j++)
    {
      if( (xmat(i,0)>ymat(j,0)) && (xmat(i,1)>ymat(j,1)) )
      {
        if(ymat(j,1)>max_y)
        {
          max_y = ymat(j,1);
          tmp = j;
          if(int(ymat(j,5)+0.5)!=4)
          {
            // Type-1 or type-2 on V_tilde: Case 2 applies directly.
            type1 = true;
            tmp_l = 1;
            break;
          }
          // Type-4 on V_tilde: collect for Case 3.
          frontier4.push_back(j);
          tmp_l++;
        }
      }
    }

    if(!type1){
      // No dominated vertex at all -> Case 1.
      if(tmp_l == 0) {
        tmp = -1;
      } else {
        // ----------------------------------------------------------------
        // Second pass (Case 3): for EACH type-4 vertex on V_tilde,
        // find its two Pareto predecessors along each axis direction,
        // apply Lemma A.3 to each predecessor found, and take the
        // minimum h0 over all frontier vertices and both directions.
        // (Eq. 60 in the paper: argmin rho1 over V_tilde.)
        //
        // Predecessors need NOT be dominated by A (FIX H6): a predecessor
        // B_x has the same lambda[1] as B4 and smaller lambda[0], which
        // may still exceed lambda_A[1] or lambda_A[0].
        // ----------------------------------------------------------------
        h0 = xmat(i,4);  // Case 1 default; overwritten if any predecessor found
        bool any_pred = false;

        for(int fi = 0; fi < (int)frontier4.size(); fi++)
        {
          int k = frontier4[fi];
          double k_lam0 = ymat(k,0);
          double k_lam1 = ymat(k,1);
          bool found_x = false, found_y = false;
          int  px = -1, py = -1;
          double best_dist_x = LARGE_init, best_dist_y = LARGE_init;

          for(int j = 0; j < yrow; j++)
          {
            // Predecessor along lambda[0] axis: same lambda[1], smaller lambda[0]
            if((abs(ymat(j,1)-k_lam1)<EPSILON) && (ymat(j,0)<k_lam0-EPSILON))
            {
              double d = max(abs(xmat(i,0)-ymat(j,0)), abs(xmat(i,1)-ymat(j,1)));
              if(d < best_dist_x) { best_dist_x = d; px = j; found_x = true; }
            }
            // Predecessor along lambda[1] axis: same lambda[0], smaller lambda[1]
            if((abs(ymat(j,0)-k_lam0)<EPSILON) && (ymat(j,1)<k_lam1-EPSILON))
            {
              double d = max(abs(xmat(i,0)-ymat(j,0)), abs(xmat(i,1)-ymat(j,1)));
              if(d < best_dist_y) { best_dist_y = d; py = j; found_y = true; }
            }
          }

          // Apply Lemma A.3 to the closer predecessor (argmin rho1)
          if(found_x || found_y)
          {
            int best = (best_dist_x <= best_dist_y) ? px : py;
            double candidate = ((xmat(i,0)+ymat(best,1))<(ymat(best,0)+xmat(i,1))) ?
                               abs(xmat(i,2)-ymat(best,2)) :
                               abs(xmat(i,3)-ymat(best,3));
            if(!any_pred || candidate < h0) { h0 = candidate; }
            any_pred = true;
          }
        }

        // If no predecessor found for any frontier vertex -> Case 1
        if(!any_pred) tmp = -1;
        else          tmp = 0;  // sentinel: h0 already set, skip distance block below
      }
    }
    // ----------------------------------------------------------------
    // Compute h0 for this x-vertex
    // ----------------------------------------------------------------
    if( tmp == -1 )
    {
      // Case 1: no dominated y-vertex, or Case 3 with no predecessor found
      h0 = xmat(i,4);
    } else if(type1) {
      // Case 2: type-1 or type-2 on V_tilde -- direct Lemma A.3
      if(int(ymat(tmp,5)+0.5)==1)
      {
        h0 = abs(xmat(i,4)-ymat(tmp,4));
      } else {
        h0 = ((xmat(i,0)+ymat(tmp,1))<(ymat(tmp,0)+xmat(i,1))) ?
             abs(xmat(i,2)-ymat(tmp,2)) :
             abs(xmat(i,3)-ymat(tmp,3));
      }
    }
    // Case 3: h0 already computed in the frontier4 loop above

    h = max(h, h0);

    type1 = false;
  }
  return h;
}
