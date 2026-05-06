#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <map>
#include <cmath>
#include "Hausdorff1d.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

// ---------------------------------------------------------------------------
// Helper
// ---------------------------------------------------------------------------

std::vector<double> compute_ecdf1d(const std::vector<double> &data,
                                   const std::vector<double> &points)
{
    // Both vectors must already be sorted by the caller.
    std::vector<double> ecdf_values;
    int n = static_cast<int>(data.size());
    int start = 0, count = 0;
    for (double p : points) {
        for (int i = start; i < n; ++i) {
            if (data[i] <= p) { count++; start++; }
        }
        ecdf_values.push_back(static_cast<double>(count) / n);
    }
    return ecdf_values;
}

// ---------------------------------------------------------------------------
// Core - transformation method
// ---------------------------------------------------------------------------

double H_stat_2s_1d_tr_cpp(const ArrayXd &a, const ArrayXd &b)
{
    // Step 1: frequency counts for unique elements
    std::map<double, int> a_counts, b_counts;
    for (int i = 0; i < a.size(); ++i) a_counts[a(i)]++;
    for (int i = 0; i < b.size(); ++i) b_counts[b(i)]++;

    int m_a = static_cast<int>(a.size());
    int m_b = static_cast<int>(b.size());

    std::vector<double> a_unique, b_unique;
    std::vector<double> ay1, by1;

    for (auto const &pair : a_counts) {
        a_unique.push_back(pair.first);
        ay1.push_back(pair.second / static_cast<double>(m_a));
    }
    for (auto const &pair : b_counts) {
        b_unique.push_back(pair.first);
        by1.push_back(pair.second / static_cast<double>(m_b));
    }

    std::sort(a_unique.begin(), a_unique.end());
    std::sort(b_unique.begin(), b_unique.end());

    // Swap so that a_unique is the shorter sequence
    bool swapped = false;
    if (a_unique.size() > b_unique.size()) {
        swapped = true;
        std::swap(a_unique, b_unique);
        std::swap(ay1, by1);
        std::swap(m_a, m_b);
    }

    // Shift so that the minimum value is 0
    double min_val = std::min(a_unique.front(), b_unique.front());
    for (auto &val : a_unique) val -= min_val;
    for (auto &val : b_unique) val -= min_val;

    double max_val = std::max(a_unique.back(), b_unique.back());

    // Build interleaved aa / bb vectors
    int a_len = static_cast<int>(a_unique.size());
    int b_len = static_cast<int>(b_unique.size());
    std::vector<double> aa(2 * a_len, 0.0), bb(2 * b_len, 0.0);

    for (int i = 0; i < a_len; ++i) {
        double diff = (i == 0) ? a_unique[i] : (a_unique[i] - a_unique[i - 1]);
        aa[2 * i]     = -diff;
        aa[2 * i + 1] =  ay1[i];
    }
    for (int i = 0; i < b_len; ++i) {
        double diff = (i == 0) ? b_unique[i] : (b_unique[i] - b_unique[i - 1]);
        bb[2 * i]     = -diff;
        bb[2 * i + 1] =  by1[i];
    }

    int m = 2 * a_len;
    if (!a_unique.empty() && a_unique.back() == max_val) {
        aa.resize(2 * a_len - 1);
        m = static_cast<int>(aa.size());
        if (static_cast<int>(bb.size()) >= 2 * b_len && bb[2 * b_len - 1] != -max_val) {
            bb.push_back(max_val);
        }
    }

    // Cumulative sums: (ax, ay) for aa and (bx, by) for bb
    std::vector<double> ax, ay, bx, by;
    {
        double sx = 0.0, sy = 0.0;
        for (double v : aa) { sx += std::abs(v); sy += v; ax.push_back(sx); ay.push_back(sy); }
    }
    {
        double sx = 0.0, sy = 0.0;
        for (double v : bb) { sx += std::abs(v); sy += v; bx.push_back(sx); by.push_back(sy); }
    }

    int j = 0;
    int begining = std::min(2, m);
    if (!ay.empty() && ay[0] != 0.0) begining = 1;
    double target_ax = (begining > 0) ? ax[begining - 1] : 0.0;
    while (j < static_cast<int>(bx.size()) && bx[j] <= target_ax) j++;

    double max_H = 0.0;
    for (int i = begining; i <= m; ++i) {
        if (i == 0) continue;
        int idx = i - 1;
        if (idx >= static_cast<int>(ay.size())) break;

        int prev_j = std::max(j - 1, 0);
        if (prev_j >= static_cast<int>(by.size()))
            prev_j = static_cast<int>(by.size()) - 1;

        int parity = (j % 2 == 0) ? 1 : -1;
        double H_temp = ay[idx] - by[prev_j] + parity * (ax[idx] - bx[prev_j]);

        if ((H_temp < 0) == (i % 2 != 0)) {
            // record candidate (matches original logic)
        }

        if (i < m) {
            double next_ax = ax[i];
            while (j < static_cast<int>(bx.size()) && bx[j] <= next_ax) j++;
        }
        max_H = std::max(std::abs(H_temp), max_H);
    }

    return max_H / 2.0;
}

// ---------------------------------------------------------------------------
// Core - projection method
// ---------------------------------------------------------------------------

double H_stat_2s_1d_p_cpp(const ArrayXd &a, const ArrayXd &b)
{
    std::vector<double> a_vec(a.data(), a.data() + a.size());
    std::vector<double> b_vec(b.data(), b.data() + b.size());

    std::sort(a_vec.begin(), a_vec.end());
    std::sort(b_vec.begin(), b_vec.end());

    // Joint sorted unique values
    std::vector<double> xjoint(a_vec);
    xjoint.insert(xjoint.end(), b_vec.begin(), b_vec.end());
    std::sort(xjoint.begin(), xjoint.end());
    xjoint.erase(std::unique(xjoint.begin(), xjoint.end()), xjoint.end());

    // ECDF values at joint points
    std::vector<double> ay1 = compute_ecdf1d(a_vec, xjoint);
    std::vector<double> ay2 = compute_ecdf1d(b_vec, xjoint);

    // Remove duplicate values from the individual sample vectors
    a_vec.erase(std::unique(a_vec.begin(), a_vec.end(),
                            [](double x, double y){ return std::abs(x - y) < 1e-15; }),
                a_vec.end());
    b_vec.erase(std::unique(b_vec.begin(), b_vec.end(),
                            [](double x, double y){ return std::abs(x - y) < 1e-15; }),
                b_vec.end());

    // Build upper envelope (x1/y1) and lower envelope (x2/y2)
    std::vector<double> x1, y1, x2, y2;
    double tempy1 = 0.0, tempy2 = 0.0;
    for (int i = 0; i < static_cast<int>(xjoint.size()); ++i) {
        double cur_max = std::max(ay1[i], ay2[i]);
        if (cur_max > tempy1) {
            tempy1 = cur_max;
            x1.push_back(xjoint[i]);
            y1.push_back(tempy1);
        }
        double cur_min = std::min(ay1[i], ay2[i]);
        if (cur_min > tempy2) {
            x2.push_back(xjoint[i]);
            y2.push_back(tempy2);
            tempy2 = cur_min;
        }
    }

    // Extend x1/y1 and x2/y2 with duplicated entries
    int m1 = static_cast<int>(x1.size());
    int m2 = static_cast<int>(x2.size());
    std::vector<double> new_x1, new_y1, new_x2, new_y2;

    for (int i = 0; i < 2 * m1; ++i) { int idx = i / 2; new_x1.push_back(x1[idx]); new_y1.push_back(y1[idx]); }
    new_x1.push_back(xjoint.back());
    new_y1.insert(new_y1.begin(), 0.0);

    new_x2.push_back(xjoint[0]);
    for (int i = 0; i < 2 * m2; ++i) { int idx = i / 2; new_x2.push_back(x2[idx]); new_y2.push_back(y2[idx]); }
    new_y2.push_back(1.0);

    m1 = static_cast<int>(new_x1.size());
    m2 = static_cast<int>(new_x2.size());

    // Ensure new_x1/y1 is the shorter sequence
    if (m1 > m2) {
        std::swap(new_x1, new_x2);
        std::swap(new_y1, new_y2);
        std::reverse(new_x1.begin(), new_x1.end());
        std::reverse(new_x2.begin(), new_x2.end());
        std::reverse(new_y1.begin(), new_y1.end());
        std::reverse(new_y2.begin(), new_y2.end());
        for (auto &x : new_x1) x = -x;
        for (auto &x : new_x2) x = -x;
        for (auto &y : new_y1) y = 1.0 - y;
        for (auto &y : new_y2) y = 1.0 - y;
        std::swap(m1, m2);
    }

    // Parameter vectors
    std::vector<double> par1, par2;
    for (int i = 1; i < m1; i += 2)
        par1.push_back(new_x1[i] + new_y1[i]);
    for (int i = 0; i < m2; ++i)
        par2.push_back(new_x2[i] + new_y2[i]);

    double h = 0.0;
    int loc = 0, i = 0, start = 0, end_idx = 0, batch = 0;
    bool newbatch = true, bool_change = false;

    for (int j = 0; (j < static_cast<int>(par2.size())) && (i < static_cast<int>(par1.size())); ) {
        if (par1[i] > par2[j]) {
            bool_change = (batch == j);
            j++;
            newbatch = true;
        } else {
            bool_change = false;
            if (newbatch) { start = i; end_idx = i; batch = j; newbatch = false; }
            else          { end_idx++; }
            i++;
        }
        if (bool_change) {
            h = std::max(h, (batch % 2 == 1
                             ? new_y1[2 * end_idx + 2] - new_y2[batch]
                             : new_x2[batch]           - new_x1[2 * start + 1]));
        }
    }
    h = std::max(h, (batch % 2 == 1
                     ? new_y1[2 * end_idx + 2] - new_y2[batch]
                     : new_x2[batch]            - new_x1[2 * start + 1]));

    return h;
}

// ---------------------------------------------------------------------------
// Rcpp wrappers
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
double H_stat_2s_1d_tr(NumericVector a, NumericVector b)
{
    ArrayXd ea = Rcpp::as<ArrayXd>(a);
    ArrayXd eb = Rcpp::as<ArrayXd>(b);
    return H_stat_2s_1d_tr_cpp(ea, eb);
}

// [[Rcpp::export]]
double H_stat_2s_1d_p(NumericVector a, NumericVector b)
{
    ArrayXd ea = Rcpp::as<ArrayXd>(a);
    ArrayXd eb = Rcpp::as<ArrayXd>(b);
    return H_stat_2s_1d_p_cpp(ea, eb);
}
