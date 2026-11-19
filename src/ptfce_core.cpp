#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <unordered_map>

using namespace Rcpp;

namespace {

const double LOG_2PI = 1.8378770664093454836;       // log(2 * pi)
const double HALF_LOG_2PI = 0.91893853320467274178; // 0.5 * log(2 * pi)
const double LGAMMA_2_5 = 0.28468287047291918075;   // lgamma(2.5)

// Compute log(Phi(-h)) + 0.5 * h^2 accurately without overflow or cancellation
// for all h >= 0 (including h >> 38 where Phi(-h) underflows in IEEE double).
inline double log_mills_norm(double h) {
  if (h < 12.0) {
    return R::pnorm5(h, 0.0, 1.0, 0, 1) + 0.5 * h * h;
  }
  // Laplace continued fraction for Mill's ratio:
  // exp(h^2 / 2) * int_h^Inf exp(-t^2 / 2) dt = 1 / (h + 1/(h + 2/(h + ...)))
  double cf = 0.0;
  for (int k = 32; k >= 1; --k) {
    cf = static_cast<double>(k) / (h + cf);
  }
  return -HALF_LOG_2PI - std::log(h + cf);
}

// Log of expected cluster size Es(h, V, Rd) from GRF theory
inline double log_es_scalar(double h, double logV, double logRd) {
  if (h >= 1.1) {
    double h2 = h * h;
    // log(V) + log(Phi(-h)) - (log(Rd) + log(h^2 - 1) - 0.5*h^2 - 2*log(2*pi))
    // = log(V) - log(Rd) + 2*log(2*pi) - log(h^2 - 1) + (log(Phi(-h)) + 0.5*h^2)
    return logV - logRd + 2.0 * LOG_2PI - std::log(h2 - 1.0) + log_mills_norm(h);
  }
  return logV + R::pnorm5(h, 0.0, 1.0, 0, 1);
}

// Log of unnormalized posterior density:
// l(h) = log(phi(h)) + log(lambda(h)) - c^(2/3) * lambda(h)
// Strictly concave on (1.1, Inf).
inline double log_post_unnorm(double h, double logV, double logRd, double c23) {
  double log_es = log_es_scalar(h, logV, logRd);
  double log_lam = (-2.0 / 3.0) * (log_es - LGAMMA_2_5);
  double lam_c23 = std::exp(log_lam) * c23;
  double log_phi = -HALF_LOG_2PI - 0.5 * h * h;
  return log_phi + log_lam - lam_c23;
}

// 15-point Gauss-Kronrod rule on [a, b] for exp(log_post_unnorm(x) - shift)
struct GKResult {
  double val;
  double err;
};

inline GKResult gk15_panel(double a, double b, double shift,
                           double logV, double logRd, double c23) {
  // Positive abscissae and weights for G7 / K15 on [-1, 1]
  static const double xgk[8] = {
    0.991455371120812639206854697526329,
    0.949107912342758524526189684047851,
    0.864864423359769072789712788640926,
    0.741531185599394439863864773280788,
    0.586087235467691130294144838258730,
    0.405845151377397166906606412076961,
    0.207784955007898467600689403773245,
    0.000000000000000000000000000000000
  };
  static const double wgk[8] = {
    0.022935322010529224963732008058970,
    0.063092092629978553290700663189204,
    0.104790010322250183839876322541518,
    0.140653259715525918745189590510238,
    0.169004726639267902826583426598550,
    0.190350578064785409913256402421014,
    0.204432940075298892414161999234649,
    0.209482141084727828012999174891714
  };
  static const double wg[4] = {
    0.129484966168869693270611432679082,
    0.279705391489276667901467771423780,
    0.381830050505118944950369775488975,
    0.417959183673469387755102040816327
  };

  double center = 0.5 * (a + b);
  double half = 0.5 * (b - a);
  double f_center = std::exp(log_post_unnorm(center, logV, logRd, c23) - shift);
  double res_g = wg[3] * f_center;
  double res_k = wgk[7] * f_center;

  for (int j = 0; j < 7; ++j) {
    double absc = half * xgk[j];
    double fval1 = std::exp(log_post_unnorm(center - absc, logV, logRd, c23) - shift);
    double fval2 = std::exp(log_post_unnorm(center + absc, logV, logRd, c23) - shift);
    double fsum = fval1 + fval2;
    res_k += wgk[j] * fsum;
    if (j % 2 == 1) {
      res_g += wg[j / 2] * fsum;
    }
  }
  res_k *= half;
  res_g *= half;
  return {res_k, std::fabs(res_k - res_g)};
}

double adaptive_gk15(double a, double b, double shift,
                     double logV, double logRd, double c23,
                     double tol, int depth = 0) {
  if (b <= a) return 0.0;
  GKResult r = gk15_panel(a, b, shift, logV, logRd, c23);
  if (r.err <= tol * std::max(1e-16, r.val) || depth >= 18) {
    return r.val;
  }
  double mid = 0.5 * (a + b);
  return adaptive_gk15(a, mid, shift, logV, logRd, c23, tol, depth + 1) +
         adaptive_gk15(mid, b, shift, logV, logRd, c23, tol, depth + 1);
}

// Find unique global mode of log_post_unnorm on [a, b_max] (b_max <= Inf)
inline double find_mode_on(double a, double b_max, double logV, double logRd, double c23, double &max_val) {
  double fa = log_post_unnorm(a, logV, logRd, c23);
  double eps = std::min(1e-5, 0.25 * (b_max - a));
  if (eps <= 0.0) {
    max_val = fa;
    return a;
  }
  double fa_eps = log_post_unnorm(a + eps, logV, logRd, c23);
  if (fa_eps <= fa) {
    // By strict concavity on (1.1, Inf), mode on [a, b_max] is at x = a
    max_val = fa;
    return a;
  }
  // Bracket peak
  double lo = a;
  double hi = std::min(b_max, a + 0.5);
  while (hi < b_max && log_post_unnorm(hi, logV, logRd, c23) >
                       log_post_unnorm(std::max(a, hi - 1e-4), logV, logRd, c23)) {
    lo = hi;
    double step = std::max(0.5, hi - a);
    hi = std::min(b_max, hi + step);
    if (hi == b_max) break;
  }
  lo = a; // safe bracket [a, hi]
  // Golden-section search on [lo, hi]
  const double invphi = 0.6180339887498948482;
  const double invphi2 = 0.3819660112501051518;
  double x1 = lo + invphi2 * (hi - lo);
  double x2 = lo + invphi * (hi - lo);
  double f1 = log_post_unnorm(x1, logV, logRd, c23);
  double f2 = log_post_unnorm(x2, logV, logRd, c23);
  for (int it = 0; it < 48; ++it) {
    if (f1 < f2) {
      lo = x1;
      x1 = x2;
      f1 = f2;
      x2 = lo + invphi * (hi - lo);
      f2 = log_post_unnorm(x2, logV, logRd, c23);
    } else {
      hi = x2;
      x2 = x1;
      f2 = f1;
      x1 = lo + invphi2 * (hi - lo);
      f1 = log_post_unnorm(x1, logV, logRd, c23);
    }
  }
  double mode = 0.5 * (lo + hi);
  max_val = log_post_unnorm(mode, logV, logRd, c23);
  double fb = (std::isfinite(b_max)) ? log_post_unnorm(b_max, logV, logRd, c23) : -INFINITY;
  if (fb > max_val) {
    max_val = fb;
    return b_max;
  }
  return mode;
}

// Compute log( int_a^Inf exp(log_post_unnorm(x)) dx )
double log_integral_tail(double a, double logV, double logRd, double c23) {
  double max_f = 0.0;
  double mode = find_mode_on(a, INFINITY, logV, logRd, c23, max_f);

  // Bracket upper cutoff b > mode where log_post_unnorm(b) <= max_f - 50
  const double drop = 50.0; // exp(-50) ~ 1.9e-22
  double step = 0.25;
  while (log_post_unnorm(mode + step, logV, logRd, c23) > max_f - drop) {
    step *= 2.0;
  }
  double lo = mode + ((step > 0.25) ? 0.5 * step : 0.0);
  double hi = mode + step;
  for (int it = 0; it < 28; ++it) {
    double mid = 0.5 * (lo + hi);
    if (log_post_unnorm(mid, logV, logRd, c23) > max_f - drop) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  double b = hi;

  double int_left = (mode > a) ? adaptive_gk15(a, mode, max_f, logV, logRd, c23, 1e-12) : 0.0;
  double int_right = adaptive_gk15(mode, b, max_f, logV, logRd, c23, 1e-12);
  return max_f + std::log(int_left + int_right);
}

// Compute log( int_a^b exp(log_post_unnorm(x)) dx ) for finite [a, b]
double log_integral_interval(double a, double b, double logV, double logRd, double c23) {
  if (b <= a) return -INFINITY;
  double max_f = 0.0;
  double mode = find_mode_on(a, b, logV, logRd, c23, max_f);
  double int_l = (mode > a) ? adaptive_gk15(a, mode, max_f, logV, logRd, c23, 1e-12) : 0.0;
  double int_r = (b > mode) ? adaptive_gk15(mode, b, max_f, logV, logRd, c23, 1e-12) : 0.0;
  return max_f + std::log(int_l + int_r);
}

// Exact -log(P(H > actH | C = c)) using optional cached log_den
inline double neglog_pvox_clust_cached(double actH, double ZestThr,
                                       double logV, double logRd, double c23,
                                       double log_den) {
  if (actH <= ZestThr) {
    return -R::pnorm5(actH, 0.0, 1.0, 0, 1);
  }
  double log_num = log_integral_tail(actH, logV, logRd, c23);
  double diff = log_den - log_num;
  if (diff < 0.05) {
    // Near ZestThr where P ~ 1: integrate [ZestThr, actH] directly so
    // -log(P) = log(1 + I_left / I_num) > 0 with zero cancellation!
    double log_left = log_integral_interval(ZestThr, actH, logV, logRd, c23);
    double ratio_log = log_left - log_num;
    if (ratio_log > 35.0) {
      return ratio_log;
    }
    return std::log1p(std::exp(ratio_log));
  }
  return diff;
}

// 4-point cubic Lagrange interpolation on an equispaced grid y[0 .. M-1] (M >= 4)
inline double interp_cubic_equispaced(const std::vector<double> &y,
                                      double x_min, double inv_dx, double x) {
  int M = static_cast<int>(y.size());
  double u = (x - x_min) * inv_dx;
  if (u <= 0.0) return y[0];
  if (u >= static_cast<double>(M - 1)) return y[M - 1];

  int i1 = static_cast<int>(std::floor(u));
  if (i1 < 1) i1 = 1;
  if (i1 > M - 3) i1 = M - 3;

  double t = u - static_cast<double>(i1);
  double w0 = -t * (t - 1.0) * (t - 2.0) * (1.0 / 6.0);
  double w1 =  (t + 1.0) * (t - 1.0) * (t - 2.0) * 0.5;
  double w2 = -(t + 1.0) * t * (t - 2.0) * 0.5;
  double w3 =  (t + 1.0) * t * (t - 1.0) * (1.0 / 6.0);

  return w0 * y[i1 - 1] + w1 * y[i1] + w2 * y[i1 + 1] + w3 * y[i1 + 2];
}

// Separable 3D Gaussian convolution (zero-padded boundaries; used for normalized mask convolution)
void smooth_3d_separable(std::vector<double> &vol, int nx, int ny, int nz, double fwhm_vox) {
  if (fwhm_vox <= 0.0) return;
  double sigma = fwhm_vox / std::sqrt(8.0 * std::log(2.0));
  if (sigma < 1e-4) return;
  int rad = std::max(1, static_cast<int>(std::ceil(3.5 * sigma)));
  std::vector<double> kern(2 * rad + 1);
  double ksum = 0.0;
  for (int k = -rad; k <= rad; ++k) {
    double w = std::exp(-0.5 * (k * k) / (sigma * sigma));
    kern[k + rad] = w;
    ksum += w;
  }
  for (double &w : kern) w /= ksum;

  int nxy = nx * ny;
  std::vector<double> tmp(vol.size(), 0.0);

  // Pass X
  if (nx > 1) {
    for (int z = 0; z < nz; ++z) {
      for (int y = 0; y < ny; ++y) {
        int base = y * nx + z * nxy;
        for (int x = 0; x < nx; ++x) {
          double acc = 0.0;
          int k_lo = std::max(-rad, -x);
          int k_hi = std::min(rad, nx - 1 - x);
          for (int k = k_lo; k <= k_hi; ++k) {
            acc += kern[k + rad] * vol[base + x + k];
          }
          tmp[base + x] = acc;
        }
      }
    }
    vol.swap(tmp);
  }

  // Pass Y
  if (ny > 1) {
    std::fill(tmp.begin(), tmp.end(), 0.0);
    for (int z = 0; z < nz; ++z) {
      int z_off = z * nxy;
      for (int y = 0; y < ny; ++y) {
        int k_lo = std::max(-rad, -y);
        int k_hi = std::min(rad, ny - 1 - y);
        for (int x = 0; x < nx; ++x) {
          double acc = 0.0;
          for (int k = k_lo; k <= k_hi; ++k) {
            acc += kern[k + rad] * vol[z_off + (y + k) * nx + x];
          }
          tmp[z_off + y * nx + x] = acc;
        }
      }
    }
    vol.swap(tmp);
  }

  // Pass Z
  if (nz > 1) {
    std::fill(tmp.begin(), tmp.end(), 0.0);
    for (int z = 0; z < nz; ++z) {
      int k_lo = std::max(-rad, -z);
      int k_hi = std::min(rad, nz - 1 - z);
      for (int y = 0; y < ny; ++y) {
        int y_off = y * nx;
        for (int x = 0; x < nx; ++x) {
          double acc = 0.0;
          for (int k = k_lo; k <= k_hi; ++k) {
            acc += kern[k + rad] * vol[(z + k) * nxy + y_off + x];
          }
          tmp[z * nxy + y_off + x] = acc;
        }
      }
    }
    vol.swap(tmp);
  }
}

} // namespace

// [[Rcpp::export]]
NumericVector es_cpp(NumericVector h, double V, double Rd) {
  int n = h.size();
  NumericVector out(n);
  double logV = std::log(V);
  double logRd = std::log(Rd);
  for (int i = 0; i < n; ++i) {
    if (NumericVector::is_na(h[i])) {
      out[i] = NA_REAL;
    } else {
      out[i] = std::exp(log_es_scalar(h[i], logV, logRd));
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericVector neglog_pvox_clust_cpp(double V, double Rd, NumericVector c,
                                    NumericVector actH, double ZestThr = 1.3) {
  int nc = c.size();
  int nh = actH.size();
  int n = std::max(nc, nh);
  if (nc == 0 || nh == 0) return NumericVector(0);

  NumericVector out(n);
  double logV = std::log(V);
  double logRd = std::log(Rd);
  std::unordered_map<double, double> den_cache;

  for (int i = 0; i < n; ++i) {
    double ci = c[i % nc];
    double hi = actH[i % nh];
    if (NumericVector::is_na(hi)) {
      out[i] = NA_REAL;
      continue;
    }
    if (NumericVector::is_na(ci) || ci <= 0.0 || hi <= ZestThr) {
      out[i] = -R::pnorm5(hi, 0.0, 1.0, 0, 1);
      continue;
    }
    double c23 = std::pow(ci, 2.0 / 3.0);
    double lden;
    auto it = den_cache.find(ci);
    if (it == den_cache.end()) {
      lden = log_integral_tail(ZestThr, logV, logRd, c23);
      den_cache[ci] = lden;
    } else {
      lden = it->second;
    }
    out[i] = neglog_pvox_clust_cached(hi, ZestThr, logV, logRd, c23, lden);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector dvox_clust_cpp(NumericVector h, double V, double Rd, double c,
                             double ZestThr = 1.3) {
  int n = h.size();
  NumericVector out(n);
  if (NumericVector::is_na(c) || c <= 0.0) {
    for (int i = 0; i < n; ++i) {
      out[i] = R::dnorm4(h[i], 0.0, 1.0, 0);
    }
    return out;
  }
  double logV = std::log(V);
  double logRd = std::log(Rd);
  double c23 = std::pow(c, 2.0 / 3.0);
  double log_den = log_integral_tail(ZestThr, logV, logRd, c23);
  for (int i = 0; i < n; ++i) {
    if (NumericVector::is_na(h[i])) {
      out[i] = NA_REAL;
    } else if (h[i] < ZestThr) {
      out[i] = 0.0;
    } else {
      out[i] = std::exp(log_post_unnorm(h[i], logV, logRd, c23) - log_den);
    }
  }
  return out;
}

// Single-tail pTFCE enhancement using single-pass descending Disjoint-Set Union (Union-Find)
// Supports both stationary integer cluster sizes and non-stationary continuous resel weights.
// [[Rcpp::export]]
NumericVector ptfce_tail_cpp(NumericVector img, LogicalVector mask,
                             int nx, int ny, int nz,
                             double V, double Rd,
                             NumericVector threshs, NumericVector weights,
                             double dh_legacy, double ZestThr,
                             int connectivity, bool analytical_subthr,
                             Nullable<NumericVector> voxel_weights = R_NilValue) {
  int n = nx * ny * nz;
  int nth = threshs.size();
  double logV = std::log(V);
  double logRd = std::log(Rd);
  double x_zest = -R::pnorm5(ZestThr, 0.0, 1.0, 0, 1);
  bool use_nonstat = voxel_weights.isNotNull();
  NumericVector vw;
  if (use_nonstat) {
    vw = NumericVector(voxel_weights);
  }

  // Build neighbor offsets for 6, 18, or 26 connectivity
  std::vector<int> dx, dy, dz;
  for (int oz = -1; oz <= 1; ++oz) {
    for (int oy = -1; oy <= 1; ++oy) {
      for (int ox = -1; ox <= 1; ++ox) {
        int manhattan = std::abs(ox) + std::abs(oy) + std::abs(oz);
        if (manhattan == 0) continue;
        if (connectivity == 6 && manhattan > 1) continue;
        if (connectivity == 18 && manhattan > 2) continue;
        dx.push_back(ox);
        dy.push_back(oy);
        dz.push_back(oz);
      }
    }
  }
  int n_nb = static_cast<int>(dx.size());

  // Identify first threshold strictly above ZestThr
  int first_cc_hi = nth;
  for (int hi = 0; hi < nth; ++hi) {
    if (threshs[hi] > ZestThr) {
      first_cc_hi = hi;
      break;
    }
  }

  NumericVector accum(n, 0.0);

  if (analytical_subthr) {
    // Exact analytical integral int_0^{min(x_i, x_ZestThr)} x dx = 0.5 * min(x_i, x_ZestThr)^2
    for (int i = 0; i < n; ++i) {
      if (!mask[i]) continue;
      double xi = -R::pnorm5(img[i], 0.0, 1.0, 0, 1);
      double base_x = std::min(xi, x_zest);
      accum[i] = 0.5 * base_x * base_x;
    }
  } else if (first_cc_hi > 0) {
    // Legacy discrete grid below ZestThr: sum -log(Phi(-threshs[hi])) via prefix sums in O(N log N_sub)
    std::vector<double> pref(first_cc_hi + 1, 0.0);
    for (int hi = 0; hi < first_cc_hi; ++hi) {
      double nlp = (hi == 0 && threshs[hi] < -1e6) ? 0.0 : (-R::pnorm5(threshs[hi], 0.0, 1.0, 0, 1));
      pref[hi + 1] = pref[hi] + weights[hi] * nlp;
    }
    for (int i = 0; i < n; ++i) {
      if (!mask[i]) continue;
      double val = img[i];
      int k = static_cast<int>(
        std::lower_bound(threshs.begin(), threshs.begin() + first_cc_hi, val) - threshs.begin()
      );
      accum[i] = pref[k];
    }
  }

  // Descending Union-Find sweep across thresholds > ZestThr
  if (first_cc_hi < nth) {
    double min_cc_th = threshs[first_cc_hi];
    std::vector<int> order;
    order.reserve(static_cast<size_t>(n / 8 + 16));
    for (int i = 0; i < n; ++i) {
      if (mask[i] && img[i] > min_cc_th) {
        order.push_back(i);
      }
    }
    std::sort(order.begin(), order.end(), [&](int a, int b) {
      return img[a] > img[b];
    });

    std::vector<int> parent(n, -1);
    std::vector<int> comp_sz(n, 0);

    auto find_root = [&](int x) -> int {
      int r = x;
      while (parent[r] != r) r = parent[r];
      while (x != r) {
        int nxt = parent[x];
        parent[x] = r;
        x = nxt;
      }
      return r;
    };

    const int nxy = nx * ny;
    size_t ptr = 0;
    const size_t n_ord = order.size();

    if (!use_nonstat) {
      // Exact stationary path (integer cluster sizes)
      std::unordered_map<int, double> log_den_cache;
      std::unordered_map<int, double> step_val_cache;
      log_den_cache.reserve(512);
      step_val_cache.reserve(256);

      for (int hi = nth - 1; hi >= first_cc_hi; --hi) {
        if (hi % 64 == 0) Rcpp::checkUserInterrupt();
        double h = threshs[hi];
        double w = weights[hi];

        while (ptr < n_ord && img[order[ptr]] > h) {
          int idx = order[ptr++];
          parent[idx] = idx;
          comp_sz[idx] = 1;

          int z = idx / nxy;
          int rem = idx % nxy;
          int y = rem / nx;
          int x = rem % nx;

          for (int d = 0; d < n_nb; ++d) {
            int xi = x + dx[d];
            int yi = y + dy[d];
            int zi = z + dz[d];
            if (xi >= 0 && xi < nx && yi >= 0 && yi < ny && zi >= 0 && zi < nz) {
              int nidx = xi + yi * nx + zi * nxy;
              if (parent[nidx] != -1) {
                int ra = find_root(idx);
                int rb = find_root(nidx);
                if (ra != rb) {
                  if (comp_sz[ra] < comp_sz[rb]) std::swap(ra, rb);
                  parent[rb] = ra;
                  comp_sz[ra] += comp_sz[rb];
                }
              }
            }
          }
        }

        if (ptr == 0) continue;

        step_val_cache.clear();
        for (size_t k = 0; k < ptr; ++k) {
          int idx = order[k];
          int r = find_root(idx);
          if (r == idx) {
            int c = comp_sz[r];
            if (step_val_cache.find(c) == step_val_cache.end()) {
              double c23 = std::pow(static_cast<double>(c), 2.0 / 3.0);
              double lden;
              auto it = log_den_cache.find(c);
              if (it == log_den_cache.end()) {
                lden = log_integral_tail(ZestThr, logV, logRd, c23);
                log_den_cache[c] = lden;
              } else {
                lden = it->second;
              }
              double nlp = neglog_pvox_clust_cached(h, ZestThr, logV, logRd, c23, lden);
              step_val_cache[c] = w * nlp;
            }
          }
        }

        for (size_t k = 0; k < ptr; ++k) {
          int idx = order[k];
          int c = comp_sz[find_root(idx)];
          accum[idx] += step_val_cache[c];
        }
      }
    } else {
      // Non-stationary path: continuous resel-weighted cluster mass sum_{v in C} vw[v]
      std::vector<double> comp_mass(n, 0.0);
      std::vector<double> root_step_val(n, 0.0);

      // Precompute 256-node cubic interpolation table for log_den(log(c)) across all valid masses
      double c_min = 1e30, c_max = 0.0;
      for (size_t k = 0; k < n_ord; ++k) {
        double wi = std::max(1e-4, vw[order[k]]);
        if (wi < c_min) c_min = wi;
        c_max += wi;
      }
      if (n_ord > 0 && c_max <= c_min) {
        c_max = c_min * 2.0;
      }

      const int N_DEN_GRID = 256;
      std::vector<double> den_grid(N_DEN_GRID, 0.0);
      double s_min = (n_ord > 0) ? std::log(c_min * 0.95) : -3.0;
      double s_max = (n_ord > 0) ? std::log(c_max * 1.05 + 1.0) : 12.0;
      double ds_den = (s_max - s_min) / static_cast<double>(N_DEN_GRID - 1);
      double inv_ds_den = 1.0 / ds_den;

      if (n_ord > 0) {
        for (int g = 0; g < N_DEN_GRID; ++g) {
          double s = s_min + g * ds_den;
          double c23 = std::exp((2.0 / 3.0) * s);
          den_grid[g] = log_integral_tail(ZestThr, logV, logRd, c23);
        }
      }

      std::vector<int> roots;
      roots.reserve(256);
      const int N_STEP_GRID = 64;
      std::vector<double> step_grid(N_STEP_GRID, 0.0);

      for (int hi = nth - 1; hi >= first_cc_hi; --hi) {
        if (hi % 64 == 0) Rcpp::checkUserInterrupt();
        double h = threshs[hi];
        double w = weights[hi];

        while (ptr < n_ord && img[order[ptr]] > h) {
          int idx = order[ptr++];
          parent[idx] = idx;
          comp_sz[idx] = 1;
          comp_mass[idx] = std::max(1e-4, vw[idx]);

          int z = idx / nxy;
          int rem = idx % nxy;
          int y = rem / nx;
          int x = rem % nx;

          for (int d = 0; d < n_nb; ++d) {
            int xi = x + dx[d];
            int yi = y + dy[d];
            int zi = z + dz[d];
            if (xi >= 0 && xi < nx && yi >= 0 && yi < ny && zi >= 0 && zi < nz) {
              int nidx = xi + yi * nx + zi * nxy;
              if (parent[nidx] != -1) {
                int ra = find_root(idx);
                int rb = find_root(nidx);
                if (ra != rb) {
                  if (comp_sz[ra] < comp_sz[rb]) std::swap(ra, rb);
                  parent[rb] = ra;
                  comp_sz[ra] += comp_sz[rb];
                  comp_mass[ra] += comp_mass[rb];
                }
              }
            }
          }
        }

        if (ptr == 0) continue;

        roots.clear();
        double r_cmin = 1e30, r_cmax = 0.0;
        for (size_t k = 0; k < ptr; ++k) {
          int idx = order[k];
          int r = find_root(idx);
          if (r == idx) {
            roots.push_back(r);
            double cm = comp_mass[r];
            if (cm < r_cmin) r_cmin = cm;
            if (cm > r_cmax) r_cmax = cm;
          }
        }

        if (static_cast<int>(roots.size()) <= N_STEP_GRID || r_cmax <= r_cmin * 1.0001) {
          // Few active clusters: evaluate exact quadrature directly for each root
          for (int r : roots) {
            double cm = comp_mass[r];
            double log_cm = std::log(cm);
            double c23 = std::exp((2.0 / 3.0) * log_cm);
            double lden = interp_cubic_equispaced(den_grid, s_min, inv_ds_den, log_cm);
            double nlp = neglog_pvox_clust_cached(h, ZestThr, logV, logRd, c23, lden);
            root_step_val[r] = w * nlp;
          }
        } else {
          // Many active clusters: build 64-node cubic spline table in log(c) for this h
          double ls_lo = std::log(r_cmin);
          double ls_hi = std::log(r_cmax);
          double ds_step = (ls_hi - ls_lo) / static_cast<double>(N_STEP_GRID - 1);
          double inv_ds_step = 1.0 / ds_step;
          for (int g = 0; g < N_STEP_GRID; ++g) {
            double s = ls_lo + g * ds_step;
            double c23 = std::exp((2.0 / 3.0) * s);
            double lden = interp_cubic_equispaced(den_grid, s_min, inv_ds_den, s);
            step_grid[g] = w * neglog_pvox_clust_cached(h, ZestThr, logV, logRd, c23, lden);
          }
          for (int r : roots) {
            double log_cm = std::log(comp_mass[r]);
            double val = interp_cubic_equispaced(step_grid, ls_lo, inv_ds_step, log_cm);
            root_step_val[r] = std::max(0.0, val);
          }
        }

        for (size_t k = 0; k < ptr; ++k) {
          int idx = order[k];
          accum[idx] += root_step_val[find_root(idx)];
        }
      }
    }
  }

  // Convert accumulated evidence into aggregated -log(p)
  NumericVector out_neglogp(n, 0.0);
  if (analytical_subthr) {
    for (int i = 0; i < n; ++i) {
      if (mask[i] && accum[i] > 0.0) {
        out_neglogp[i] = std::sqrt(2.0 * accum[i]);
      }
    }
  } else {
    for (int i = 0; i < n; ++i) {
      if (mask[i] && accum[i] > 0.0) {
        out_neglogp[i] = 0.5 * (std::sqrt(dh_legacy * (8.0 * accum[i] + dh_legacy)) - dh_legacy);
      }
    }
  }

  return out_neglogp;
}

// Fast cache-friendly C++ smoothness estimation for 2D/3D Z-maps and 4D residuals
// Optionally computes spatially regularized local dLh(v) and normalized voxel weights
// for non-stationary cluster size adjustment (Worsley et al., 1999; Salimi-Khorshidi et al., 2011).
// [[Rcpp::export]]
List smoothest_cpp(NumericVector img, NumericVector mask,
                   int nx, int ny, int nz, int nt,
                   bool compute_local = false, double smooth_fwhm = 3.0,
                   double dof_correction = 1.0) {
  int n_vox = nx * ny * nz;
  bool usez = (nz > 1);
  int nxy = nx * ny;

  std::vector<unsigned char> in_mask(n_vox, 0);
  std::vector<unsigned char> active(n_vox, 0);
  for (int i = 0; i < n_vox; ++i) {
    double m = mask[i];
    if (!NumericVector::is_na(m) && std::isfinite(m) && m > 0.5) {
      in_mask[i] = 1;
      active[i] = 1;
    }
  }

  std::vector<double> mean_v;
  std::vector<double> inv_sd_v;

  if (nt > 1) {
    // 4D residual mode: standardize each voxel along time without copying 4D array
    mean_v.assign(n_vox, 0.0);
    inv_sd_v.assign(n_vox, 0.0);

    // Pass 1: compute temporal means
    for (int t = 0; t < nt; ++t) {
      R_xlen_t off = static_cast<R_xlen_t>(t) * n_vox;
      for (int i = 0; i < n_vox; ++i) {
        if (!active[i]) continue;
        double val = img[off + i];
        if (!std::isfinite(val)) {
          active[i] = 0;
        } else {
          mean_v[i] += val;
        }
      }
    }
    double inv_nt = 1.0 / static_cast<double>(nt);
    for (int i = 0; i < n_vox; ++i) {
      if (active[i]) mean_v[i] *= inv_nt;
    }

    // Pass 2: compute centered sum of squares for numerically stable sample SD
    std::vector<double> ssq(n_vox, 0.0);
    for (int t = 0; t < nt; ++t) {
      R_xlen_t off = static_cast<R_xlen_t>(t) * n_vox;
      for (int i = 0; i < n_vox; ++i) {
        if (!active[i]) continue;
        double diff = img[off + i] - mean_v[i];
        ssq[i] += diff * diff;
      }
    }
    double df_t = static_cast<double>(nt - 1);
    for (int i = 0; i < n_vox; ++i) {
      if (!active[i]) continue;
      double sd = std::sqrt(ssq[i] / df_t);
      if (!std::isfinite(sd) || sd <= 0.0) {
        active[i] = 0;
      } else {
        inv_sd_v[i] = 1.0 / sd;
      }
    }
  } else {
    // 2D/3D mode: guard non-finite image values
    for (int i = 0; i < n_vox; ++i) {
      if (active[i] && !std::isfinite(img[i])) {
        active[i] = 0;
      }
    }
  }

  int mask_volume = 0;
  for (int i = 0; i < n_vox; ++i) {
    if (active[i]) mask_volume++;
  }

  // Collect interior voxels where (x,y,z) and all backward neighbors are active (FSL global rule)
  std::vector<int> interior;
  interior.reserve(mask_volume);
  int zstart = usez ? 1 : 0;
  for (int z = zstart; z < nz; ++z) {
    for (int y = 1; y < ny; ++y) {
      int base = y * nx + z * nxy;
      for (int x = 1; x < nx; ++x) {
        int idx = base + x;
        if (active[idx] &&
            active[idx - 1] &&
            active[idx - nx] &&
            (!usez || active[idx - nxy])) {
          interior.push_back(idx);
        }
      }
    }
  }

  R_xlen_t N = static_cast<R_xlen_t>(interior.size());
  double SSminus[3] = {0.0, 0.0, 0.0};
  double S2[3] = {0.0, 0.0, 0.0};

  // Optional per-axis valid edge lists & voxel accumulators for local smoothness
  std::vector<int> edges_x, edges_y, edges_z;
  std::vector<double> vox_ss_x, vox_ss_y, vox_ss_z;
  std::vector<double> vox_s2_x, vox_s2_y, vox_s2_z;

  if (compute_local) {
    vox_ss_x.assign(n_vox, 0.0);
    vox_ss_y.assign(n_vox, 0.0);
    vox_s2_x.assign(n_vox, 0.0);
    vox_s2_y.assign(n_vox, 0.0);
    if (usez) {
      vox_ss_z.assign(n_vox, 0.0);
      vox_s2_z.assign(n_vox, 0.0);
    }
    edges_x.reserve(mask_volume);
    edges_y.reserve(mask_volume);
    if (usez) edges_z.reserve(mask_volume);

    for (int z = 0; z < nz; ++z) {
      for (int y = 0; y < ny; ++y) {
        int base = y * nx + z * nxy;
        for (int x = 0; x < nx; ++x) {
          int idx = base + x;
          if (!active[idx]) continue;
          if (x >= 1 && active[idx - 1]) edges_x.push_back(idx);
          if (y >= 1 && active[idx - nx]) edges_y.push_back(idx);
          if (usez && z >= 1 && active[idx - nxy]) edges_z.push_back(idx);
        }
      }
    }

    if (nt > 1) {
      // Standardized time series has exact sum_t z^2 == nt - 1 at every active voxel
      double pair_s2 = static_cast<double>(nt - 1);
      for (int idx : edges_x) {
        vox_s2_x[idx] += pair_s2;
        vox_s2_x[idx - 1] += pair_s2;
      }
      for (int idx : edges_y) {
        vox_s2_y[idx] += pair_s2;
        vox_s2_y[idx - nx] += pair_s2;
      }
      if (usez) {
        for (int idx : edges_z) {
          vox_s2_z[idx] += pair_s2;
          vox_s2_z[idx - nxy] += pair_s2;
        }
      }
    }
  }

  if (nt > 1) {
    for (int t = 0; t < nt; ++t) {
      if (t % 32 == 0) Rcpp::checkUserInterrupt();
      R_xlen_t off = static_cast<R_xlen_t>(t) * n_vox;

      for (R_xlen_t k = 0; k < N; ++k) {
        int idx = interior[k];
        double v0 = (img[off + idx] - mean_v[idx]) * inv_sd_v[idx];
        double vx = (img[off + idx - 1] - mean_v[idx - 1]) * inv_sd_v[idx - 1];
        double vy = (img[off + idx - nx] - mean_v[idx - nx]) * inv_sd_v[idx - nx];

        double v0_sq = v0 * v0;
        SSminus[0] += v0 * vx;
        SSminus[1] += v0 * vy;
        S2[0] += 0.5 * (v0_sq + vx * vx);
        S2[1] += 0.5 * (v0_sq + vy * vy);

        if (usez) {
          double vz = (img[off + idx - nxy] - mean_v[idx - nxy]) * inv_sd_v[idx - nxy];
          SSminus[2] += v0 * vz;
          S2[2] += 0.5 * (v0_sq + vz * vz);
        }
      }

      if (compute_local) {
        for (int idx : edges_x) {
          double p = ((img[off + idx] - mean_v[idx]) * inv_sd_v[idx]) *
                     ((img[off + idx - 1] - mean_v[idx - 1]) * inv_sd_v[idx - 1]);
          vox_ss_x[idx] += p;
          vox_ss_x[idx - 1] += p;
        }
        for (int idx : edges_y) {
          double p = ((img[off + idx] - mean_v[idx]) * inv_sd_v[idx]) *
                     ((img[off + idx - nx] - mean_v[idx - nx]) * inv_sd_v[idx - nx]);
          vox_ss_y[idx] += p;
          vox_ss_y[idx - nx] += p;
        }
        if (usez) {
          for (int idx : edges_z) {
            double p = ((img[off + idx] - mean_v[idx]) * inv_sd_v[idx]) *
                       ((img[off + idx - nxy] - mean_v[idx - nxy]) * inv_sd_v[idx - nxy]);
            vox_ss_z[idx] += p;
            vox_ss_z[idx - nxy] += p;
          }
        }
      }
    }
  } else {
    for (R_xlen_t k = 0; k < N; ++k) {
      int idx = interior[k];
      double v0 = img[idx];
      double vx = img[idx - 1];
      double vy = img[idx - nx];

      double v0_sq = v0 * v0;
      SSminus[0] += v0 * vx;
      SSminus[1] += v0 * vy;
      S2[0] += 0.5 * (v0_sq + vx * vx);
      S2[1] += 0.5 * (v0_sq + vy * vy);

      if (usez) {
        double vz = img[idx - nxy];
        SSminus[2] += v0 * vz;
        S2[2] += 0.5 * (v0_sq + vz * vz);
      }
    }

    if (compute_local) {
      for (int idx : edges_x) {
        double v0 = img[idx], v1 = img[idx - 1];
        double p = v0 * v1, s = 0.5 * (v0 * v0 + v1 * v1);
        vox_ss_x[idx] += p; vox_ss_x[idx - 1] += p;
        vox_s2_x[idx] += s; vox_s2_x[idx - 1] += s;
      }
      for (int idx : edges_y) {
        double v0 = img[idx], v1 = img[idx - nx];
        double p = v0 * v1, s = 0.5 * (v0 * v0 + v1 * v1);
        vox_ss_y[idx] += p; vox_ss_y[idx - nx] += p;
        vox_s2_y[idx] += s; vox_s2_y[idx - nx] += s;
      }
      if (usez) {
        for (int idx : edges_z) {
          double v0 = img[idx], v1 = img[idx - nxy];
          double p = v0 * v1, s = 0.5 * (v0 * v0 + v1 * v1);
          vox_ss_z[idx] += p; vox_ss_z[idx - nxy] += p;
          vox_s2_z[idx] += s; vox_s2_z[idx - nxy] += s;
        }
      }
    }
  }

  if (!compute_local) {
    return List::create(
      Named("volume") = mask_volume,
      Named("N") = static_cast<double>(N),
      Named("usez") = usez,
      Named("SSminus") = NumericVector::create(SSminus[0], SSminus[1], SSminus[2]),
      Named("S2") = NumericVector::create(S2[0], S2[1], S2[2])
    );
  }

  // Normalized within-mask Gaussian smoothing of sufficient statistics + empirical-Bayes shrinkage
  smooth_3d_separable(vox_ss_x, nx, ny, nz, smooth_fwhm);
  smooth_3d_separable(vox_s2_x, nx, ny, nz, smooth_fwhm);
  smooth_3d_separable(vox_ss_y, nx, ny, nz, smooth_fwhm);
  smooth_3d_separable(vox_s2_y, nx, ny, nz, smooth_fwhm);
  if (usez) {
    smooth_3d_separable(vox_ss_z, nx, ny, nz, smooth_fwhm);
    smooth_3d_separable(vox_s2_z, nx, ny, nz, smooth_fwhm);
  }

  auto clamp_rho = [](double r) -> double {
    return std::min(0.99999, std::max(1e-4, r));
  };
  double rho_glob_x = (S2[0] > 0.0) ? clamp_rho(std::fabs(SSminus[0] / S2[0])) : 0.5;
  double rho_glob_y = (S2[1] > 0.0) ? clamp_rho(std::fabs(SSminus[1] / S2[1])) : 0.5;
  double rho_glob_z = (usez && S2[2] > 0.0) ? clamp_rho(std::fabs(SSminus[2] / S2[2])) : 0.5;

  double mean_s2_x = 0.0, mean_s2_y = 0.0, mean_s2_z = 0.0;
  if (mask_volume > 0) {
    for (int i = 0; i < n_vox; ++i) {
      if (!active[i]) continue;
      mean_s2_x += vox_s2_x[i];
      mean_s2_y += vox_s2_y[i];
      if (usez) mean_s2_z += vox_s2_z[i];
    }
    mean_s2_x /= mask_volume;
    mean_s2_y /= mask_volume;
    if (usez) mean_s2_z /= mask_volume;
  }
  double tau_x = 0.05 * std::max(1e-8, mean_s2_x);
  double tau_y = 0.05 * std::max(1e-8, mean_s2_y);
  double tau_z = 0.05 * std::max(1e-8, mean_s2_z);

  NumericVector dlh_map(n_vox, 0.0);
  NumericVector voxel_weights(n_vox, 0.0);
  const double scale_dim = usez ? std::pow(8.0, -0.5) : std::pow(4.0, -0.5);
  double sum_dlh = 0.0;

  for (int i = 0; i < n_vox; ++i) {
    if (!active[i]) continue;
    double rx = clamp_rho((std::max(0.0, vox_ss_x[i]) + tau_x * rho_glob_x) / (vox_s2_x[i] + tau_x));
    double ry = clamp_rho((std::max(0.0, vox_ss_y[i]) + tau_y * rho_glob_y) / (vox_s2_y[i] + tau_y));
    double sx2 = -1.0 / (4.0 * std::log(rx));
    double sy2 = -1.0 / (4.0 * std::log(ry));
    double dlh_i;
    if (usez) {
      double rz = clamp_rho((std::max(0.0, vox_ss_z[i]) + tau_z * rho_glob_z) / (vox_s2_z[i] + tau_z));
      double sz2 = -1.0 / (4.0 * std::log(rz));
      dlh_i = std::pow(sx2 * sy2 * sz2, -0.5) * scale_dim * dof_correction;
    } else {
      dlh_i = std::pow(sx2 * sy2, -0.5) * scale_dim * dof_correction;
    }
    dlh_map[i] = dlh_i;
    sum_dlh += dlh_i;
  }

  double mean_dlh = (mask_volume > 0 && sum_dlh > 0.0) ? (sum_dlh / mask_volume) : 1.0;
  double sum_w = 0.0;
  for (int i = 0; i < n_vox; ++i) {
    if (active[i]) {
      // Winsorize extreme local roughness ratios to [0.02, 50.0] for numerical stability
      double w = std::min(50.0, std::max(0.02, dlh_map[i] / mean_dlh));
      voxel_weights[i] = w;
      sum_w += w;
    } else if (in_mask[i]) {
      // Fallback weight 1.0 for mask voxels that had zero temporal variance in 4D residuals
      voxel_weights[i] = 1.0;
    }
  }
  if (mask_volume > 0 && sum_w > 0.0) {
    double norm_w = static_cast<double>(mask_volume) / sum_w;
    for (int i = 0; i < n_vox; ++i) {
      if (active[i]) {
        voxel_weights[i] *= norm_w;
      }
    }
  }

  return List::create(
    Named("volume") = mask_volume,
    Named("N") = static_cast<double>(N),
    Named("usez") = usez,
    Named("SSminus") = NumericVector::create(SSminus[0], SSminus[1], SSminus[2]),
    Named("S2") = NumericVector::create(S2[0], S2[1], S2[2]),
    Named("dLh_map") = dlh_map,
    Named("voxel_weights") = voxel_weights
  );
}
