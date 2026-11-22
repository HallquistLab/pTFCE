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

// Log of expected cluster size Es(h, V, Rd) from D-dimensional GRF theory (D in {2, 3})
inline double log_es_scalar(double h, double logV, double logRd, int D = 3) {
  if (D == 2) {
    if (h >= 0.5) {
      // 2D Euler density: Rd * (2*pi)^(-3/2) * h * exp(-0.5*h^2)
      return logV - logRd + 1.5 * LOG_2PI - std::log(h) + log_mills_norm(h);
    }
    return logV + R::pnorm5(h, 0.0, 1.0, 0, 1);
  }
  if (h >= 1.1) {
    double h2 = h * h;
    // log(V) + log(Phi(-h)) - (log(Rd) + log(h^2 - 1) - 0.5*h^2 - 2*log(2*pi))
    // = log(V) - log(Rd) + 2*log(2*pi) - log(h^2 - 1) + (log(Phi(-h)) + 0.5*h^2)
    return logV - logRd + 2.0 * LOG_2PI - std::log(h2 - 1.0) + log_mills_norm(h);
  }
  return logV + R::pnorm5(h, 0.0, 1.0, 0, 1);
}

// Log of unnormalized posterior density:
// l_D(h) = log(phi(h)) + log(lambda_D(h)) - c^(2/D) * lambda_D(h),
// where c_pow = c^(2/D) (c^1 for D=2, c^(2/3) for D=3).
// Strictly concave on (0.5, Inf) for D=2 and on (h_0, Inf) (h_0 ~ 1.2610) for D=3.
inline double log_post_unnorm(double h, double logV, double logRd, double c_pow, int D = 3) {
  double log_es = log_es_scalar(h, logV, logRd, D);
  double log_lam = (D == 2) ? (-log_es) : ((-2.0 / 3.0) * (log_es - LGAMMA_2_5));
  double lam_cpow = std::exp(log_lam) * c_pow;
  double log_phi = -HALF_LOG_2PI - 0.5 * h * h;
  return log_phi + log_lam - lam_cpow;
}

// 15-point Gauss-Kronrod rule on [a, b] for exp(log_post_unnorm(x) - shift)
struct GKResult {
  double val;
  double err;
};

inline GKResult gk15_panel(double a, double b, double shift,
                           double logV, double logRd, double c_pow, int D = 3) {
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
  double f_center = std::exp(log_post_unnorm(center, logV, logRd, c_pow, D) - shift);
  double res_g = wg[3] * f_center;
  double res_k = wgk[7] * f_center;

  for (int j = 0; j < 7; ++j) {
    double absc = half * xgk[j];
    double fval1 = std::exp(log_post_unnorm(center - absc, logV, logRd, c_pow, D) - shift);
    double fval2 = std::exp(log_post_unnorm(center + absc, logV, logRd, c_pow, D) - shift);
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
                     double logV, double logRd, double c_pow,
                     double tol, int D = 3, int depth = 0) {
  if (b <= a) return 0.0;
  GKResult r = gk15_panel(a, b, shift, logV, logRd, c_pow, D);
  if (r.err <= tol * std::max(1e-16, r.val) || depth >= 18) {
    return r.val;
  }
  double mid = 0.5 * (a + b);
  return adaptive_gk15(a, mid, shift, logV, logRd, c_pow, tol, D, depth + 1) +
         adaptive_gk15(mid, b, shift, logV, logRd, c_pow, tol, D, depth + 1);
}

// Find unique global mode of log_post_unnorm on [a, b_max] (b_max <= Inf)
inline double find_mode_on(double a, double b_max, double logV, double logRd, double c_pow,
                           double &max_val, int D = 3) {
  double fa = log_post_unnorm(a, logV, logRd, c_pow, D);
  double eps = std::min(1e-5, 0.25 * (b_max - a));
  if (eps <= 0.0) {
    max_val = fa;
    return a;
  }
  double fa_eps = log_post_unnorm(a + eps, logV, logRd, c_pow, D);
  if (fa_eps <= fa) {
    // By strict concavity on support, mode on [a, b_max] is at x = a
    max_val = fa;
    return a;
  }
  // Bracket peak
  double lo = a;
  double hi = std::min(b_max, a + 0.5);
  while (hi < b_max && log_post_unnorm(hi, logV, logRd, c_pow, D) >
                       log_post_unnorm(std::max(a, hi - 1e-4), logV, logRd, c_pow, D)) {
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
  double f1 = log_post_unnorm(x1, logV, logRd, c_pow, D);
  double f2 = log_post_unnorm(x2, logV, logRd, c_pow, D);
  for (int it = 0; it < 48; ++it) {
    if (f1 < f2) {
      lo = x1;
      x1 = x2;
      f1 = f2;
      x2 = lo + invphi * (hi - lo);
      f2 = log_post_unnorm(x2, logV, logRd, c_pow, D);
    } else {
      hi = x2;
      x2 = x1;
      f2 = f1;
      x1 = lo + invphi2 * (hi - lo);
      f1 = log_post_unnorm(x1, logV, logRd, c_pow, D);
    }
  }
  double mode = 0.5 * (lo + hi);
  max_val = log_post_unnorm(mode, logV, logRd, c_pow, D);
  double fb = (std::isfinite(b_max)) ? log_post_unnorm(b_max, logV, logRd, c_pow, D) : -INFINITY;
  if (fb > max_val) {
    max_val = fb;
    return b_max;
  }
  return mode;
}

// Compute log( int_a^Inf exp(log_post_unnorm(x)) dx )
double log_integral_tail(double a, double logV, double logRd, double c_pow, int D = 3) {
  double max_f = 0.0;
  double mode = find_mode_on(a, INFINITY, logV, logRd, c_pow, max_f, D);

  // Bracket upper cutoff b > mode where log_post_unnorm(b) <= max_f - 50
  const double drop = 50.0; // exp(-50) ~ 1.9e-22
  double step = 0.25;
  while (log_post_unnorm(mode + step, logV, logRd, c_pow, D) > max_f - drop) {
    step *= 2.0;
  }
  double lo = mode + ((step > 0.25) ? 0.5 * step : 0.0);
  double hi = mode + step;
  for (int it = 0; it < 28; ++it) {
    double mid = 0.5 * (lo + hi);
    if (log_post_unnorm(mid, logV, logRd, c_pow, D) > max_f - drop) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  double b = hi;

  double int_left = (mode > a) ? adaptive_gk15(a, mode, max_f, logV, logRd, c_pow, 1e-12, D) : 0.0;
  double int_right = adaptive_gk15(mode, b, max_f, logV, logRd, c_pow, 1e-12, D);
  return max_f + std::log(int_left + int_right);
}

// Compute log( int_a^b exp(log_post_unnorm(x)) dx ) for finite [a, b]
double log_integral_interval(double a, double b, double logV, double logRd, double c_pow, int D = 3) {
  if (b <= a) return -INFINITY;
  double max_f = 0.0;
  double mode = find_mode_on(a, b, logV, logRd, c_pow, max_f, D);
  double int_l = (mode > a) ? adaptive_gk15(a, mode, max_f, logV, logRd, c_pow, 1e-12, D) : 0.0;
  double int_r = (b > mode) ? adaptive_gk15(mode, b, max_f, logV, logRd, c_pow, 1e-12, D) : 0.0;
  return max_f + std::log(int_l + int_r);
}

// Exact -log(P(H > actH | C = c)) using optional cached log_den
inline double neglog_pvox_clust_cached(double actH, double ZestThr,
                                       double logV, double logRd, double c_pow,
                                       double log_den, int D = 3) {
  if (actH <= ZestThr) {
    return -R::pnorm5(actH, 0.0, 1.0, 0, 1);
  }
  double log_num = log_integral_tail(actH, logV, logRd, c_pow, D);
  double diff = log_den - log_num;
  if (diff < 0.05) {
    // Near ZestThr where P ~ 1: integrate [ZestThr, actH] directly so
    // -log(P) = log(1 + I_left / I_num) > 0 with zero cancellation!
    double log_left = log_integral_interval(ZestThr, actH, logV, logRd, c_pow, D);
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
NumericVector es_cpp(NumericVector h, double V, double Rd, int D = 3) {
  int dim_eff = (D == 2) ? 2 : 3;
  int n = h.size();
  NumericVector out(n);
  double logV = std::log(V);
  double logRd = std::log(Rd);
  for (int i = 0; i < n; ++i) {
    if (NumericVector::is_na(h[i])) {
      out[i] = NA_REAL;
    } else {
      out[i] = std::exp(log_es_scalar(h[i], logV, logRd, dim_eff));
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericVector neglog_pvox_clust_cpp(double V, double Rd, NumericVector c,
                                    NumericVector actH, double ZestThr = 1.3,
                                    int D = 3) {
  int dim_eff = (D == 2) ? 2 : 3;
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
    double c_pow = (dim_eff == 2) ? ci : std::pow(ci, 2.0 / 3.0);
    double lden;
    auto it = den_cache.find(ci);
    if (it == den_cache.end()) {
      lden = log_integral_tail(ZestThr, logV, logRd, c_pow, dim_eff);
      den_cache[ci] = lden;
    } else {
      lden = it->second;
    }
    out[i] = neglog_pvox_clust_cached(hi, ZestThr, logV, logRd, c_pow, lden, dim_eff);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector dvox_clust_cpp(NumericVector h, double V, double Rd, double c,
                             double ZestThr = 1.3, int D = 3) {
  int dim_eff = (D == 2) ? 2 : 3;
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
  double c_pow = (dim_eff == 2) ? c : std::pow(c, 2.0 / 3.0);
  double log_den = log_integral_tail(ZestThr, logV, logRd, c_pow, dim_eff);
  for (int i = 0; i < n; ++i) {
    if (NumericVector::is_na(h[i])) {
      out[i] = NA_REAL;
    } else if (h[i] < ZestThr) {
      out[i] = 0.0;
    } else {
      out[i] = std::exp(log_post_unnorm(h[i], logV, logRd, c_pow, dim_eff) - log_den);
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
  const int dim_eff = (nz > 1) ? 3 : 2;
  const double beta_d = (dim_eff == 2) ? 1.0 : (2.0 / 3.0);
  double logV = std::log(V);
  double logRd = std::log(Rd);
  double x_zest = -R::pnorm5(ZestThr, 0.0, 1.0, 0, 1);
  bool use_nonstat = voxel_weights.isNotNull();
  NumericVector vw;
  if (use_nonstat) {
    vw = NumericVector(voxel_weights);
  }

  // Build neighbor offsets for 6, 18, or 26 connectivity (4 or 8 in 2D when nz == 1)
  std::vector<int> dx, dy, dz;
  int oz_lim = (nz > 1) ? 1 : 0;
  for (int oz = -oz_lim; oz <= oz_lim; ++oz) {
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
              double c_pow = (dim_eff == 2) ? static_cast<double>(c) : std::pow(static_cast<double>(c), 2.0 / 3.0);
              double lden;
              auto it = log_den_cache.find(c);
              if (it == log_den_cache.end()) {
                lden = log_integral_tail(ZestThr, logV, logRd, c_pow, dim_eff);
                log_den_cache[c] = lden;
              } else {
                lden = it->second;
              }
              double nlp = neglog_pvox_clust_cached(h, ZestThr, logV, logRd, c_pow, lden, dim_eff);
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
          double c_pow = std::exp(beta_d * s);
          den_grid[g] = log_integral_tail(ZestThr, logV, logRd, c_pow, dim_eff);
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
            double c_pow = std::exp(beta_d * log_cm);
            double lden = interp_cubic_equispaced(den_grid, s_min, inv_ds_den, log_cm);
            double nlp = neglog_pvox_clust_cached(h, ZestThr, logV, logRd, c_pow, lden, dim_eff);
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
            double c_pow = std::exp(beta_d * s);
            double lden = interp_cubic_equispaced(den_grid, s_min, inv_ds_den, s);
            step_grid[g] = w * neglog_pvox_clust_cached(h, ZestThr, logV, logRd, c_pow, lden, dim_eff);
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
      // For a single 3D Z-map (nt == 1), estimate global noise variance robustly via
      // median(Z^2) / qchisq(0.5, 1) and form local lag-1 covariances from 0.5*(Z_u - Z_v)^2
      // so non-zero mean signal clusters (Z = mu + eps) do not inflate S2 by mu^2 and
      // artificially crush their own local roughness weights.
      std::vector<double> z2_in;
      z2_in.reserve(mask_volume);
      for (int i = 0; i < n_vox; ++i) {
        if (active[i]) z2_in.push_back(img[i] * img[i]);
      }
      double sigma2_rob = 1.0;
      if (!z2_in.empty()) {
        size_t mid_idx = z2_in.size() / 2;
        std::nth_element(z2_in.begin(), z2_in.begin() + mid_idx, z2_in.end());
        double med_z2 = z2_in[mid_idx];
        if (med_z2 > 1e-12) {
          sigma2_rob = med_z2 / 0.4549364231195728; // qchisq(0.5, df = 1)
        }
      }
      for (int idx : edges_x) {
        double d = img[idx] - img[idx - 1];
        double p = sigma2_rob - 0.5 * d * d, s = sigma2_rob;
        vox_ss_x[idx] += p; vox_ss_x[idx - 1] += p;
        vox_s2_x[idx] += s; vox_s2_x[idx - 1] += s;
      }
      for (int idx : edges_y) {
        double d = img[idx] - img[idx - nx];
        double p = sigma2_rob - 0.5 * d * d, s = sigma2_rob;
        vox_ss_y[idx] += p; vox_ss_y[idx - nx] += p;
        vox_s2_y[idx] += s; vox_s2_y[idx - nx] += s;
      }
      if (usez) {
        for (int idx : edges_z) {
          double d = img[idx] - img[idx - nxy];
          double p = sigma2_rob - 0.5 * d * d, s = sigma2_rob;
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

namespace {

// Solve for L2-unit 1D Gaussian kernel whose discrete lag-1 autocorrelation matches
// exp(-1 / (4 * target_sigmasq)), ensuring exact unit variance and target roughness
// across all smoothness regimes (from sub-voxel to ultra-smooth).
void solve_discrete_kernel_1d(double target_sigmasq, int axis_len,
                              std::vector<double> &kern, int &rad) {
  if (axis_len <= 1 || !std::isfinite(target_sigmasq) || target_sigmasq <= 1e-4) {
    rad = 0;
    kern.assign(1, 1.0);
    return;
  }
  double rho_target = std::exp(-1.0 / (4.0 * target_sigmasq));
  if (rho_target <= 1e-7) {
    rad = 0;
    kern.assign(1, 1.0);
    return;
  }
  rho_target = std::min(0.99995, rho_target);

  auto eval_rho = [](double s) -> double {
    int r = std::max(1, std::min(128, static_cast<int>(std::ceil(4.0 * s + 1.0))));
    double inv_s2 = 1.0 / (s * s);
    double s0 = 1.0;
    double s1 = 0.0;
    double w_prev = 1.0; // m = 0
    for (int m = 1; m <= r; ++m) {
      double wm = std::exp(-0.5 * static_cast<double>(m * m) * inv_s2);
      s0 += 2.0 * wm * wm;
      s1 += 2.0 * wm * w_prev;
      w_prev = wm;
    }
    return s1 / s0;
  };

  double s_lo = 0.05;
  double s_hi = std::max(2.0 * std::sqrt(target_sigmasq) + 2.0, 4.0);
  while (s_hi < 64.0 && eval_rho(s_hi) < rho_target) {
    s_hi *= 2.0;
  }
  for (int it = 0; it < 42; ++it) {
    double mid = 0.5 * (s_lo + s_hi);
    if (eval_rho(mid) < rho_target) {
      s_lo = mid;
    } else {
      s_hi = mid;
    }
  }
  double s_star = 0.5 * (s_lo + s_hi);
  rad = std::max(1, std::min(128, static_cast<int>(std::ceil(3.75 * s_star))));
  kern.resize(2 * rad + 1);
  double inv_s2 = 1.0 / (s_star * s_star);
  double l2_sq = 0.0;
  for (int m = -rad; m <= rad; ++m) {
    double wm = std::exp(-0.5 * static_cast<double>(m * m) * inv_s2);
    kern[m + rad] = wm;
    l2_sq += wm * wm;
  }
  double inv_l2 = 1.0 / std::sqrt(l2_sq);
  for (double &wm : kern) wm *= inv_l2;
}

// Convolve 3D white-noise array `w_pad` of dims (px, py, pz) centered around a bounding box
// of dims (bx, by, bz) having symmetric padding (pad_x, pad_y, pad_z) >= (rad_x, rad_y, rad_z),
// writing unit-variance output into `out_bbox` of size bx * by * bz.
void convolve_padded_bbox_3d(const std::vector<double> &w_pad,
                             int px, int py, int pz,
                             int pad_x, int pad_y, int pad_z,
                             int bx, int by, int bz,
                             const std::vector<double> &kx, int rx,
                             const std::vector<double> &ky, int ry,
                             const std::vector<double> &kz, int rz,
                             std::vector<double> &buf1,
                             std::vector<double> &buf2,
                             std::vector<double> &out_bbox) {
  const int pxy = px * py;
  const int y_keep = by + 2 * rz; // wait: keep full py in pass 1, reduce to by in pass 2, bz in pass 3
  (void)y_keep;

  // Pass X: (px, py, pz) -> buf1 of (bx, py, pz) using z in [pad_z - rz, pad_z + bz - 1 + rz],
  // y in [pad_y - ry, pad_y + by - 1 + ry]
  const int wy = by + 2 * ry;
  const int wz = bz + 2 * rz;
  const int y0 = pad_y - ry;
  const int z0 = pad_z - rz;

  buf1.resize(static_cast<size_t>(bx) * wy * wz);
  for (int z = 0; z < wz; ++z) {
    int zp = z0 + z;
    for (int y = 0; y < wy; ++y) {
      int yp = y0 + y;
      const double *src_row = &w_pad[static_cast<size_t>(zp) * pxy + yp * px + pad_x];
      double *dst_row = &buf1[(static_cast<size_t>(z) * wy + y) * bx];
      if (rx == 0) {
        std::copy(src_row, src_row + bx, dst_row);
      } else {
        for (int x = 0; x < bx; ++x) {
          double acc = 0.0;
          const double *s = src_row + x - rx;
          for (int k = 0; k <= 2 * rx; ++k) {
            acc += kx[k] * s[k];
          }
          dst_row[x] = acc;
        }
      }
    }
  }

  // Pass Y: buf1 (bx, wy, wz) -> buf2 (bx, by, wz)
  buf2.resize(static_cast<size_t>(bx) * by * wz);
  for (int z = 0; z < wz; ++z) {
    const double *src_slice = &buf1[static_cast<size_t>(z) * wy * bx];
    double *dst_slice = &buf2[static_cast<size_t>(z) * by * bx];
    if (ry == 0) {
      std::copy(src_slice, src_slice + static_cast<size_t>(bx) * by, dst_slice);
    } else {
      for (int y = 0; y < by; ++y) {
        double *dst_row = dst_slice + y * bx;
        std::fill(dst_row, dst_row + bx, 0.0);
        for (int k = 0; k <= 2 * ry; ++k) {
          const double wk = ky[k];
          const double *src_row = src_slice + (y + k) * bx;
          for (int x = 0; x < bx; ++x) {
            dst_row[x] += wk * src_row[x];
          }
        }
      }
    }
  }

  // Pass Z: buf2 (bx, by, wz) -> out_bbox (bx, by, bz)
  const size_t n_bxy = static_cast<size_t>(bx) * by;
  out_bbox.resize(n_bxy * bz);
  if (rz == 0) {
    std::copy(buf2.begin(), buf2.begin() + n_bxy * bz, out_bbox.begin());
  } else {
    for (int z = 0; z < bz; ++z) {
      double *dst_slice = &out_bbox[static_cast<size_t>(z) * n_bxy];
      std::fill(dst_slice, dst_slice + n_bxy, 0.0);
      for (int k = 0; k <= 2 * rz; ++k) {
        const double wk = kz[k];
        const double *src_slice = &buf2[static_cast<size_t>(z + k) * n_bxy];
        for (size_t i = 0; i < n_bxy; ++i) {
          dst_slice[i] += wk * src_slice[i];
        }
      }
    }
  }
}

// Fast reusable workspace and lazy 2D posterior LUT for computing max(enhanced Z)
// across hundreds of null replicates without numerical quadrature or heap churn inside the loop.
class NullPtfceEngine {
public:
  int nx, ny, nz, n, nxy, dim_eff;
  double beta_d;
  double V, Rd, logV, logRd, ZestThr, x_zest, base_acc_zest;
  int connectivity;
  int sided_mode; // 0 = positive, 1 = negative, 2 = two.sided
  bool use_nonstat;
  bool use_legacy;
  int Nh_eff;
  double logpmin, logpmax_fixed, dh;

  LogicalVector mask_r;
  Nullable<NumericVector> vw_nullable;
  std::vector<int> mask_idx;
  std::vector<double> vw_vec;
  std::vector<int> dx, dy, dz;
  int n_nb;

  // Modern fixed-dh grid and lazy 2D posterior LUT
  std::vector<double> threshs;
  std::vector<double> weights;
  double top_x_grid;

  enum { SMALL_C_MAX = 64, N_C_LUT = 96 };
  double s_min, s_max, ds_lut, inv_ds_lut;
  std::vector<double> small_den; // [1 .. SMALL_C_MAX]
  std::vector<double> lut_den;   // [0 .. N_C_LUT - 1]
  std::vector<unsigned char> row_ready;
  std::vector<double> small_step_lut; // row_hi * (SMALL_C_MAX + 1) + c
  std::vector<double> cubic_step_lut; // row_hi * N_C_LUT + g

  // Reusable Union-Find buffers (reset selectively in O(|order|))
  std::vector<int> order;
  std::vector<int> parent;
  std::vector<int> comp_sz;
  std::vector<double> comp_mass;
  std::vector<double> root_step;
  std::vector<double> accum;

  NullPtfceEngine(const LogicalVector &mask_,
                  int nx_, int ny_, int nz_,
                  double V_, double Rd_, double ZestThr_,
                  int connectivity_, int sided_mode_,
                  bool use_legacy_, int Nh_eff_,
                  double logpmin_, double logpmax_fixed_, double dh_,
                  Nullable<NumericVector> voxel_weights_)
      : nx(nx_), ny(ny_), nz(nz_), n(nx_ * ny_ * nz_), nxy(nx_ * ny_),
        dim_eff((nz_ > 1) ? 3 : 2),
        beta_d((nz_ > 1) ? (2.0 / 3.0) : 1.0),
        V(V_), Rd(Rd_), logV(std::log(V_)), logRd(std::log(Rd_)),
        ZestThr(ZestThr_),
        x_zest(-R::pnorm5(ZestThr_, 0.0, 1.0, 0, 1)),
        base_acc_zest(0.5 * x_zest * x_zest),
        connectivity(connectivity_), sided_mode(sided_mode_),
        use_nonstat(voxel_weights_.isNotNull()),
        use_legacy(use_legacy_), Nh_eff(Nh_eff_),
        logpmin(logpmin_), logpmax_fixed(logpmax_fixed_), dh(dh_),
        mask_r(mask_), vw_nullable(voxel_weights_),
        parent(n, -1), comp_sz(n, 0), accum(n, 0.0) {

    mask_idx.reserve(n / 2 + 16);
    for (int i = 0; i < n; ++i) {
      if (mask_r[i]) mask_idx.push_back(i);
    }

    if (use_nonstat) {
      NumericVector vw_r(voxel_weights_);
      vw_vec.assign(n, 0.0);
      comp_mass.assign(n, 0.0);
      double c_min = 1e30, c_max = 0.0;
      for (int idx : mask_idx) {
        double wi = std::max(1e-4, static_cast<double>(vw_r[idx]));
        vw_vec[idx] = wi;
        if (wi < c_min) c_min = wi;
        c_max += wi;
      }
      if (c_max <= c_min) c_max = c_min * 2.0;
      s_min = std::log(c_min * 0.95);
      s_max = std::log(c_max * 1.05 + 1.0);
    } else {
      s_min = std::log(static_cast<double>(SMALL_C_MAX) * 0.95);
      s_max = std::log(std::max(static_cast<double>(SMALL_C_MAX) * 2.0, V * 1.05 + 2.0));
    }
    ds_lut = (s_max - s_min) / static_cast<double>(N_C_LUT - 1);
    inv_ds_lut = 1.0 / ds_lut;

    root_step.assign(n, 0.0);
    order.reserve(mask_idx.size() / 4 + 32);

    int oz_lim = (nz > 1) ? 1 : 0;
    for (int oz = -oz_lim; oz <= oz_lim; ++oz) {
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
    n_nb = static_cast<int>(dx.size());

    if (!use_legacy) {
      // Precompute normalizing denominators log_den once
      if (!use_nonstat) {
        small_den.assign(SMALL_C_MAX + 1, 0.0);
        for (int c = 1; c <= SMALL_C_MAX; ++c) {
          double c_pow = (dim_eff == 2) ? static_cast<double>(c) : std::pow(static_cast<double>(c), 2.0 / 3.0);
          small_den[c] = log_integral_tail(ZestThr, logV, logRd, c_pow, dim_eff);
        }
      }
      lut_den.assign(N_C_LUT, 0.0);
      for (int g = 0; g < N_C_LUT; ++g) {
        double s = s_min + g * ds_lut;
        double c_pow = std::exp(beta_d * s);
        lut_den[g] = log_integral_tail(ZestThr, logV, logRd, c_pow, dim_eff);
      }
      top_x_grid = x_zest;
      extend_modern_grid_to(30.0); // initially cover up to -log(p) = 30 (Z ~ 7.4)
    }
  }

  void extend_modern_grid_to(double x_max) {
    if (x_max <= top_x_grid) return;
    // Rebuild modern grid up to max(x_max + 2.0, top_x_grid) using exact ptfce() rule
    double target_x = std::max(x_max + 2.0, 30.0);
    std::vector<double> new_m;
    std::vector<double> new_w;

    double x_fine_cap = std::min(target_x, 100.0);
    int n_lo = std::max(1, static_cast<int>(std::ceil((x_fine_cap - x_zest) / dh)));
    new_m.reserve(n_lo + 205);
    new_w.reserve(n_lo + 205);
    for (int i = 1; i <= n_lo; ++i) {
      new_m.push_back(x_zest + (static_cast<double>(i) - 0.5) * dh);
      new_w.push_back(dh);
    }
    double top_lo = x_zest + n_lo * dh;
    if (target_x > top_lo) {
      const int n_hi = 200;
      double dh_hi = (target_x - top_lo) / static_cast<double>(n_hi);
      for (int i = 1; i <= n_hi; ++i) {
        new_m.push_back(top_lo + (static_cast<double>(i) - 0.5) * dh_hi);
        new_w.push_back(dh_hi);
      }
    }

    int old_nth = static_cast<int>(threshs.size());
    int new_nth = static_cast<int>(new_m.size());
    // Note: if previous grid stayed <= 100.0, all existing prefix points are identical!
    bool prefix_identical = (top_x_grid <= 100.0) && (old_nth <= new_nth);

    threshs.resize(new_nth);
    weights = new_w;
    int start_conv = prefix_identical ? old_nth : 0;
    for (int i = start_conv; i < new_nth; ++i) {
      threshs[i] = R::qnorm5(-new_m[i], 0.0, 1.0, 0, 1);
    }

    if (prefix_identical) {
      row_ready.resize(new_nth, 0);
      if (!use_nonstat) {
        small_step_lut.resize(static_cast<size_t>(new_nth) * (SMALL_C_MAX + 1), 0.0);
      }
      cubic_step_lut.resize(static_cast<size_t>(new_nth) * N_C_LUT, 0.0);
    } else {
      row_ready.assign(new_nth, 0);
      if (!use_nonstat) {
        small_step_lut.assign(static_cast<size_t>(new_nth) * (SMALL_C_MAX + 1), 0.0);
      }
      cubic_step_lut.assign(static_cast<size_t>(new_nth) * N_C_LUT, 0.0);
    }
    top_x_grid = new_m.back() + 0.5 * new_w.back();
  }

  inline void ensure_lut_row(int hi) {
    if (row_ready[hi]) return;
    double h = threshs[hi];
    double w = weights[hi];
    if (!use_nonstat) {
      double *s_row = &small_step_lut[static_cast<size_t>(hi) * (SMALL_C_MAX + 1)];
      for (int c = 1; c <= SMALL_C_MAX; ++c) {
        double c_pow = (dim_eff == 2) ? static_cast<double>(c) : std::pow(static_cast<double>(c), 2.0 / 3.0);
        s_row[c] = w * neglog_pvox_clust_cached(h, ZestThr, logV, logRd, c_pow, small_den[c], dim_eff);
      }
    }
    double *c_row = &cubic_step_lut[static_cast<size_t>(hi) * N_C_LUT];
    for (int g = 0; g < N_C_LUT; ++g) {
      double s = s_min + g * ds_lut;
      double c_pow = std::exp(beta_d * s);
      c_row[g] = w * neglog_pvox_clust_cached(h, ZestThr, logV, logRd, c_pow, lut_den[g], dim_eff);
    }
    row_ready[hi] = 1;
  }

  inline double interp_row_cubic(int hi, double log_cm) const {
    const double *y = &cubic_step_lut[static_cast<size_t>(hi) * N_C_LUT];
    double u = (log_cm - s_min) * inv_ds_lut;
    if (u <= 0.0) return y[0];
    if (u >= static_cast<double>(N_C_LUT - 1)) return y[N_C_LUT - 1];
    int i1 = static_cast<int>(std::floor(u));
    if (i1 < 1) i1 = 1;
    if (i1 > N_C_LUT - 3) i1 = N_C_LUT - 3;
    double t = u - static_cast<double>(i1);
    double w0 = -t * (t - 1.0) * (t - 2.0) * (1.0 / 6.0);
    double w1 =  (t + 1.0) * (t - 1.0) * (t - 2.0) * 0.5;
    double w2 = -(t + 1.0) * t * (t - 2.0) * 0.5;
    double w3 =  (t + 1.0) * t * (t - 1.0) * (1.0 / 6.0);
    return std::max(0.0, w0 * y[i1 - 1] + w1 * y[i1] + w2 * y[i1 + 1] + w3 * y[i1 + 2]);
  }

  inline int find_root(int x) {
    int r = x;
    while (parent[r] != r) r = parent[r];
    while (x != r) {
      int nxt = parent[x];
      parent[x] = r;
      x = nxt;
    }
    return r;
  }

  // Evaluate max enhanced Z on a single tail (`sign_mult` = +1.0 for positive tail, -1.0 for negative tail).
  // `candidates` is an optional pre-filtered list of in-mask indices with sign_mult * z_vol[idx] > ZestThr;
  // if nullptr, all `mask_idx` are scanned.
  double eval_single_tail_modern(const std::vector<double> &z_vol, double sign_mult,
                                 const std::vector<int> *candidates, double sub_fallback) {
    double z_peak = -INFINITY;
    if (candidates == nullptr) {
      for (int idx : mask_idx) {
        double v = sign_mult * z_vol[idx];
        if (v > z_peak) z_peak = v;
      }
    } else if (!candidates->empty()) {
      for (int idx : *candidates) {
        double v = sign_mult * z_vol[idx];
        if (v > z_peak) z_peak = v;
      }
    } else {
      return sub_fallback;
    }

    if (z_peak <= ZestThr) {
      return z_peak;
    }

    double x_max = -R::pnorm5(z_peak, 0.0, 1.0, 0, 1);
    if (x_max > top_x_grid) {
      extend_modern_grid_to(x_max);
    }

    // Highest threshold strictly below z_peak
    int top_hi = static_cast<int>(
      std::lower_bound(threshs.begin(), threshs.end(), z_peak) - threshs.begin()
    ) - 1;

    if (top_hi < 0) {
      return ZestThr;
    }

    double min_cc_th = threshs[0];
    order.clear();
    if (candidates == nullptr) {
      for (int idx : mask_idx) {
        if (sign_mult * z_vol[idx] > min_cc_th) {
          order.push_back(idx);
        }
      }
    } else {
      for (int idx : *candidates) {
        if (sign_mult * z_vol[idx] > min_cc_th) {
          order.push_back(idx);
        }
      }
    }

    if (order.empty()) {
      return ZestThr;
    }

    if (sign_mult > 0.0) {
      std::sort(order.begin(), order.end(), [&](int a, int b) {
        return z_vol[a] > z_vol[b];
      });
    } else {
      std::sort(order.begin(), order.end(), [&](int a, int b) {
        return z_vol[a] < z_vol[b];
      });
    }

    const size_t n_ord = order.size();
    for (size_t k = 0; k < n_ord; ++k) {
      accum[order[k]] = base_acc_zest;
    }

    size_t ptr = 0;
    if (!use_nonstat) {
      for (int hi = top_hi; hi >= 0; --hi) {
        double h = threshs[hi];
        while (ptr < n_ord && (sign_mult * z_vol[order[ptr]]) > h) {
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

        ensure_lut_row(hi);
        const double *s_row = &small_step_lut[static_cast<size_t>(hi) * (SMALL_C_MAX + 1)];

        for (size_t k = 0; k < ptr; ++k) {
          int idx = order[k];
          int r = find_root(idx); // compresses parent[idx] -> r directly
          if (r == idx) {
            int c = comp_sz[r];
            root_step[r] = (c <= SMALL_C_MAX) ? s_row[c] : interp_row_cubic(hi, std::log(static_cast<double>(c)));
          }
        }
        for (size_t k = 0; k < ptr; ++k) {
          int idx = order[k];
          accum[idx] += root_step[parent[idx]];
        }
      }
    } else {
      for (int hi = top_hi; hi >= 0; --hi) {
        double h = threshs[hi];
        while (ptr < n_ord && (sign_mult * z_vol[order[ptr]]) > h) {
          int idx = order[ptr++];
          parent[idx] = idx;
          comp_sz[idx] = 1;
          comp_mass[idx] = vw_vec[idx];

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

        ensure_lut_row(hi);
        for (size_t k = 0; k < ptr; ++k) {
          int idx = order[k];
          int r = find_root(idx);
          if (r == idx) {
            root_step[r] = interp_row_cubic(hi, std::log(comp_mass[r]));
          }
        }
        for (size_t k = 0; k < ptr; ++k) {
          int idx = order[k];
          accum[idx] += root_step[parent[idx]];
        }
      }
    }

    double max_acc = base_acc_zest;
    for (size_t k = 0; k < n_ord; ++k) {
      int idx = order[k];
      if (accum[idx] > max_acc) max_acc = accum[idx];
      parent[idx] = -1; // selective O(|order|) reset
    }
    double max_nlp = std::sqrt(2.0 * max_acc);
    return R::qnorm5(-max_nlp, 0.0, 1.0, 0, 1);
  }

  // Exact legacy fallback for `use_legacy = true`
  double eval_legacy_max(const std::vector<double> &z_vol) {
    double z_peak = -INFINITY;
    for (int idx : mask_idx) {
      double v = (sided_mode == 0) ? z_vol[idx] :
                 (sided_mode == 1) ? -z_vol[idx] : std::fabs(z_vol[idx]);
      if (v > z_peak) z_peak = v;
    }
    double lpmax = R_FINITE(logpmax_fixed) ? logpmax_fixed : (-R::pnorm5(z_peak, 0.0, 1.0, 0, 1));
    if (!std::isfinite(lpmax) || lpmax <= logpmin) lpmax = logpmin + 1e-6;
    NumericVector th(Nh_eff);
    NumericVector wt(Nh_eff, 1.0);
    double dh_leg = (lpmax - logpmin) / static_cast<double>(Nh_eff - 1);
    for (int i = 0; i < Nh_eff; ++i) {
      double lp = logpmin + i * dh_leg;
      th[i] = (i == 0 && logpmin == 0.0) ? -9999999999.0 : R::qnorm5(-lp, 0.0, 1.0, 0, 1);
    }

    auto max_from_tail = [&](double sgn) -> double {
      NumericVector img_r(n, 0.0);
      for (int idx : mask_idx) img_r[idx] = sgn * z_vol[idx];
      NumericVector nlp = ptfce_tail_cpp(img_r, mask_r, nx, ny, nz, V, Rd,
                                         th, wt, dh_leg, ZestThr, connectivity,
                                         false, vw_nullable);
      double max_nlp = -1.0;
      for (int idx : mask_idx) {
        if (sided_mode == 2 && (sgn * z_vol[idx]) <= 0.0) continue;
        if (nlp[idx] > max_nlp) max_nlp = nlp[idx];
      }
      if (max_nlp <= 0.0) {
        return (sided_mode == 2) ? 0.0 : R::qnorm5(1.0 - 1.11e-16, 0.0, 1.0, 0, 0);
      }
      double z = R::qnorm5(-max_nlp, 0.0, 1.0, 0, 1);
      return (sided_mode == 2) ? std::max(0.0, z) : z;
    };

    if (sided_mode == 0) return max_from_tail(1.0);
    if (sided_mode == 1) return max_from_tail(-1.0);
    return std::max(max_from_tail(1.0), max_from_tail(-1.0));
  }

  double eval_max(const std::vector<double> &z_vol,
                  const std::vector<int> *cand_pos = nullptr,
                  const std::vector<int> *cand_neg = nullptr,
                  double sub_pos = 0.0, double sub_neg = 0.0) {
    if (use_legacy) {
      return eval_legacy_max(z_vol);
    }
    if (sided_mode == 0) {
      return eval_single_tail_modern(z_vol, 1.0, cand_pos, sub_pos);
    }
    if (sided_mode == 1) {
      return eval_single_tail_modern(z_vol, -1.0, cand_neg, sub_neg);
    }
    // Two-sided max(|Z|): note if all voxels have |Z| <= ZestThr, max(|Z|) >= 0
    double mp = eval_single_tail_modern(z_vol, 1.0, cand_pos, std::max(0.0, sub_pos));
    double mn = eval_single_tail_modern(z_vol, -1.0, cand_neg, std::max(0.0, sub_neg));
    return std::max(0.0, std::max(mp, mn));
  }
};

} // namespace

// Parametric Gaussian Random Field Monte Carlo FWER calibration (stationary & non-stationary)
// [[Rcpp::export]]
NumericVector calibrate_fwer_grf_cpp(LogicalVector mask,
                                     int nx, int ny, int nz,
                                     double V, double Rd,
                                     NumericVector sigmasq_base,
                                     double ZestThr,
                                     int connectivity,
                                     int sided_mode,
                                     bool use_legacy,
                                     int Nh_eff,
                                     double logpmin,
                                     double logpmax_fixed,
                                     double dh,
                                     int n_perm,
                                     Nullable<NumericVector> voxel_weights = R_NilValue) {
  RNGScope rng_scope;

  NullPtfceEngine engine(mask, nx, ny, nz, V, Rd, ZestThr,
                         connectivity, sided_mode, use_legacy, Nh_eff,
                         logpmin, logpmax_fixed, dh, voxel_weights);

  if (engine.mask_idx.empty() || n_perm <= 0) {
    return NumericVector(std::max(0, n_perm), NA_REAL);
  }

  const int nxy = nx * ny;
  // Tight 3D bounding box around active mask voxels
  int xmin = nx, xmax = -1, ymin = ny, ymax = -1, zmin = nz, zmax = -1;
  for (int idx : engine.mask_idx) {
    int z = idx / nxy;
    int rem = idx % nxy;
    int y = rem / nx;
    int x = rem % nx;
    if (x < xmin) xmin = x;
    if (x > xmax) xmax = x;
    if (y < ymin) ymin = y;
    if (y > ymax) ymax = y;
    if (z < zmin) zmin = z;
    if (z > zmax) zmax = z;
  }
  const int bx = xmax - xmin + 1;
  const int by = ymax - ymin + 1;
  const int bz = zmax - zmin + 1;
  const int n_mask = static_cast<int>(engine.mask_idx.size());

  // Map each active mask voxel to its linear offset inside the tight bounding box
  std::vector<size_t> mask_bbox_off(n_mask);
  for (int m = 0; m < n_mask; ++m) {
    int idx = engine.mask_idx[m];
    int z = idx / nxy - zmin;
    int rem = idx % nxy;
    int y = rem / nx - ymin;
    int x = rem % nx - xmin;
    mask_bbox_off[m] = (static_cast<size_t>(z) * by + y) * bx + x;
  }

  // Scale base sigmasq so its exact determinant roughness matches Rd / V
  bool usez = (nz > 1);
  double dim_d = usez ? 3.0 : 2.0;
  double sx2 = std::max(1e-6, static_cast<double>(sigmasq_base[0]));
  double sy2 = std::max(1e-6, static_cast<double>(sigmasq_base[1]));
  double sz2 = usez ? std::max(1e-6, static_cast<double>(sigmasq_base[2])) : 0.0;
  double dlh_target = Rd / V;
  double dlh_raw = usez
      ? (std::pow(sx2 * sy2 * sz2, -0.5) * std::pow(8.0, -0.5))
      : (std::pow(sx2 * sy2, -0.5) * std::pow(4.0, -0.5));
  if (std::isfinite(dlh_raw) && dlh_raw > 0.0 && dlh_target > 0.0) {
    double scale_s2 = std::pow(dlh_raw / dlh_target, 2.0 / dim_d);
    sx2 *= scale_s2;
    sy2 *= scale_s2;
    if (usez) sz2 *= scale_s2;
  }

  // Determine whether multi-scale non-stationary synthesis is required
  bool multi_scale = false;
  double u_min = 0.0, u_max = 0.0;
  std::vector<double> u_vox;
  if (engine.use_nonstat) {
    u_vox.resize(n_mask);
    for (int m = 0; m < n_mask; ++m) {
      double wi = std::min(50.0, std::max(0.02, engine.vw_vec[engine.mask_idx[m]]));
      u_vox[m] = -(1.0 / dim_d) * std::log(wi); // log(sigma(v) / sigma_0)
    }
    // Use [0.5%, 99.5%] quantiles of local log-scale so isolated boundary-winsorized
    // voxels do not blow up 3D padding radius across the entire brain bounding box
    std::vector<double> u_sort = u_vox;
    int idx_lo = static_cast<int>(0.005 * static_cast<double>(n_mask));
    int idx_hi = std::max(idx_lo, n_mask - 1 - idx_lo);
    std::nth_element(u_sort.begin(), u_sort.begin() + idx_lo, u_sort.end());
    u_min = u_sort[idx_lo];
    std::nth_element(u_sort.begin() + idx_lo, u_sort.begin() + idx_hi, u_sort.end());
    u_max = u_sort[idx_hi];
    // Cap extreme kernel bloat at +/- log(2.5) around median while preserving 99%+ range
    u_min = std::max(u_min, -1.10);
    u_max = std::min(u_max,  1.10);
    if (u_max - u_min > 1e-3) {
      multi_scale = true;
    }
  }

  std::vector<double> z_vol(engine.n, 0.0);
  std::vector<double> buf1, buf2;
  std::vector<int> cand_pos, cand_neg, touched;
  if (!use_legacy) {
    cand_pos.reserve(n_mask / 6 + 32);
    cand_neg.reserve(n_mask / 6 + 32);
    touched.reserve(n_mask / 3 + 64);
  }
  const double z_cut = 0.99999 * ZestThr;
  NumericVector max_dist(n_perm);

  auto eval_from_bbox_vals = [&](auto &&get_val) -> double {
    if (use_legacy) {
      for (int m = 0; m < n_mask; ++m) {
        z_vol[engine.mask_idx[m]] = get_val(m);
      }
      return engine.eval_max(z_vol);
    }
    cand_pos.clear();
    cand_neg.clear();
    touched.clear();
    double sub_pos = -INFINITY;
    double sub_neg = -INFINITY;
    for (int m = 0; m < n_mask; ++m) {
      double val = get_val(m);
      if (sided_mode != 1) {
        if (val > sub_pos) sub_pos = val;
        if (val > z_cut) {
          int idx = engine.mask_idx[m];
          z_vol[idx] = val;
          cand_pos.push_back(idx);
          touched.push_back(idx);
        }
      }
      if (sided_mode != 0) {
        double nval = -val;
        if (nval > sub_neg) sub_neg = nval;
        if (nval > z_cut) {
          int idx = engine.mask_idx[m];
          z_vol[idx] = val;
          cand_neg.push_back(idx);
          touched.push_back(idx);
        }
      }
    }
    double out_m = engine.eval_max(z_vol, &cand_pos, &cand_neg, sub_pos, sub_neg);
    for (int idx : touched) {
      z_vol[idx] = 0.0;
    }
    return out_m;
  };

  if (!multi_scale) {
    // Stationary anisotropic 3D Gaussian filter on tight padded mask bounding box
    std::vector<double> kx, ky, kz;
    int rx = 0, ry = 0, rz = 0;
    solve_discrete_kernel_1d(sx2, nx, kx, rx);
    solve_discrete_kernel_1d(sy2, ny, ky, ry);
    solve_discrete_kernel_1d(sz2, nz, kz, rz);

    const int px = bx + 2 * rx;
    const int py = by + 2 * ry;
    const int pz = bz + 2 * rz;
    const size_t n_pad = static_cast<size_t>(px) * py * pz;
    std::vector<double> w_pad(n_pad);
    std::vector<double> g_bbox;

    for (int b = 0; b < n_perm; ++b) {
      if (b % 32 == 0) Rcpp::checkUserInterrupt();
      for (size_t i = 0; i < n_pad; ++i) {
        w_pad[i] = R::norm_rand();
      }
      convolve_padded_bbox_3d(w_pad, px, py, pz, rx, ry, rz, bx, by, bz,
                              kx, rx, ky, ry, kz, rz, buf1, buf2, g_bbox);
      max_dist[b] = eval_from_bbox_vals([&](int m) -> double {
        return g_bbox[mask_bbox_off[m]];
      });
    }
  } else {
    // Multi-scale separable Gaussian filter bank (K = 4 log-spaced scales) driven by
    // shared white noise with exact closed-form cross-scale variance normalization
    const int K = 4;
    double du = (u_max - u_min) / static_cast<double>(K - 1);
    std::vector<std::vector<double>> kx_k(K), ky_k(K), kz_k(K);
    std::vector<int> rx_k(K), ry_k(K), rz_k(K);
    int max_rx = 0, max_ry = 0, max_rz = 0;

    for (int k = 0; k < K; ++k) {
      double alpha2 = std::exp(2.0 * (u_min + k * du));
      solve_discrete_kernel_1d(sx2 * alpha2, nx, kx_k[k], rx_k[k]);
      solve_discrete_kernel_1d(sy2 * alpha2, ny, ky_k[k], ry_k[k]);
      solve_discrete_kernel_1d(sz2 * alpha2, nz, kz_k[k], rz_k[k]);
      if (rx_k[k] > max_rx) max_rx = rx_k[k];
      if (ry_k[k] > max_ry) max_ry = ry_k[k];
      if (rz_k[k] > max_rz) max_rz = rz_k[k];
    }

    auto kern_dot = [](const std::vector<double> &ka, int ra,
                       const std::vector<double> &kb, int rb) -> double {
      int r_common = std::min(ra, rb);
      double acc = 0.0;
      for (int m = -r_common; m <= r_common; ++m) {
        acc += ka[m + ra] * kb[m + rb];
      }
      return acc;
    };

    std::vector<double> cross_corr(K - 1, 1.0);
    for (int k = 0; k < K - 1; ++k) {
      cross_corr[k] = kern_dot(kx_k[k], rx_k[k], kx_k[k + 1], rx_k[k + 1]) *
                      kern_dot(ky_k[k], ry_k[k], ky_k[k + 1], ry_k[k + 1]) *
                      kern_dot(kz_k[k], rz_k[k], kz_k[k + 1], rz_k[k + 1]);
    }

    // Precompute exact unit-variance blending coefficients (k0[m], c0[m], c1[m]) per mask voxel
    std::vector<int> knot_lo(n_mask);
    std::vector<double> coef0(n_mask), coef1(n_mask);
    double inv_du = 1.0 / du;
    for (int m = 0; m < n_mask; ++m) {
      double u_clamped = std::min(u_max, std::max(u_min, u_vox[m]));
      double pos = (u_clamped - u_min) * inv_du;
      int k = static_cast<int>(std::floor(pos));
      if (k < 0) k = 0;
      if (k > K - 2) k = K - 2;
      double t = std::min(1.0, std::max(0.0, pos - static_cast<double>(k)));
      double om_t = 1.0 - t;
      double var_blend = om_t * om_t + t * t + 2.0 * t * om_t * cross_corr[k];
      double inv_sd = 1.0 / std::sqrt(std::max(1e-12, var_blend));
      knot_lo[m] = k;
      coef0[m] = om_t * inv_sd;
      coef1[m] = t * inv_sd;
    }

    const int px = bx + 2 * max_rx;
    const int py = by + 2 * max_ry;
    const int pz = bz + 2 * max_rz;
    const size_t n_pad = static_cast<size_t>(px) * py * pz;
    std::vector<double> w_pad(n_pad);
    std::vector<std::vector<double>> g_scales(K);

    for (int b = 0; b < n_perm; ++b) {
      if (b % 16 == 0) Rcpp::checkUserInterrupt();
      for (size_t i = 0; i < n_pad; ++i) {
        w_pad[i] = R::norm_rand();
      }
      for (int k = 0; k < K; ++k) {
        convolve_padded_bbox_3d(w_pad, px, py, pz, max_rx, max_ry, max_rz,
                                bx, by, bz,
                                kx_k[k], rx_k[k], ky_k[k], ry_k[k], kz_k[k], rz_k[k],
                                buf1, buf2, g_scales[k]);
      }
      max_dist[b] = eval_from_bbox_vals([&](int m) -> double {
        int k = knot_lo[m];
        size_t boff = mask_bbox_off[m];
        return coef0[m] * g_scales[k][boff] + coef1[m] * g_scales[k + 1][boff];
      });
    }
  }

  return max_dist;
}

// Fast Rademacher sign-flipping FWER calibration for 4D subject COPEs or 4D GLM residuals
// [[Rcpp::export]]
List calibrate_fwer_signflip_cpp(NumericVector data4d,
                                 LogicalVector mask,
                                 int nx, int ny, int nz, int nt,
                                 bool center_residuals,
                                 double V, double Rd, double ZestThr,
                                 int connectivity, int sided_mode,
                                 bool use_legacy, int Nh_eff,
                                 double logpmin, double logpmax_fixed,
                                 double dh, int n_perm,
                                 Nullable<NumericVector> voxel_weights = R_NilValue,
                                 double dof_hint = NA_REAL) {
  RNGScope rng_scope;

  NullPtfceEngine engine(mask, nx, ny, nz, V, Rd, ZestThr,
                         connectivity, sided_mode, use_legacy, Nh_eff,
                         logpmin, logpmax_fixed, dh, voxel_weights);

  const R_xlen_t n_vox = static_cast<R_xlen_t>(nx) * ny * nz;
  const int n_mask = static_cast<int>(engine.mask_idx.size());

  // Extract active finite non-constant mask voxels and unit-normalize each voxel's nt-vector
  std::vector<int> active_vox;
  active_vox.reserve(n_mask);
  std::vector<double> tmp_t(nt);
  std::vector<double> rows_tmp;
  rows_tmp.reserve(static_cast<size_t>(n_mask) * nt);

  double ar1_sum = 0.0;
  int ar1_cnt = 0;

  for (int m = 0; m < n_mask; ++m) {
    int idx = engine.mask_idx[m];
    bool ok = true;
    double sum_v = 0.0;
    for (int t = 0; t < nt; ++t) {
      double val = data4d[static_cast<R_xlen_t>(t) * n_vox + idx];
      if (!std::isfinite(val)) {
        ok = false;
        break;
      }
      tmp_t[t] = val;
      sum_v += val;
    }
    if (!ok) continue;

    if (center_residuals) {
      double mu = sum_v / static_cast<double>(nt);
      for (int t = 0; t < nt; ++t) tmp_t[t] -= mu;
    }

    double ssq = 0.0;
    for (int t = 0; t < nt; ++t) ssq += tmp_t[t] * tmp_t[t];
    if (!std::isfinite(ssq) || ssq <= 1e-30) continue;

    double inv_norm = 1.0 / std::sqrt(ssq);
    for (int t = 0; t < nt; ++t) tmp_t[t] *= inv_norm;

    // Track lag-1 autocorrelation along 4th dimension (centered) to diagnose unwhitened TRs
    double mu_u = 0.0;
    if (!center_residuals) {
      for (int t = 0; t < nt; ++t) mu_u += tmp_t[t];
      mu_u /= static_cast<double>(nt);
    }
    double num_ar = 0.0, den_ar = 0.0;
    for (int t = 0; t < nt; ++t) {
      double d0 = tmp_t[t] - mu_u;
      den_ar += d0 * d0;
      if (t > 0) num_ar += d0 * (tmp_t[t - 1] - mu_u);
    }
    if (den_ar > 1e-16) {
      ar1_sum += num_ar / den_ar;
      ar1_cnt++;
    }

    active_vox.push_back(idx);
    for (int t = 0; t < nt; ++t) rows_tmp.push_back(tmp_t[t]);
  }

  const int n_act = static_cast<int>(active_vox.size());
  if (n_act == 0) {
    stop("4D input for `sign_flip` contains zero finite non-constant voxels inside `mask`.");
  }
  double mean_temp_ar1 = (ar1_cnt > 0) ? (ar1_sum / static_cast<double>(ar1_cnt)) : 0.0;

  // Transpose into contiguous column-major `U_col[t * n_act + a]` for SIMD DAXPY speed
  std::vector<double> U_col(static_cast<size_t>(nt) * n_act);
  for (int a = 0; a < n_act; ++a) {
    const double *row = &rows_tmp[static_cast<size_t>(a) * nt];
    for (int t = 0; t < nt; ++t) {
      U_col[static_cast<size_t>(t) * n_act + a] = row[t];
    }
  }

  // If `center_residuals`, identify an exact orthonormal basis Q (nt x nu_sub) of the
  // empirical GLM residual row-space via two-pass Modified Gram-Schmidt
  std::vector<double> Q;
  int nu_sub = nt;
  if (center_residuals) {
    Q.reserve(static_cast<size_t>(nt) * nt);
    nu_sub = 0;
    int max_rank = nt - 1;
    int stride = std::max(1, n_act / std::max(8 * nt, 512));
    std::vector<double> v(nt);

    for (int pass = 0; pass < 2 && nu_sub < max_rank; ++pass) {
      int step = (pass == 0) ? stride : 1;
      for (int a = 0; a < n_act && nu_sub < max_rank; a += step) {
        const double *row = &rows_tmp[static_cast<size_t>(a) * nt];
        std::copy(row, row + nt, v.begin());
        // Two-pass Gram-Schmidt orthogonalization for machine-precision orthogonality
        for (int gs = 0; gs < 2; ++gs) {
          for (int j = 0; j < nu_sub; ++j) {
            const double *qj = &Q[static_cast<size_t>(j) * nt];
            double dot = 0.0;
            for (int t = 0; t < nt; ++t) dot += qj[t] * v[t];
            for (int t = 0; t < nt; ++t) v[t] -= dot * qj[t];
          }
        }
        double nrm2 = 0.0;
        for (int t = 0; t < nt; ++t) nrm2 += v[t] * v[t];
        if (nrm2 > 1e-9) {
          double inv_n = 1.0 / std::sqrt(nrm2);
          for (int t = 0; t < nt; ++t) Q.push_back(v[t] * inv_n);
          nu_sub++;
        }
      }
    }
    if (n_act < nt && nu_sub < nt - 1) {
      Rcpp::warning(
        "Active mask voxel count (%d) is smaller than 4D residual volume count (%d); "
        "empirical residual-subspace rank is capped at %d. Consider supplying 4D `copes` "
        "or `fwe_method = \"grf_sim\"` for very small ROIs.",
        n_act, nt, nu_sub
      );
    } else if (R_FINITE(dof_hint) && dof_hint > 0.0 && std::fabs(static_cast<double>(nu_sub) - dof_hint) > 2.5) {
      Rcpp::warning(
        "Empirical residual-subspace rank (%d) differs from supplied `dof` (%.1f); "
        "using empirical residual rank for sign-flip studentization.",
        nu_sub, dof_hint
      );
    }
  }

  std::vector<double>().swap(rows_tmp); // free temporary row buffer

  double dof_t = static_cast<double>(nu_sub - 1);
  if (dof_t < 4.0) {
    stop("Too few residual degrees of freedom (eff_df < 4) for Rademacher sign-flipping calibration.");
  }

  // Build 512-node cubic Lagrange LUT for t -> Z on [0, 24] + exact log-space fallback above 24
  const int N_T_LUT = 512;
  const double T_LUT_MAX = 24.0;
  const double dt_lut = T_LUT_MAX / static_cast<double>(N_T_LUT - 1);
  const double inv_dt_lut = 1.0 / dt_lut;
  std::vector<double> t2z_lut(N_T_LUT, 0.0);
  for (int g = 1; g < N_T_LUT; ++g) {
    double tv = g * dt_lut;
    double lp = R::pt(tv, dof_t, 0, 1);
    t2z_lut[g] = R::qnorm5(lp, 0.0, 1.0, 0, 1);
  }

  auto fast_r_to_z_pos = [&](double r_pos) -> double {
    if (r_pos <= 0.0) return 0.0;
    double r_clamped = std::min(0.99999999, r_pos);
    double om_r2 = std::max(1e-15, (1.0 - r_clamped) * (1.0 + r_clamped));
    double tv = std::sqrt(dof_t / om_r2) * r_clamped;
    if (tv >= T_LUT_MAX - 2.0 * dt_lut) {
      double lp = R::pt(tv, dof_t, 0, 1);
      return R::qnorm5(lp, 0.0, 1.0, 0, 1);
    }
    return interp_cubic_equispaced(t2z_lut, 0.0, inv_dt_lut, tv);
  };

  auto exact_r_to_z_signed = [&](double r_val) -> double {
    if (r_val == 0.0) return 0.0;
    double sgn = (r_val > 0.0) ? 1.0 : -1.0;
    return sgn * fast_r_to_z_pos(std::fabs(r_val));
  };

  // Cutoff on |r| corresponding to ZestThr (with 0.99999 safety margin)
  double lp_zest = R::pnorm5(ZestThr, 0.0, 1.0, 0, 1);
  double t_zest = R::qt(lp_zest, dof_t, 0, 1);
  double r_zest = 0.99999 * (t_zest / std::sqrt(dof_t + t_zest * t_zest));

  std::vector<double> z_vol(engine.n, 0.0);
  std::vector<double> r_vec(n_act, 0.0);
  std::vector<double> s_raw(nt, 1.0);
  std::vector<double> q_flip(nt, 0.0);
  std::vector<double> q_coef(nu_sub, 0.0);
  std::vector<int> cand_pos, cand_neg, touched;
  cand_pos.reserve(n_act / 6 + 32);
  cand_neg.reserve(n_act / 6 + 32);
  touched.reserve(n_act / 3 + 64);

  NumericVector max_dist(n_perm);
  const double inv_sqrt_nt = 1.0 / std::sqrt(static_cast<double>(nt));

  for (int b = 0; b < n_perm; ++b) {
    if (b % 32 == 0) Rcpp::checkUserInterrupt();

    // Draw valid non-degenerate Rademacher vector and normalize in subspace if needed
    for (int attempt = 0; attempt < 64; ++attempt) {
      int n_pos = 0;
      for (int t = 0; t < nt; ++t) {
        if (R::unif_rand() < 0.5) {
          s_raw[t] = 1.0;
          n_pos++;
        } else {
          s_raw[t] = -1.0;
        }
      }
      if (!center_residuals) {
        // Exclude unflipped all-+1 (and all--1 in two-sided mode) if nt is large enough
        if (n_pos == nt || (sided_mode == 2 && n_pos == 0)) {
          if (attempt < 63) continue;
        }
        for (int t = 0; t < nt; ++t) q_flip[t] = s_raw[t] * inv_sqrt_nt;
        break;
      } else {
        // Project s_raw onto orthonormal residual subspace Q: q_flip = Q (Q^T s_raw)
        for (int j = 0; j < nu_sub; ++j) {
          const double *qj = &Q[static_cast<size_t>(j) * nt];
          double dot = 0.0;
          for (int t = 0; t < nt; ++t) dot += qj[t] * s_raw[t];
          q_coef[j] = dot;
        }
        std::fill(q_flip.begin(), q_flip.end(), 0.0);
        for (int j = 0; j < nu_sub; ++j) {
          const double *qj = &Q[static_cast<size_t>(j) * nt];
          double cj = q_coef[j];
          for (int t = 0; t < nt; ++t) q_flip[t] += cj * qj[t];
        }
        double nrm2 = 0.0;
        for (int t = 0; t < nt; ++t) nrm2 += q_flip[t] * q_flip[t];
        if (nrm2 >= 0.15 * static_cast<double>(nu_sub) || attempt == 63) {
          double inv_n = 1.0 / std::sqrt(std::max(1e-30, nrm2));
          for (int t = 0; t < nt; ++t) q_flip[t] *= inv_n;
          break;
        }
      }
    }

    // Contiguous SIMD DAXPY: r_vec = U_col * q_flip
    std::fill(r_vec.begin(), r_vec.end(), 0.0);
    for (int t = 0; t < nt; ++t) {
      const double qt = q_flip[t];
      const double *col = &U_col[static_cast<size_t>(t) * n_act];
      for (int a = 0; a < n_act; ++a) {
        r_vec[a] += qt * col[a];
      }
    }

    if (use_legacy) {
      for (int a = 0; a < n_act; ++a) {
        z_vol[active_vox[a]] = exact_r_to_z_signed(r_vec[a]);
      }
      max_dist[b] = engine.eval_max(z_vol);
    } else {
      cand_pos.clear();
      cand_neg.clear();
      touched.clear();
      double max_r_pos = -2.0;
      double max_r_neg = -2.0;

      for (int a = 0; a < n_act; ++a) {
        double rv = r_vec[a];
        if (sided_mode != 1) {
          if (rv > max_r_pos) max_r_pos = rv;
          if (rv > r_zest) {
            int idx = active_vox[a];
            z_vol[idx] = fast_r_to_z_pos(rv);
            cand_pos.push_back(idx);
            touched.push_back(idx);
          }
        }
        if (sided_mode != 0) {
          double nrv = -rv;
          if (nrv > max_r_neg) max_r_neg = nrv;
          if (nrv > r_zest) {
            int idx = active_vox[a];
            z_vol[idx] = -fast_r_to_z_pos(nrv);
            cand_neg.push_back(idx);
            touched.push_back(idx);
          }
        }
      }

      double sub_pos = (sided_mode != 1 && cand_pos.empty()) ? exact_r_to_z_signed(max_r_pos) : 0.0;
      double sub_neg = (sided_mode != 0 && cand_neg.empty()) ? exact_r_to_z_signed(max_r_neg) : 0.0;

      max_dist[b] = engine.eval_max(z_vol, &cand_pos, &cand_neg, sub_pos, sub_neg);

      for (int idx : touched) {
        z_vol[idx] = 0.0;
      }
    }
  }

  return List::create(
    Named("null_max_Z") = max_dist,
    Named("eff_df") = dof_t,
    Named("mean_temp_ar1") = mean_temp_ar1,
    Named("n_active") = n_act
  );
}

