#include <cassert>
#include <cmath>
#include <cassert>
#include <algorithm>

using namespace std;

#define check(x, y) \
  do { \
    assert((x) >= 0.0); \
    assert((x) <= 1.0); \
    assert((y) >= 0.0); \
    assert((y) <= 1.0); \
    if ((x) == 0.0) return 0.0; \
    if ((y) == 0.0) return 0.0; \
    if ((x) == 1.0) return (y); \
    if ((y) == 1.0) return (x); \
  } while (0)

//
// non-parameterized t-norms
//
double logical_product(double x, double y) {
  check(x, y);
  return min(x, y);
}

double hamacher_product(double x, double y) {
  check(x, y);
  return x * y / (x + y - x * y);
}

double algebraic_product(double x, double y) {
  check(x, y);
  return x * y;
}

double einstein_product(double x, double y) {
  check(x, y);
  return x * y / (1.0 + (1.0 - x) * (1.0 - y));
}

double bounded_product(double x, double y) {
  check(x, y);
  return max(0.0, x + y - 1.0);
}

double drastic_product(double x, double y) {
  check(x, y);
  return 0.0;
}

//
// parameterized t-norms
//
double yager_tnorm(double x, double y, double p) {
  assert(p > 0.0);
  check(x, y);
  return 1.0 - (min(1.0, pow(pow(1.0 - x, p) + pow(1.0 - y, p), 1.0 / p)));
}

double schweizer1_tnorm(double x, double y, double p) {
  assert(p > 0.0);
  check(x, y);
  return pow(max(0.0, pow(x, p) + pow(y, p) - 1.0), 1.0 / p);
}

double schweizer2_tnorm(double x, double y, double p) {
  assert(p > 0.0);
  check(x, y);
  return 1.0 / pow(1.0 / pow(x, p) + 1.0 / pow(y, p) - 1.0, 1.0 / p);
}

double schweizer3_tnorm(double x, double y, double p) {
  assert(p > 0.0);
  check(x, y);

  double tmp1 = pow(1.0 - x, p);
  double tmp2 = pow(1.0 - y, p);
  return 1.0 - pow(tmp1 + tmp2 - tmp1 * tmp2, 1.0 / p);
}

double hamacher_tnorm(double x, double y, double p) {
  assert(p >= 0.0);
  check(x, y);
  return x * y / (p + (1.0 - p) * (x + y - x * y));
}

double frank_tnorm(double x, double y, double p) {
  assert((p > 0.0) && (p != 1.0));
  check(x, y);
  return log(1.0 + (pow(p, x) - 1.0) * (pow(p, y) - 1.0) / (p - 1.0)) / log(p);
}

double dombi_tnorm(double x, double y, double p) {
  assert(p > 0.0);
  check(x, y);
  return 1.0 / (1.0 + pow(pow((1.0 - x) / x, p) + pow((1.0 - y) / y, p), 1.0 / p));
}

double weber_tnorm(double x, double y, double p) {
  assert(p >= -1.0);
  check(x, y);
  return max((1.0 + p) * (x + y - 1.0) - p * x * y, 0.0);
}

double dubois_tnorm(double x, double y, double p) {
  assert((p >= 0.0) && (p <= 1.0));
  check(x, y);
  return x * y / max(max(x, y), p);
}
