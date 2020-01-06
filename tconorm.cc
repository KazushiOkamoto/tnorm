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
    if ((x) == 0.0) return (y); \
    if ((y) == 0.0) return (x); \
    if ((x) == 1.0) return 1.0; \
    if ((y) == 1.0) return 1.0; \
  } while (0)

//
// non-parameterized t-conorms
//
double logical_sum(double x, double y) {
  check(x, y);
  return max(x, y);
}

double hamacher_sum(double x, double y) {
  check(x, y);
  return (x + y - 2.0 * x * y) / (1.0 - x * y);
}

double algebraic_sum(double x, double y) {
  check(x, y);
  return x + y - x * y;
}

double einstein_sum(double x, double y) {
  check(x, y);
  return (x + y) / (1.0 + x * y);
}

double bounded_sum(double x, double y) {
  check(x, y);
  return max(1.0, x + y);
}

double drastic_sum(double x, double y) {
  check(x, y);
  return 1.0;
}

//
// parameterized t-conorms
//
double yager_tconorm(double x, double y, double p) {
  assert(p > 0.0);
  check(x, y);
  return min(1.0, pow(pow(x, p) + pow(y, p), 1.0 / p));
}

double schweizer1_tconorm(double x, double y, double p) {
  assert(p > 0.0);
  check(x, y);
  return 1.0 - pow(max(0.0, pow(1.0 - x, p) + pow(1.0 - y, p) - 1.0), 1.0 / p);
}

double schweizer2_tconorm(double x, double y, double p) {
  assert(p > 0.0);
  check(x, y);
  return 1.0 - 1.0 / pow(1.0 / pow(1.0 - x, p) + 1.0 / pow(1.0 - y, p) - 1.0, 1.0 / p);
}

double schweizer3_tconorm(double x, double y, double p) {
  assert(p > 0.0);
  check(x, y);

  double tmp1 = pow(x, p);
  double tmp2 = pow(y, p);
  return pow(tmp1 + tmp2 - tmp1 * tmp2, 1.0 / p);
}

double hamacher_tconorm(double x, double y, double p) {
  assert(p >= 0.0);
  check(x, y);
  return (x + y - x * y - (1.0 - p) * x * y) / (1.0 - (1.0 - p) * x * y);
}

double frank_tconorm(double x, double y, double p) {
  assert((p > 0.0) && (p != 1.0));
  check(x, y);
  return 1.0 - log(1.0 + (pow(p, 1.0 - x) - 1.0) * (pow(p, 1.0 - y) - 1.0) / (p - 1.0)) / log(p);
}

double dombi_tconorm(double x, double y, double p) {
  assert(p > 0.0);
  check(x, y);
  return 1.0 - 1.0 / (1.0 + pow(pow(x/ (1.0 - x), p) + pow(y / (1.0 - y), p), 1.0 / p));
}

double weber_tconorm(double x, double y, double p) {
  assert(p >= -1.0);
  check(x, y);
  return min(x + y + p * x * y, 1.0);
}

double dubois_tconorm(double x, double y, double p) {
  assert((p >= 0.0) && (p <= 1.0));
  check(x, y);
  return 1.0 - (1.0 - x) * (1.0 - y) / max(max(1.0 - x, 1.0 - y), p);
}
