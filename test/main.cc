#include <cstdio>
#include <vector>
#include "../tnorm.h"
#include "../tconorm.h"

using namespace std;

struct Point {
  double x;
  double y;
};

const double EPS = 1e-13;

vector<Point> generate(int n) {
  vector<Point> points;

  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      Point point;
      point.x = i / double(n);
      point.y = j / double(n);
      points.push_back(point);
    }

  return points;
}

template<class T, class S>
void test(T tnorm, S tconorm, const vector<Point> &points, const char *label) {
  for (auto point: points) {
    double tmp1 = 1.0 - tnorm(1.0 - point.x, 1.0 - point.y);
    double tmp2 = tconorm(point.x, point.y);
    double diff = tmp1 - tmp2;

    if (diff > EPS)
      printf("%s: x = %f, y = %f, diff = %e\n", label, point.x, point.y, diff);
  }
}

template<class T, class S>
void test(T tnorm, S tconorm, double p, const vector<Point> &points, const char *label) {
  for (auto point: points) {
    double tmp1 = 1.0 - tnorm(1.0 - point.x, 1.0 - point.y, p);
    double tmp2 = tconorm(point.x, point.y, p);
    double diff = tmp1 - tmp2;

    if (diff > EPS)
      printf("%s: x = %f, y = %f, p = %f, diff = %e\n", label, point.x, point.y, p, diff);
  }
}

int main(int argc, char **argv) {
  int n = 1000;
  vector<Point> points = generate(n);

  // non-parameterized t-norms and t-conorms
  test(logical_product, logical_sum, points, "logical product & sum");
  test(hamacher_product, hamacher_sum, points, "hamacher product & sum");
  test(algebraic_product, algebraic_sum, points, "algebraic product & sum");
  test(einstein_product, einstein_sum, points, "einstein product & sum");
  test(bounded_product, bounded_sum, points, "bounded product & sum");
  test(drastic_product, drastic_sum, points, "drastic product & sum");

  // p = 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.5, 9.0, 10.0
  vector<double> parameters1{0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.5, 9.0, 10.0};

  for (auto p: parameters1) {
    test(yager_tnorm, yager_tconorm, p, points, "yager tnorm & tconorm");
    test(schweizer1_tnorm, schweizer1_tconorm, p, points, "schweizer1 tnorm & tconorm");
    test(schweizer2_tnorm, schweizer2_tconorm, p, points, "schweizer2 tnorm & tconorm");
    test(schweizer3_tnorm, schweizer3_tconorm, p, points, "schweizer3 tnorm & tconorm");
    test(hamacher_tnorm, hamacher_tconorm, p, points, "hamacher tnorm & tconorm");

    if (p != 1.0)
      test(frank_tnorm, frank_tconorm, p, points, "frank tnorm & tconorm");

    test(dombi_tnorm, dombi_tconorm, p, points, "dombi tnorm & tconorm");
  }

  // p = -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.5, 9.0, 10.0
  vector<double> parameters2{-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.5, 9.0, 10.0};

  for (auto p: parameters2)
    test(weber_tnorm, weber_tconorm, p, points, "weber tnorm & tconorm");

  // p = 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
  vector<double> parameters3{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

  for (auto p: parameters3)
    test(dubois_tnorm, dubois_tconorm, p, points, "dubois tnorm & tconorm");

  return 0;
}
