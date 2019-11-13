#ifndef TNORM_H_
#define TNORM_H_

//
// reference:
//   - Masaharu Mizumoto,
//     Pictorial representations of fuzzy connectives, Part I: Cases of t-norms, t-conorms and averaging operators,
//     Fuzzy Sets and Systems, vol.31, no.2, pp.217-242, 1989
//

// t-norms
extern double logical_product(double x, double y);
extern double hamacher_product(double x, double y);
extern double algebraic_product(double x, double y);
extern double einstein_product(double x, double y);
extern double bounded_product(double x, double y);
extern double drastic_product(double x, double y);

// parameterized t-norms
extern double yager_tnorm(double x, double y, double p);
extern double schweizer1_tnorm(double x, double y, double p);
extern double schweizer2_tnorm(double x, double y, double p);
extern double schweizer3_tnorm(double x, double y, double p);
extern double hamacher_tnorm(double x, double y, double p);
extern double frank_tnorm(double x, double y, double p);
extern double dombi_tnorm(double x, double y, double p);
extern double weber_tnorm(double x, double y, double p);
extern double dubois_tnorm(double x, double y, double p);

#endif // TNORM_H_
