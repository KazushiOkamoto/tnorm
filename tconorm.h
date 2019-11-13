#ifndef TCONORM_H_
#define TCONORM_H_

//
// reference:
//   - Masaharu Mizumoto,
//     Pictorial representations of fuzzy connectives, Part I: Cases of t-norms, t-conorms and averaging operators,
//     Fuzzy Sets and Systems, vol.31, no.2, pp.217-242, 1989
//

// t-conorms
extern double logical_sum(double x, double y);
extern double hamacher_sum(double x, double y);
extern double algebraic_sum(double x, double y);
extern double einstein_sum(double x, double y);
extern double bounded_sum(double x, double y);
extern double drastic_sum(double x, double y);

// parameterized t-conorms
extern double yager_tconorm(double x, double y, double p);
extern double schweizer1_tconorm(double x, double y, double p);
extern double schweizer2_tconorm(double x, double y, double p);
extern double schweizer3_tconorm(double x, double y, double p);
extern double hamacher_tconorm(double x, double y, double p);
extern double frank_tconorm(double x, double y, double p);
extern double dombi_tconorm(double x, double y, double p);
extern double weber_tconorm(double x, double y, double p);
extern double dubois_tconorm(double x, double y, double p);

#endif // TCONORM_H_
