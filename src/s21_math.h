#ifndef S21_MATH_H
#define S21_MATH_H
#define s21_e 2.7182818284590452353602874713527
#define s21_PI 3.1415926535897932384626433832795
#define epsilon 1e-20
#define s21_nan 0. / 0.
#define s21_POS_INF 1. / 0.
#define s21_NEG_INF -1. / 0.
long double s21_pow(double base, double exp);
/* This function returns value of BASE raised to the power of EXP.
Example:
    s21_pow(23.2, 2.5) == 2592.507146;
    s21_pow(2.235, 4.1) == 27.041941;
    s21_pow(-346, 0) == 1.000000; */
long double s21_fmod(double x, double y);
/* This function computes the floating-point remainder of dividing X by Y.
Example:
    s21_fmod(-7.2342, 92) == -7.234200;
    s21_fmod(953, -2.42352) == 0.556640;
    s21_fmod(953, 0) == -nan; */
long double s21_sqrt(double x);
/* This function returns the nonnegative square root of X.
Example:
    s21_sqrt(144) == 12.000000;
    s21_sqrt(242.2352) 15.563907;
    s21_sqrt(-78) == -nan; */
long double s21_fabs(double x);
/* This function returns the absolute value of the floating-point number X.
Example:
    s21_fabs(-14.25) == 14.250000;
    s21_fabs(234.2) == 234.200000;
    s21_fabs(-8234.325) == 8234.325000; */
long double s21_asin(double x);
/* This function calculates the principal value of the arc sine of X; that is
the value whose sine is X. Example: s21_asin(0.3324) == 0.338847;
    s21_asin(0.1242) == 0.124522;
    s21_asin(-0.666) == -0.728834; */
long double s21_acos(double x);
/* This function calculates the arc cosine of X; that is the value whose cosine
is X. Example: s21_acos(-0.5) == 2.094395; s21_acos(0.848) == 0.558596;
    s21_acos(0.1) == 1.470629; */
long double s21_floor(double x);
/* This function returns the largest integral value that is not greater than X.
Example:
    s21_floor(65.5) == 65.000000;
    s21_floor(838.999934) == 838.000000;
    s21_floor(-9.000012) == -10.000000; */
long double s21_ceil(double x);
/* This function returns the smallest integral value that is not less than X.
Example:
    s21_ceil(65.5) == 66.000000;
    s21_ceil(838.999934) == 839.000000;
    s21_ceil(-9.000012) == -9.000000; */
long double s21_atan(double x);
/* This function calculate the principal value of the arc tangent of X, where X
is given in radians. Example: s21_atan(0.562362) == 0.512285; s21_atan(0.2513)
== 0.246202; s21_atan(-0.923) == -0.745378; */
long double s21_tan(double x);
/* This function returns the tangent of X, where X is given in radians.
Example:
    s21_tan(3.5) == 0.374586;
    s21_tan(-6) == 0.291006;
    s21_tan(0.92369825) == 1.323389; */
long double s21_sin(double x);
/* This function returns the sine of X, where X is given in radians.
Example:
    s21_sin(2) == 0.909297;
    s21_sin(-4.44444444) == 0.964317;
    s21_sin(0.1234567890987654321) == 0.123143; */
long double s21_cos(double x);
/* This function returns the cosine of X, where X is given in radians.
Example:
    s21_cos(0.10239745) == 0.994762;
    s21_cos(-500.235235) == -0.750483;
    s21_cos(82347125) == 0.999076; */
long double s21_exp(double x);
/* This function returns the value of e (the base of natural logarithms) raised
to the power of X. Example: s21_exp(8.888) == 7244.515653;
    s21_exp(-0.0000098464) == 0.999990;
    s21_exp(0) == 1.000000; */
long double s21_log(double x);
/* This function returns the natural logarithm of X.
Example:
    s21_log(0.5) == -0.693147;
    s21_log(55.66666) == 4.019381;
    s21_log(-7.9234) == -nan; */
int s21_abs(int x);
/* This function returns the absolute value of the integer number X.
Example:
    s21_abs(-14) == 14;
    s21_abs(234) == 234;
    s21_abs(-74) == 74; */
#endif