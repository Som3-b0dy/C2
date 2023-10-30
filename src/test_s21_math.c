#include <check.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "s21_math.h"
#include "unit_test.h"

START_TEST(test_fmod) {
  double x = -7.2342, y = 92;
  double res_1 = s21_fmod(x, y);
  double res_2 = fmod(x, y);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = 953, y = -2.42352;
  res_1 = s21_fmod(x, y);
  res_2 = fmod(x, y);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = 953, y = 0;
  res_1 = s21_fmod(x, y);
  ck_assert_ldouble_nan(res_1);

  x = s21_POS_INF, y = s21_nan;
  res_1 = s21_fmod(x, y);
  ck_assert_ldouble_nan(res_1);

  x = s21_POS_INF, y = -s21_nan;
  res_1 = s21_fmod(x, y);
  ck_assert_ldouble_nan(res_1);

  x = s21_POS_INF, y = s21_POS_INF;
  res_1 = s21_fmod(x, y);
  ck_assert_ldouble_nan(res_1);

  x = s21_POS_INF, y = s21_NEG_INF;
  res_1 = s21_fmod(x, y);
  ck_assert_ldouble_nan(res_1);

  x = s21_NEG_INF, y = s21_nan;
  res_1 = s21_fmod(x, y);
  ck_assert_ldouble_nan(res_1);

  x = s21_NEG_INF, y = -s21_nan;
  res_1 = s21_fmod(x, y);
  ck_assert_ldouble_nan(res_1);

  x = s21_NEG_INF, y = s21_POS_INF;
  res_1 = s21_fmod(x, y);
  ck_assert_ldouble_nan(res_1);

  x = s21_NEG_INF, y = s21_NEG_INF;
  res_1 = s21_fmod(x, y);
  ck_assert_ldouble_nan(res_1);

  x = s21_nan, y = s21_nan;
  res_1 = s21_fmod(x, y);
  ck_assert_ldouble_nan(res_1);

  x = s21_nan, y = -s21_nan;
  res_1 = s21_fmod(x, y);
  ck_assert_ldouble_nan(res_1);

  x = s21_nan, y = s21_POS_INF;
  res_1 = s21_fmod(x, y);
  ck_assert_ldouble_nan(res_1);

  x = s21_nan, y = s21_NEG_INF;
  res_1 = s21_fmod(x, y);
  ck_assert_ldouble_nan(res_1);

  x = -s21_nan, y = s21_nan;
  res_1 = s21_fmod(x, y);
  ck_assert_ldouble_nan(res_1);

  x = -s21_nan, y = -s21_nan;
  res_1 = s21_fmod(x, y);
  ck_assert_ldouble_nan(res_1);

  x = -s21_nan, y = s21_POS_INF;
  res_1 = s21_fmod(x, y);
  ck_assert_ldouble_nan(res_1);

  x = -s21_nan, y = s21_NEG_INF;
  res_1 = s21_fmod(x, y);
  ck_assert_ldouble_nan(res_1);

  x = 23.2, y = 823.2;
  res_1 = s21_fmod(x, y);
  res_2 = fmod(x, y);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = -52.3, y = -9.1;
  res_1 = s21_fmod(x, y);
  res_2 = fmod(x, y);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = 0.0, y = 0.0;
  res_1 = s21_fmod(x, y);
  ck_assert_ldouble_nan(res_1);

  x = -36, y = 0.0;
  res_1 = s21_fmod(x, y);
  ck_assert_ldouble_nan(res_1);

  x = 0.0, y = 24.2;
  res_1 = s21_fmod(x, y);
  res_2 = fmod(x, y);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = 0.0, y = -2.1;
  res_1 = s21_fmod(x, y);
  res_2 = fmod(x, y);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);
}
END_TEST

START_TEST(test_sqrt) {
  double x = 144;
  double res_1 = s21_sqrt(x);
  double res_2 = sqrt(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = 242.2352;
  res_1 = s21_sqrt(x);
  res_2 = sqrt(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = -78;
  res_1 = s21_sqrt(x);
  ck_assert_ldouble_nan(res_1);

  x = s21_POS_INF;
  res_1 = s21_sqrt(x);
  ck_assert_ldouble_infinite(res_1);

  x = s21_nan;
  res_1 = s21_sqrt(x);
  ck_assert_ldouble_nan(res_1);
}
END_TEST

START_TEST(test_asin) {
  double x = 0.3324;
  double res_1 = s21_asin(x);
  double res_2 = asin(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = 0.1242;
  res_1 = s21_asin(x);
  res_2 = asin(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = -0.666;
  res_1 = s21_asin(x);
  res_2 = asin(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = 1;
  res_1 = s21_asin(x);
  res_2 = asin(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = -1;
  res_1 = s21_asin(x);
  res_2 = asin(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = -5.25;
  res_1 = s21_asin(x);
  ck_assert_ldouble_nan(res_1);

  x = 7.213;
  res_1 = s21_asin(x);
  ck_assert_ldouble_nan(res_1);

  x = s21_nan;
  res_1 = s21_asin(x);
  ck_assert_ldouble_nan(res_1);
}
END_TEST

START_TEST(test_tan) {
  double x = 3.5;
  double res_1 = s21_tan(x);
  double res_2 = tan(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = -6;
  res_1 = s21_tan(x);
  res_2 = tan(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = 0.92369825;
  res_1 = s21_tan(x);
  res_2 = tan(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = s21_nan;
  res_1 = s21_tan(x);
  ck_assert_ldouble_nan(res_1);

  x = s21_POS_INF;
  res_1 = s21_tan(x);
  ck_assert_ldouble_nan(res_1);

  x = s21_NEG_INF;
  res_1 = s21_tan(x);
  ck_assert_ldouble_nan(res_1);
}
END_TEST

START_TEST(test_sin) {
  double x = 2;
  double res_1 = s21_sin(x);
  double res_2 = sin(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = -4.44444444;
  res_1 = s21_sin(x);
  res_2 = sin(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = 0.1234567890987654321;
  res_1 = s21_sin(x);
  res_2 = sin(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = s21_nan;
  res_1 = s21_sin(x);
  ck_assert_ldouble_nan(res_1);

  x = s21_POS_INF;
  res_1 = s21_sin(x);
  ck_assert_ldouble_nan(res_1);

  x = s21_NEG_INF;
  res_1 = s21_sin(x);
  ck_assert_ldouble_nan(res_1);
}
END_TEST

START_TEST(test_pow) {
  double base = 23.2, exp = 2.5;
  double res_1 = s21_pow(base, exp);
  double res_2 = pow(base, exp);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  base = 2.235, exp = 4.1;
  res_1 = s21_pow(base, exp);
  res_2 = pow(base, exp);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  base = -346, exp = 0;
  res_1 = s21_pow(base, exp);
  res_2 = pow(base, exp);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  base = s21_nan, exp = 1;
  res_1 = s21_pow(base, exp);
  ck_assert_ldouble_nan(res_1);

  base = 0, exp = s21_NEG_INF;
  res_1 = s21_pow(base, exp);
  ck_assert_ldouble_infinite(res_1);

  base = s21_POS_INF, exp = -200;
  res_1 = s21_pow(base, exp);
  res_2 = pow(base, exp);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  base = s21_POS_INF, exp = 3.12;
  res_1 = s21_pow(base, exp);
  ck_assert_ldouble_infinite(res_1);

  base = s21_NEG_INF, exp = -200;
  res_1 = s21_pow(base, exp);
  res_2 = pow(base, exp);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  base = s21_NEG_INF, exp = s21_POS_INF;
  res_1 = s21_pow(base, exp);
  ck_assert_ldouble_infinite(res_1);

  base = s21_NEG_INF, exp = 3;
  res_1 = s21_pow(base, exp);
  ck_assert_ldouble_infinite(res_1);

  base = 1.5, exp = s21_POS_INF;
  res_1 = s21_pow(base, exp);
  ck_assert_ldouble_infinite(res_1);

  base = -0.3, exp = s21_POS_INF;
  res_1 = s21_pow(base, exp);
  res_2 = pow(base, exp);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  base = 0, exp = 3.3123123;
  res_1 = s21_pow(base, exp);
  res_2 = pow(base, exp);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  base = s21_NEG_INF, exp = 4;
  res_1 = s21_pow(base, exp);
  ck_assert_ldouble_infinite(res_1);
}
END_TEST

START_TEST(test_fabs) {
  double x = -14.25;
  double res_1 = s21_fabs(x);
  double res_2 = fabs(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = 234.2;
  res_1 = s21_fabs(x);
  res_2 = fabs(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = -8234.325;
  res_1 = s21_fabs(x);
  res_2 = fabs(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = s21_nan;
  res_1 = s21_fabs(x);
  ck_assert_ldouble_nan(res_1);

  x = s21_POS_INF;
  res_1 = s21_fabs(x);
  ck_assert_ldouble_infinite(res_1);
}
END_TEST

START_TEST(test_acos) {
  double x = -0.5;
  double res_1 = s21_acos(x);
  double res_2 = acos(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = 0.848;
  res_1 = s21_acos(x);
  res_2 = acos(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = 0.1;
  res_1 = s21_acos(x);
  res_2 = acos(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = s21_nan;
  res_1 = s21_acos(x);
  ck_assert_ldouble_nan(res_1);
}
END_TEST

START_TEST(test_floor) {
  double x = 65.5;
  double res_1 = s21_floor(x);
  double res_2 = floor(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = 838.999934;
  res_1 = s21_floor(x);
  res_2 = floor(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = -9.000012;
  res_1 = s21_floor(x);
  res_2 = floor(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = s21_nan;
  res_1 = s21_floor(x);
  ck_assert_ldouble_nan(res_1);

  x = s21_POS_INF;
  res_1 = s21_floor(x);
  ck_assert_ldouble_infinite(res_1);

  x = s21_NEG_INF;
  res_1 = s21_floor(x);
  ck_assert_ldouble_infinite(res_1);

  x = 0;
  res_1 = s21_floor(x);
  res_2 = floor(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);
}
END_TEST

START_TEST(test_ceil) {
  double x = 65.5;
  double res_1 = s21_ceil(x);
  double res_2 = ceil(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = 838.999934;
  res_1 = s21_ceil(x);
  res_2 = ceil(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = -9.000012;
  res_1 = s21_ceil(x);
  res_2 = ceil(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = s21_nan;
  res_1 = s21_ceil(x);
  ck_assert_ldouble_nan(res_1);

  x = s21_POS_INF;
  res_1 = s21_ceil(x);
  ck_assert_ldouble_infinite(res_1);

  x = s21_NEG_INF;
  res_1 = s21_ceil(x);
  ck_assert_ldouble_infinite(res_1);

  x = 0;
  res_1 = s21_ceil(x);
  res_2 = ceil(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);
}
END_TEST

START_TEST(test_atan) {
  double x = 0.562362;
  double res_1 = s21_atan(x);
  double res_2 = atan(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = 0.2513;
  res_1 = s21_atan(x);
  res_2 = atan(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = -0.923;
  res_1 = s21_atan(x);
  res_2 = atan(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = s21_nan;
  res_1 = s21_atan(x);
  ck_assert_ldouble_nan(res_1);

  x = s21_POS_INF;
  res_1 = s21_atan(x);
  res_2 = atan(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = s21_NEG_INF;
  res_1 = s21_atan(x);
  res_2 = atan(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = 3.1223;
  res_1 = s21_atan(x);
  res_2 = atan(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = -3.1223;
  res_1 = s21_atan(x);
  res_2 = atan(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);
}
END_TEST

START_TEST(test_cos) {
  double x = 0.10239745;
  double res_1 = s21_cos(x);
  double res_2 = cos(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = -500.235235;
  res_1 = s21_cos(x);
  res_2 = cos(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = 2347125;
  res_1 = s21_cos(x);
  res_2 = cos(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = s21_nan;
  res_1 = s21_cos(x);
  ck_assert_ldouble_nan(res_1);

  x = s21_POS_INF;
  res_1 = s21_cos(x);
  ck_assert_ldouble_nan(res_1);
}
END_TEST

START_TEST(test_exp) {
  double x = 8.888;
  double res_1 = s21_exp(x);
  double res_2 = exp(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = -0.0000098464;
  res_1 = s21_exp(x);
  res_2 = exp(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = 0;
  res_1 = s21_exp(x);
  res_2 = exp(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = s21_nan;
  res_1 = s21_exp(x);
  ck_assert_ldouble_nan(res_1);

  x = s21_POS_INF;
  res_1 = s21_exp(x);
  ck_assert_ldouble_infinite(res_1);

  x = s21_NEG_INF;
  res_1 = s21_exp(x);
  res_2 = exp(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);
}
END_TEST

START_TEST(test_log) {
  double x = 0.5;
  double res_1 = s21_log(x);
  double res_2 = log(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = 55.66666;
  res_1 = s21_log(x);
  res_2 = log(x);
  ck_assert_ldouble_eq_tol(res_1, res_2, 0.0000001);

  x = -7.9234;
  res_1 = s21_log(x);
  ck_assert_ldouble_nan(res_1);

  x = s21_nan;
  res_1 = s21_log(x);
  ck_assert_ldouble_nan(res_1);

  x = s21_POS_INF;
  res_1 = s21_log(x);
  ck_assert_ldouble_infinite(res_1);

  x = 0;
  res_1 = s21_log(x);
  ck_assert_ldouble_infinite(res_1);
}
END_TEST

START_TEST(test_abs) {
  int x = -14;
  int res_1 = s21_abs(x);
  int res_2 = abs(x);
  ck_assert_int_eq(res_1, res_2);

  x = 234;
  res_1 = s21_abs(x);
  res_2 = abs(x);
  ck_assert_int_eq(res_1, res_2);

  x = -74;
  res_1 = s21_abs(x);
  res_2 = abs(x);
  ck_assert_int_eq(res_1, res_2);
}
END_TEST

Suite *suite_pow(void) {
  Suite *s;
  TCase *tc_pow;

  s = suite_create("s21_pow");
  tc_pow = tcase_create("s21_pow tcase");

  tcase_add_test(tc_pow, test_pow);
  suite_add_tcase(s, tc_pow);
  tcase_set_timeout(tc_pow, 100);

  return s;
}

Suite *suite_fmod(void) {
  Suite *s;
  TCase *tc_fmod;

  s = suite_create("s21_fmod");
  tc_fmod = tcase_create("s21_fmod tcase");

  tcase_add_test(tc_fmod, test_fmod);
  suite_add_tcase(s, tc_fmod);
  tcase_set_timeout(tc_fmod, 100);

  return s;
}

Suite *suite_sqrt(void) {
  Suite *s;
  TCase *tc_sqrt;

  s = suite_create("s21_sqrt");
  tc_sqrt = tcase_create("s21_sqrt tcase");

  tcase_add_test(tc_sqrt, test_sqrt);
  suite_add_tcase(s, tc_sqrt);
  tcase_set_timeout(tc_sqrt, 100);

  return s;
}

Suite *suite_fabs(void) {
  Suite *s;
  TCase *tc_fabs;

  s = suite_create("s21_fabs");
  tc_fabs = tcase_create("s21_fabs tcase");
  tcase_add_test(tc_fabs, test_fabs);
  suite_add_tcase(s, tc_fabs);
  tcase_set_timeout(tc_fabs, 100);

  return s;
}

Suite *suite_asin(void) {
  Suite *s;
  TCase *tc_asin;

  s = suite_create("s21_asin");
  tc_asin = tcase_create("s21_asin tcase");
  tcase_add_test(tc_asin, test_asin);
  suite_add_tcase(s, tc_asin);
  tcase_set_timeout(tc_asin, 100);

  return s;
}

Suite *suite_acos(void) {
  Suite *s;
  TCase *tc_acos;

  s = suite_create("s21_acos");
  tc_acos = tcase_create("s21_acos tcase");
  tcase_add_test(tc_acos, test_acos);
  suite_add_tcase(s, tc_acos);
  tcase_set_timeout(tc_acos, 100);

  return s;
}

Suite *suite_floor(void) {
  Suite *s;
  TCase *tc_floor;

  s = suite_create("s21_floor");
  tc_floor = tcase_create("s21_floor tcase");
  tcase_add_test(tc_floor, test_floor);
  suite_add_tcase(s, tc_floor);
  tcase_set_timeout(tc_floor, 100);

  return s;
}

Suite *suite_ceil(void) {
  Suite *s;
  TCase *tc_ceil;

  s = suite_create("s21_ceil");
  tc_ceil = tcase_create("s21_ceil tcase");
  tcase_add_test(tc_ceil, test_ceil);
  suite_add_tcase(s, tc_ceil);
  tcase_set_timeout(tc_ceil, 100);

  return s;
}

Suite *suite_atan(void) {
  Suite *s;
  TCase *tc_atan;

  s = suite_create("s21_atan");
  tc_atan = tcase_create("s21_atan tcase");
  tcase_add_test(tc_atan, test_atan);
  suite_add_tcase(s, tc_atan);
  tcase_set_timeout(tc_atan, 100);

  return s;
}

Suite *suite_tan(void) {
  Suite *s;
  TCase *tc_tan;

  s = suite_create("s21_tan");
  tc_tan = tcase_create("s21_tan tcase");
  tcase_add_test(tc_tan, test_tan);
  suite_add_tcase(s, tc_tan);
  tcase_set_timeout(tc_tan, 100);

  return s;
}

Suite *suite_sin(void) {
  Suite *s;
  TCase *tc_sin;

  s = suite_create("s21_sin");
  tc_sin = tcase_create("s21_sin tcase");
  tcase_add_test(tc_sin, test_sin);
  suite_add_tcase(s, tc_sin);
  tcase_set_timeout(tc_sin, 100);

  return s;
}

Suite *suite_cos(void) {
  Suite *s;
  TCase *tc_cos;

  s = suite_create("s21_cos");
  tc_cos = tcase_create("s21_cos tcase");
  tcase_add_test(tc_cos, test_cos);
  suite_add_tcase(s, tc_cos);
  tcase_set_timeout(tc_cos, 100);

  return s;
}

Suite *suite_exp(void) {
  Suite *s;
  TCase *tc_exp;

  s = suite_create("s21_exp");
  tc_exp = tcase_create("s21_exp tcase");
  tcase_add_test(tc_exp, test_exp);
  suite_add_tcase(s, tc_exp);
  tcase_set_timeout(tc_exp, 100);

  return s;
}

Suite *suite_log(void) {
  Suite *s;
  TCase *tc_log;

  s = suite_create("s21_log");
  tc_log = tcase_create("s21_log tcase");
  tcase_add_test(tc_log, test_log);
  suite_add_tcase(s, tc_log);
  tcase_set_timeout(tc_log, 100);

  return s;
}

Suite *suite_abs(void) {
  Suite *s;
  TCase *tc_abs;

  s = suite_create("s21_abs");
  tc_abs = tcase_create("s21_abs tcase");
  tcase_add_test(tc_abs, test_abs);
  suite_add_tcase(s, tc_abs);
  tcase_set_timeout(tc_abs, 100);

  return s;
}