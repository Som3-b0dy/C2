#include <check.h>
#include <stdio.h>
#include <stdlib.h>

#include "unit_test.h"

int main(void) {
  int failed = 0;
  Suite *s_1, *s_2, *s_3, *s_4, *s_5, *s_6, *s_7, *s_8, *s_9, *s_10, *s_11,
      *s_12, *s_13, *s_14, *s_15;
  SRunner *runner;

  s_1 = suite_pow();
  s_2 = suite_fmod();
  s_3 = suite_sqrt();
  s_4 = suite_fabs();
  s_5 = suite_asin();
  s_6 = suite_acos();
  s_7 = suite_floor();
  s_8 = suite_ceil();
  s_9 = suite_atan();
  s_10 = suite_tan();
  s_11 = suite_sin();
  s_12 = suite_cos();
  s_13 = suite_exp();
  s_14 = suite_log();
  s_15 = suite_abs();

  runner = srunner_create(s_1);
  srunner_add_suite(runner, s_2);
  srunner_add_suite(runner, s_3);
  srunner_add_suite(runner, s_4);
  srunner_add_suite(runner, s_5);
  srunner_add_suite(runner, s_6);
  srunner_add_suite(runner, s_7);
  srunner_add_suite(runner, s_8);
  srunner_add_suite(runner, s_9);
  srunner_add_suite(runner, s_10);
  srunner_add_suite(runner, s_11);
  srunner_add_suite(runner, s_12);
  srunner_add_suite(runner, s_13);
  srunner_add_suite(runner, s_14);
  srunner_add_suite(runner, s_15);

  srunner_run_all(runner, CK_NORMAL);

  failed = srunner_ntests_failed(runner);
  srunner_free(runner);
  return (failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}