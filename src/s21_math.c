#include "s21_math.h"

long double fact(double x);

long double s21_atan(double x) {
  long double res = 0;
  if (x != x)
    res = s21_nan;
  else if (x == s21_POS_INF)
    res = s21_PI / 2;
  else if (x == s21_NEG_INF)
    res = -s21_PI / 2;
  else {
    if (s21_fabs(x) < 1) {
      double tmp = x / s21_sqrt(1 + x * x);
      res = s21_asin(tmp);
    } else if (x > 1) {
      res = s21_PI / 2 - s21_atan(1 / x);
    } else {
      res = -s21_PI / 2 - s21_atan(1 / x);
    }
  }
  return res;
}

long double s21_fmod(double x, double y) {
  long double res = 0;
  if (y == 0 || x == s21_NEG_INF || x == s21_POS_INF || x != x) {
    res = s21_nan;
  } else {
    int minus = 1;
    if ((x < 0 && y > 0) || (x < 0 && y < 0)) {
      minus = -1;
    }
    x = s21_fabs(x);
    y = s21_fabs(y);
    while (x > y) {
      x -= y;
    }
    res = x * minus;
  }
  return res;
}

long double s21_pow(double base, double exp) {
  long double res = 0;
  if (base == 1 || exp == 0)
    res = 1;
  else if ((exp != exp) || (base != base))
    res = s21_nan;
  else if (base == 0) {
    if (exp == s21_NEG_INF)
      res = s21_POS_INF;
    else
      res = 0;
  } else if (base == s21_POS_INF) {
    if (exp == s21_NEG_INF || exp < 0)
      res = 0;
    else
      res = s21_POS_INF;
  } else if (base == s21_NEG_INF) {
    if (exp == s21_NEG_INF || exp < 0)
      res = 0;
    else if (exp == s21_POS_INF)
      res = s21_POS_INF;
    else {
      int tmp = (int)exp;
      if ((tmp - exp) == 0 && tmp % 2 != 0)
        res = s21_NEG_INF;
      else
        res = s21_POS_INF;
    }
  } else {
    if (((s21_fabs(base) > 1) && exp == s21_POS_INF) ||
        ((s21_fabs(base) < 1) && exp == s21_NEG_INF))
      res = s21_POS_INF;
    else if (((s21_fabs(base) < 1) && exp == s21_POS_INF) ||
             ((s21_fabs(base) > 1) && exp == s21_NEG_INF))
      res = 0;
    else
      res = s21_exp(exp * s21_log(base));
  }
  return res;
}

long double s21_sqrt(double x) {
  long double res = 0;
  if (x == s21_POS_INF)
    res = s21_POS_INF;
  else if (x != x)
    res = s21_nan;
  else if (x >= 0)
    res = s21_pow(x, 0.5);
  else
    res = -s21_nan;
  return res;
}

long double s21_tan(double x) {
  long double res = 0;
  if (x == s21_POS_INF || x == s21_NEG_INF)
    res = -s21_nan;
  else if (x != x)
    res = s21_nan;
  else
    res = s21_sin(x) / s21_cos(x);
  return res;
}

long double s21_acos(double x) {
  long double res = 0;
  if (x != x || x > 1 || x < -1)
    res = s21_nan;
  else
    res = s21_PI / 2 - s21_asin(x);
  return res;
}

long double s21_asin(double x) {
  long double res = 0;
  if (x != x || x > 1 || x < -1)
    res = s21_nan;
  else {
    if (x == 1)
      res = s21_PI / 2;
    else if (x == -1)
      res = -s21_PI / 2;
    else {
      long double numerator = 1;
      long double denumerator = 2;
      long double elem = 1;
      long double poly = x * x * x;
      int current_degree = 3;
      while (s21_fabs(elem) >= epsilon) {
        elem = numerator * poly / denumerator / current_degree;
        res += elem;
        numerator *= current_degree;
        denumerator *= current_degree + 1;
        poly *= x * x;
        current_degree += 2;
      }
      res += x;
    }
  }
  return res;
}

long double s21_sin(double x) {
  long double res = 0;
  if (x != x)
    res = s21_nan;
  else if (x == s21_POS_INF || x == s21_NEG_INF)
    res = -s21_nan;
  else {
    x = s21_fmod(x, 100000 * s21_PI);
    long double numerator = x;
    long double denumerator = 1;
    long double elem = 1;
    int n = 1;
    while (s21_fabs(elem) >= epsilon) {
      elem = numerator / denumerator;
      res += elem;
      numerator *= -x * x;
      denumerator *= (n + 1) * (n + 2);
      n += 2;
    }
  }
  return res;
}

long double s21_cos(double x) {
  long double res = 0;
  if (x != x)
    res = s21_nan;
  else if (x == s21_POS_INF || x == s21_NEG_INF)
    res = -s21_nan;
  else {
    x = s21_fmod(x, 100000 * s21_PI);
    long double element = 0;
    long double stepen = x * x;
    long double factorial = 2;
    int count = 1;
    if (s21_fabs(x) <= 2 * s21_PI) {
      do {
        element = stepen / fact(factorial);
        if (count % 2 != 0) {
          res -= element;
        } else {
          res += element;
        }
        stepen *= x * x;
        factorial += 2;
        count++;
      } while (element >= epsilon);
      res += 1;
    } else if (x > 2 * s21_PI)
      res = s21_cos(x - 2 * s21_PI);
    else
      res = s21_cos(-x - 2 * s21_PI);
  }
  return res;
}

long double s21_ceil(double x) {
  long double result = 0;
  if (x != x)
    result = s21_nan;
  else if (x == s21_POS_INF)
    result = s21_POS_INF;
  else if (x == s21_NEG_INF)
    result = s21_NEG_INF;
  else {
    if (x > 0) {
      for (int i = 0; i < x; i++) result = i;
      result = result + 1;
    } else if (x < 0)
      for (int i = 0; i > x; i--) result = i;
    else
      result = 0;
  }
  return result;
}

long double s21_floor(double x) {
  long double result = 0;
  if (x != x)
    result = s21_nan;
  else if (x == s21_POS_INF)
    result = s21_POS_INF;
  else if (x == s21_NEG_INF)
    result = s21_NEG_INF;
  else {
    if (x > 0)
      for (int i = 0; i < x; i++) result = i;
    else if (x < 0) {
      for (int i = 0; i > x; i--) result = i;
      result = result - 1;
    } else
      result = 0;
  }
  return result;
}

long double s21_fabs(double x) {
  long double res = 0;
  if (x != x)
    res = s21_nan;
  else if (x == s21_POS_INF || x == s21_NEG_INF)
    res = s21_POS_INF;
  else if (x < 0)
    res = -x;
  else
    res = x;
  return res;
}

int s21_abs(int x) {
  if (x < 0) x = -x;
  return x;
}

long double fact(double x) {
  long double result = 1;
  for (int i = x; i > 0; i--) {
    result *= i;
  }
  return result;
}

long double s21_exp(double x) {
  long double res = 0;
  if (x != x)
    res = s21_nan;
  else if (x == s21_POS_INF)
    res = s21_POS_INF;
  else if (x == s21_NEG_INF)
    res = 0;
  else {
    long double element = 0;
    long double numerator = x;
    long double factorial = 1;
    long double denominator = 1;
    do {
      element = numerator / factorial;
      res += element;
      numerator *= x;
      denominator++;
      factorial *= denominator;
    } while (element >= epsilon);
    res += 1;
  }
  return res;
}

long double s21_log(double x) {
  long double res = 0;
  if (x != x)
    res = s21_nan;
  else if (x < 0)
    res = -s21_nan;
  else if (x == s21_POS_INF)
    res = s21_POS_INF;
  else if (x == 0)
    res = s21_NEG_INF;
  else {
    long double temp = 1;
    long double xx = (x - 1) / (x + 1);
    long double stepen = xx;
    long double sum = 0;
    long double element = 0;
    do {
      element = stepen / temp;
      sum += element;
      temp += 2;
      stepen *= xx * xx;
    } while (s21_fabs(element) >= epsilon);
    res = 2 * sum;
  }
  return res;
}