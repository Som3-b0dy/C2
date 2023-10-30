SYSSHELL:=$(shell uname)
CC=gcc
CFLAGS=-Wall -Wextra -Werror -std=c11
ifeq ($(SYSSHELL), Linux)
TFLAGS=-lm -lcheck -lpthread -lrt -lsubunit --coverage 
else
TFLAGS=$(shell pkg-config --cflags --libs check) --coverage
endif
SRC1=test_s21_math.c
SRC2=test_main.c
SRC3=s21_math.c

.PHONY: all clean

all: s21_math.a test gcov_report

s21_math.a:
	$(CC) $(CFLAGS) $(SRC3) -c
	ar rc s21_math.a s21_math.o

test: s21_math.a
	$(CC) $(CFLAGS) $(SRC1) $(SRC2) $(TFLAGS) -c
	$(CC) $(CFLAGS) $(SRC1) $(SRC2) $(SRC3) $(TFLAGS) -o test -L. s21_math.a
	./test

gcov_report: test
	lcov -t "s21_math" -o s21_math.info -c -d . 
	genhtml -o report s21_math.info

clean:
	rm -rf *.gcno *.info *.gcda *.o *.a test report