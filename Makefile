CPP=g++
CPPFLAGS=-g -Wall -DCHECK

CC=gcc
CFLAGS=-std=c99 -O3 -fomit-frame-pointer
LFLAGS=-lm

all: as check c libdclxvipairing.so

libdclxvipairing.so:
	gcc -shared -Wl,-soname=libdclxvipairing.so -o libdclxvipairing.so -O3 -fomit-frame-pointer -fPIC -DQHASM linefunction.c optate.c fpe.c fp2e.c fp6e.c fp12e.c curvepoint_fp.c twistpoint_fp2.c final_expo.c scalar.c parameters.c mul.c mydouble.c fp2e_add2.s fp2e_sub2.s fp2e_double2.s fp2e_triple2.s fp2e_neg2.s fp2e_mul.s fp2e_mul_fpe.s fp2e_short_coeffred.s fp2e_add.s fp2e_sub.s fp2e_parallel_coeffmul.s fp2e_mulxi.s fp2e_double.s fp2e_triple.s fp2e_neg.s fp2e_conjugate.s fpe_mul.s fp2e_square.s consts.s

c: bilintest-c \
	 speedtest-c

as: bilintest-as \
		speedtest-as

check: bilintest-check \
	     speedtest-check

bilintest-check: bilintest.c linefunction.c optate.c fpe.c fp2e.c fp6e.c fp12e.c curvepoint_fp.c twistpoint_fp2.c final_expo.c scalar.c parameters.c mul.c mydouble.c
	$(CPP) $(CPPFLAGS) -DNTESTS=20 -o $@ $^

bilintest-c: bilintest.c linefunction.c optate.c fpe.c fp2e.c fp6e.c fp12e.c curvepoint_fp.c twistpoint_fp2.c final_expo.c scalar.c parameters.c mul.c mydouble.c
	$(CC) $(CFLAGS) $(LFLAGS) -DNTESTS=1000 -o $@ $^

bilintest-as: bilintest.c linefunction.c optate.c fpe.c fp2e.c fp6e.c fp12e.c curvepoint_fp.c twistpoint_fp2.c final_expo.c scalar.c parameters.c mul.c mydouble.c asfunctions.a
	$(CC) $(CFLAGS) $(LFLAGS) -DQHASM -DNTESTS=1000000 -o $@ $^

speedtest-check: speedtest.c linefunction.c optate.c fpe.c fp2e.c fp6e.c fp12e.c curvepoint_fp.c twistpoint_fp2.c final_expo.c scalar.c parameters.c mul.c mydouble.c
	$(CPP) $(CPPFLAGS) -o $@ $^

speedtest-c: speedtest.c linefunction.c optate.c fpe.c fp2e.c fp6e.c fp12e.c curvepoint_fp.c twistpoint_fp2.c final_expo.c scalar.c parameters.c mul.c mydouble.c
	$(CC) $(CFLAGS) $(LFLAGS) -o $@ $^

speedtest-as: speedtest.c linefunction.c optate.c fpe.c fp2e.c fp6e.c fp12e.c curvepoint_fp.c twistpoint_fp2.c final_expo.c scalar.c parameters.c mul.c mydouble.c asfunctions.a
	$(CC) $(CFLAGS) $(LFLAGS) -DQHASM -o $@ $^

%.o: %.s
	$(CC) $(CFLAGS) -c -o $@ $^

asfunctions.a: fp2e_add2.o fp2e_sub2.o \
	fp2e_double2.o fp2e_triple2.o fp2e_neg2.o \
	fp2e_mul.o fp2e_mul_fpe.o fp2e_short_coeffred.o \
	fp2e_add.o fp2e_sub.o fp2e_parallel_coeffmul.o fp2e_mulxi.o\
	fp2e_double.o fp2e_triple.o fp2e_neg.o fp2e_conjugate.o \
	fpe_mul.o fp2e_square.o \
	consts.o
	rm -f asfunctions.a
	ar cr asfunctions.a $^

.PHONY: clean

clean:
	-rm bilintest-check 
	-rm speedtest-check 
	-rm bilintest-c
	-rm speedtest-c
	-rm bilintest-as 
	-rm speedtest-as 
	-rm *.o
	-rm asfunctions.a
	-rm libdclxvipairing.so
