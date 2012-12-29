CPP=g++
CPPFLAGS=-g -Wall -DCHECK

CC=gcc
CFLAGS=-std=c99 -O3 -fomit-frame-pointer
LFLAGS=-lm

all: as check c

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
