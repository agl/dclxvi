/*
 * File:   dclxvi-20110718/speedtest.c
 * Author: Ruben Niederhagen, Peter Schwabe
 * Public Domain
 */

#include <stdio.h>
#include <stdlib.h>

#include "curvepoint_fp.h"
#include "twistpoint_fp2.h"
#include "fp12e.h"
#include "optate.h"
#include "linefunction.h"

#define NTESTS 100

#define REP 50

#ifdef __x86_64__
#define mycpucycles(RES) \
  __asm__ volatile("rdtsc;shlq $32,%%rdx;orq %%rdx,%%rax" : "=a" (RES) ::  "%rdx");
#else
#define mycpucycles(RES) \
  __asm__ volatile(".byte 15;.byte 49" : "=A" (RES));
#endif

extern const curvepoint_fp_t bn_curvegen;
extern const twistpoint_fp2_t bn_twistgen;
extern const scalar_t bn_n;

static int cmp_llu(const void *a, const void *b)
{
  if (*(unsigned long long *)a < *(unsigned long long *)b)
    return -1;
  if (*(unsigned long long *)a > *(unsigned long long *)b)
    return 1;
  return 0;
}

static unsigned long long median(unsigned long long *l, size_t llen)
{
  qsort(l, llen, sizeof(unsigned long long), cmp_llu);

  if (llen % 2)
    return l[llen / 2];
  else
    return (l[llen / 2 - 1] + l[llen / 2]) / 2;
}

static void print_bench(unsigned long long *l, size_t llen)
{
  int i;
  for (i = 0; i < (int)llen - 1; i++) {
    l[i] = l[i + 1] - l[i];
  }
  printf("Median of averages: %llu\n", median(l, llen - 1) / REP);
  printf("\n");
}

fp12e_t e1;
curvepoint_fp_t p1;
twistpoint_fp2_t p2;
twistpoint_fp2_t p3;
twistpoint_fp2_t rq;

fp2e_t rop11, rop12, rop13, r2;
fpe_t fpe1;

scalar_t s1, s2, s3;
unsigned int s1_size, s2_size, s3_size;

unsigned long long t[NTESTS];

int main(int argc, char *argv[])
{
  int i, j;

  scalar_setrandom(s1, bn_n);
  s1_size = scalar_scanb(s1) + 1;
  scalar_setrandom(s2, bn_n);
  s2_size = scalar_scanb(s2) + 1;
  scalar_setrandom(s3, bn_n);
  s3_size = scalar_scanb(s3) + 1;

  curvepoint_fp_mul(p1, bn_curvegen, s1, s1_size);
  curvepoint_fp_makeaffine(p1);
  twistpoint_fp2_mul(p2, bn_twistgen, s2, s2_size);
  twistpoint_fp2_makeaffine(p2);
  twistpoint_fp2_mul(p3, bn_twistgen, s3, s3_size);
  twistpoint_fp2_makeaffine(p3);
  fp2e_setone(rop11);
  fp2e_setone(rop12);
  fp2e_setone(rop13);
  fp2e_setone(r2);
  fpe_setone(fpe1);
  fp12e_setone(e1);

  printf("Fp2 multiplication:\n");
  for (i = 0; i < NTESTS; i++) {
    mycpucycles(t[i]);
    for (j = 0; j < REP; j++)
      fp2e_mul(r2, rop11, rop12);
  }
  print_bench(t, NTESTS);

  printf("Fp2 squaring:\n");
  for (i = 0; i < NTESTS; i++) {
    mycpucycles(t[i]);
    for (j = 0; j < REP; j++)
      fp2e_square(r2, rop11);
  }
  print_bench(t, NTESTS);

  printf("Fp2xFp multiplication:\n");
  for (i = 0; i < NTESTS; i++) {
    mycpucycles(t[i]);
    for (j = 0; j < REP; j++)
      fp2e_mul_fpe(r2, rop11, fpe1);
  }
  print_bench(t, NTESTS);

  printf("Fp2 short coeffred:\n");
  for (i = 0; i < NTESTS; i++) {
    mycpucycles(t[i]);
    for (j = 0; j < REP; j++)
      fp2e_short_coeffred(r2);
  }
  print_bench(t, NTESTS);

  printf("linefunction evaluation (addition):\n");
  for (i = 0; i < NTESTS; i++) {
    mycpucycles(t[i]);
    for (j = 0; j < REP; j++)
      linefunction_add_ate(rop11, rop12, rop13, rq, p2, p3, p1, r2);
  }
  print_bench(t, NTESTS);

  printf("linefunction evaluation (doubling):\n");
  for (i = 0; i < NTESTS; i++) {
    mycpucycles(t[i]);
    for (j = 0; j < REP; j++)
      linefunction_double_ate(rop11, rop12, rop13, rq, p2, p1);
  }
  print_bench(t, NTESTS);

  printf("Fp12 multiplication:\n");
  for (i = 0; i < NTESTS; i++) {
    mycpucycles(t[i]);
    for (j = 0; j < REP; j++)
      fp12e_mul(e1, e1, e1);
  }
  print_bench(t, NTESTS);

  printf("Fp12 squaring:\n");
  for (i = 0; i < NTESTS; i++) {
    mycpucycles(t[i]);
    for (j = 0; j < REP; j++)
      fp12e_square(e1, e1);
  }
  print_bench(t, NTESTS);

  printf("Fp12 linefunction multiplication:\n");
  for (i = 0; i < NTESTS; i++) {
    mycpucycles(t[i]);
    for (j = 0; j < REP; j++)
      fp12e_mul_line(e1, e1, rop11, rop12, rop13);
  }
  print_bench(t, NTESTS);

  printf("Fp12 inversion:\n");
  for (i = 0; i < NTESTS; i++) {
    mycpucycles(t[i]);
    for (j = 0; j < REP; j++)
      fp12e_invert(e1, e1);
  }
  print_bench(t, NTESTS);

  printf("Miller loop:\n");
  for (i = 0; i < NTESTS; i++) {
    mycpucycles(t[i]);
    for (j = 0; j < REP; j++)
      optate_miller(e1, p2, p1);
  }
  print_bench(t, NTESTS);

  printf("Optimal ate pairing:\n");
  for (i = 0; i < NTESTS; i++) {
    mycpucycles(t[i]);
    for (j = 0; j < REP; j++)
      optate(e1, p2, p1);
  }
  print_bench(t, NTESTS);
}
