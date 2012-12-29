/*
 * File:   dclxvi-20110718/bilintest.c
 * Author: Ruben Niederhagen, Peter Schwabe
 * Public Domain
 */

#include <stdio.h>

#include "curvepoint_fp.h"
#include "twistpoint_fp2.h"
#include "fp12e.h"
#include "optate.h"

extern const curvepoint_fp_t bn_curvegen;
extern const twistpoint_fp2_t bn_twistgen;
extern const scalar_t bn_n;

int main(int argc, char* argv[])
{
  fp12e_t e1, e2, e3;

  curvepoint_fp_t p1;
  twistpoint_fp2_t p2;

  scalar_t s1, s2;
  unsigned int s1_size, s2_size;

  int i;

  // Test with neutral element as argument
  scalar_setrandom(s1, bn_n);
  s1_size = scalar_scanb(s1)+1;
  curvepoint_fp_set(p1, bn_curvegen);
  twistpoint_fp2_setneutral(p2);
  fpe_isreduced(p1->m_x);
  fpe_isreduced(p1->m_y);

  curvepoint_fp_mul(p1, p1, s1, s1_size);
  curvepoint_fp_makeaffine(p1);

  optate(e1, p2, p1);

  if(!fp12e_isone(e1))
    printf("Error in optimal ate: e(infty,P) != 1\n");

  scalar_setrandom(s2, bn_n);
  s2_size = scalar_scanb(s2)+1;
  curvepoint_fp_setneutral(p1);
  twistpoint_fp2_set(p2, bn_twistgen);
  fp2e_isreduced(p2->m_x);
  fp2e_isreduced(p2->m_y);

  twistpoint_fp2_mul(p2, p2, s2, s2_size);
  twistpoint_fp2_makeaffine(p2);

  optate(e1, p2, p1);
  
  if(!fp12e_isone(e1))
    printf("Error in optimal ate: e(Q,infty) != 1\n");


  // Bilinearity test of optimal ate Pairing:
  for(i=0;i<NTESTS;i++)
  {
#if (NTESTS > 100)
    if(!(i%(NTESTS/100)) && i!=0) printf("Number of tests: %d\n",i);
#else
    if(i!=0) printf("Number of tests: %d\n",i);
#endif
    scalar_setrandom(s1, bn_n);
    scalar_setrandom(s2, bn_n);
    s1_size = scalar_scanb(s1)+1;
    s2_size = scalar_scanb(s2)+1;
    curvepoint_fp_set(p1, bn_curvegen);
    twistpoint_fp2_set(p2, bn_twistgen);
    fpe_isreduced(p1->m_x);
    fpe_isreduced(p1->m_y);
    fp2e_isreduced(p2->m_x);
    fp2e_isreduced(p2->m_y);

    curvepoint_fp_mul(p1, p1, s1, s1_size);
    curvepoint_fp_makeaffine(p1);
    twistpoint_fp2_mul(p2, p2, s2, s2_size);
    twistpoint_fp2_makeaffine(p2);

    optate(e1, p2, p1);
    
    curvepoint_fp_set(p1, bn_curvegen);
    twistpoint_fp2_set(p2, bn_twistgen);
    fpe_isreduced(p1->m_x);
    fpe_isreduced(p1->m_y);
    fp2e_isreduced(p2->m_x);
    fp2e_isreduced(p2->m_y);
    curvepoint_fp_mul(p1, p1, s2, s2_size);
    curvepoint_fp_makeaffine(p1);
    twistpoint_fp2_mul(p2, p2, s1, s1_size);
    twistpoint_fp2_makeaffine(p2);

    optate(e2, p2, p1);
    
    curvepoint_fp_set(p1, bn_curvegen);
    twistpoint_fp2_set(p2, bn_twistgen);

    optate(e3, p2, p1);

    fp12e_pow(e3, e3, s1, s1_size);
    fp12e_pow(e3, e3, s2, s2_size);

    if(!fp12e_iseq(e1,e2))
    {
      printf("Error in optimal ate: e1 != e2\n");
      printf("e1: ");
      fp12e_print(stdout, e1);
      printf("\ne2: ");
      fp12e_print(stdout, e2);
      printf("\nScalars:\n");
      printf("s1: ");
      scalar_print(stdout, s1); 
      printf("\ns2: ");
      scalar_print(stdout, s2); 
      printf("\n");
    }
    else if(!fp12e_iseq(e2,e3))
    {
      printf("Error in optimal ate: e2 != e3\n");
      printf("e2: ");
      fp12e_print(stdout, e2);
      printf("\ne3: ");
      fp12e_print(stdout, e3);
      printf("\nScalars:\n");
      printf("s1: ");
      scalar_print(stdout, s1); 
      printf("\ns2: ");
      scalar_print(stdout, s2); 
      printf("\n");
    }
    else if(fp12e_iszero(e2))
      printf("Error: Pairing value is zero\n");
    else if(fp12e_isone(e2))
      printf("Warning: Pairing value is one\n");
  }
  return 0;
}
