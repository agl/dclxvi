/*
 * File:   dclxvi-20110718/twistpoint_fp2.c
 * Author: Ruben Niederhagen, Peter Schwabe
 * Public Domain
 */

#include "twistpoint_fp2.h"

void twistpoint_fp2_set(twistpoint_fp2_t rop, const twistpoint_fp2_t op)
{
  fp2e_set(rop->m_x, op->m_x);
  fp2e_set(rop->m_y, op->m_y);
  fp2e_set(rop->m_z, op->m_z);
  fp2e_setzero(rop->m_t);
}

void twistpoint_fp2_setneutral(twistpoint_fp2_t rop)
{
  fp2e_setone(rop->m_x);
  fp2e_setone(rop->m_y);
  fp2e_setzero(rop->m_z);
  fp2e_setzero(rop->m_t);
}

void twistpoint_fp2_neg(twistpoint_fp2_t rop, const twistpoint_fp2_t op)
{
  fp2e_set(rop->m_x, op->m_x);
  fp2e_neg(rop->m_y, op->m_y);
  fp2e_set(rop->m_z, op->m_z);
  fp2e_setzero(rop->m_t);
}

void twistpoint_fp2_set_fp2e(twistpoint_fp2_t rop, const fp2e_t x, const fp2e_t y, const fp2e_t z)
{
  fp2e_set(rop->m_x, x);
  fp2e_set(rop->m_y, y);
  fp2e_set(rop->m_z, z);
  fp2e_setzero(rop->m_t);
}

void twistpoint_fp2_affineset_fp2e(twistpoint_fp2_t rop, const fp2e_t x, const fp2e_t y)
{
  fp2e_set(rop->m_x, x);
  fp2e_set(rop->m_y, y);
  fp2e_setone(rop->m_z);
  fp2e_setzero(rop->m_t);
}

void twistpoint_fp2_mixadd(twistpoint_fp2_t rop, const twistpoint_fp2_t op1, const twistpoint_fp2_t op2)
{
  fp2e_t tfp2e1, tfp2e2, tfp2e3, tfp2e4, tfp2e5, tfp2e6, tfp2e7, tfp2e8, tfp2e9;	// Temporary variables needed for intermediary results
  fp2e_square(tfp2e1, op1->m_z);
  fp2e_mul(tfp2e2, op1->m_z, tfp2e1);
  fp2e_mul(tfp2e3, op2->m_x, tfp2e1);
  fp2e_mul(tfp2e4, op2->m_y, tfp2e2);
  fp2e_sub(tfp2e5, tfp2e3, op1->m_x);
  fp2e_short_coeffred(tfp2e5);
  fp2e_sub(tfp2e6, tfp2e4, op1->m_y);
  fp2e_short_coeffred(tfp2e6);
  fp2e_square(tfp2e7, tfp2e5);
  fp2e_mul(tfp2e8, tfp2e7, tfp2e5);
  fp2e_mul(tfp2e9, op1->m_x, tfp2e7);

  fp2e_double(tfp2e1, tfp2e9);
  fp2e_add(tfp2e1, tfp2e1, tfp2e8);
  fp2e_square(rop->m_x, tfp2e6);
  fp2e_sub(rop->m_x, rop->m_x, tfp2e1);
  fp2e_short_coeffred(rop->m_x);
  fp2e_sub(tfp2e1, tfp2e9, rop->m_x);
  fp2e_short_coeffred(tfp2e1);
  fp2e_mul(tfp2e2, tfp2e1, tfp2e6);
  fp2e_mul(tfp2e3, op1->m_y, tfp2e8);
  fp2e_sub(rop->m_y, tfp2e2, tfp2e3);
  fp2e_short_coeffred(rop->m_y);
  fp2e_mul(rop->m_z, op1->m_z, tfp2e5);
}

void twistpoint_fp2_double(twistpoint_fp2_t rop, const twistpoint_fp2_t op)
{
  fp2e_t tfp2e1, tfp2e2, tfp2e3, tfp2e4;	// Temporary variables needed for intermediary results
  fp2e_square(tfp2e1, op->m_y);
  fp2e_mul(tfp2e2, tfp2e1, op->m_x);
  fp2e_double(tfp2e2, tfp2e2);
  fp2e_double(tfp2e2, tfp2e2);
  fp2e_square(tfp2e3, tfp2e1);
  fp2e_double(tfp2e3, tfp2e3);
  fp2e_double(tfp2e3, tfp2e3);
  fp2e_double(tfp2e3, tfp2e3);
  fp2e_square(tfp2e4, op->m_x);
  fp2e_triple(tfp2e4, tfp2e4);
  fp2e_short_coeffred(tfp2e4);
  fp2e_square(rop->m_x, tfp2e4);
  fp2e_double(tfp2e1, tfp2e2);
  fp2e_sub(rop->m_x, rop->m_x, tfp2e1);
  fp2e_short_coeffred(rop->m_x);
  fp2e_sub(tfp2e1, tfp2e2, rop->m_x);
  fp2e_short_coeffred(tfp2e1);
  fp2e_mul(rop->m_z, op->m_y, op->m_z);
  fp2e_double(rop->m_z, rop->m_z);
  fp2e_short_coeffred(rop->m_z);
  fp2e_mul(rop->m_y, tfp2e4, tfp2e1);
  fp2e_sub(rop->m_y, rop->m_y, tfp2e3);
  fp2e_short_coeffred(rop->m_y);
}

void twistpoint_fp2_mul(twistpoint_fp2_t rop, const twistpoint_fp2_t op, const scalar_t scalar, const unsigned int scalar_bitsize)
{
  size_t i;
  twistpoint_fp2_t r;
  twistpoint_fp2_set(r, op);
  for (i = scalar_bitsize - 1; i > 0; i--) {
    twistpoint_fp2_double(r, r);
    if (scalar_getbit(scalar, i - 1))
      twistpoint_fp2_mixadd(r, r, op);
  }
  twistpoint_fp2_set(rop, r);
}

void twistpoint_fp2_makeaffine(twistpoint_fp2_t op)
{
  fp2e_invert(op->m_z, op->m_z);
  fp2e_mul(op->m_y, op->m_y, op->m_z);
  fp2e_square(op->m_z, op->m_z);
  fp2e_mul(op->m_x, op->m_x, op->m_z);
  fp2e_mul(op->m_y, op->m_y, op->m_z);
  fp2e_setone(op->m_z);
}

void twistpoint_fp2_print(FILE * outfile, const twistpoint_fp2_t op)
{
  fprintf(outfile, "[");
  fp2e_print(outfile, op->m_x);
  fprintf(outfile, ", ");
  fp2e_print(outfile, op->m_y);
  fprintf(outfile, ", ");
  fp2e_print(outfile, op->m_z);
  fprintf(outfile, "]");
}
