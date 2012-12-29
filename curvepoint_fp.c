/*
 * File:   dclxvi-20110718/curvepoint_fp.c
 * Author: Ruben Niederhagen, Peter Schwabe
 * Public Domain
 */

#include <stdio.h>
#include <stdlib.h>
#include "fpe.h"
#include "curvepoint_fp.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//            Point initialization and deletion functions
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Global dummies usable by all curvepoints:
fpe_t curvepoint_dummy_fpe1;

// Set the coordinates of a curvepoint_fp_t by copying the coordinates from another curvepoint_fp
void curvepoint_fp_set(curvepoint_fp_t rop, const curvepoint_fp_t op)
{
  fpe_set(rop->m_x, op->m_x);
  fpe_set(rop->m_y, op->m_y);
  fpe_set(rop->m_z, op->m_z);
  fpe_setzero(rop->m_t);
}

void curvepoint_fp_setneutral(curvepoint_fp_t rop)
{
  fpe_setone(rop->m_x);
  fpe_setone(rop->m_y);
  fpe_setzero(rop->m_z);
  fpe_setzero(rop->m_t);
}

// Addition of two points, op2 is assumed to be in affine coordinates 
// For the algorithm see e.g. DA Peter Schwabe
void curvepoint_fp_mixadd(curvepoint_fp_t rop, const curvepoint_fp_t op1, const curvepoint_fp_t op2)
{
  fpe_t tfpe1, tfpe2, tfpe3, tfpe4, tfpe5, tfpe6, tfpe7, tfpe8, tfpe9;	// Temporary variables needed for intermediary results
  fpe_square(tfpe1, op1->m_z);
  fpe_mul(tfpe2, op1->m_z, tfpe1);
  fpe_mul(tfpe3, op2->m_x, tfpe1);
  fpe_mul(tfpe4, op2->m_y, tfpe2);
  fpe_sub(tfpe5, tfpe3, op1->m_x);
  fpe_short_coeffred(tfpe5);
  fpe_sub(tfpe6, tfpe4, op1->m_y);
  fpe_square(tfpe7, tfpe5);
  fpe_mul(tfpe8, tfpe7, tfpe5);
  fpe_mul(tfpe9, op1->m_x, tfpe7);

  fpe_double(tfpe1, tfpe9);
  fpe_add(tfpe1, tfpe1, tfpe8);
  fpe_square(rop->m_x, tfpe6);
  fpe_sub(rop->m_x, rop->m_x, tfpe1);
  fpe_short_coeffred(rop->m_x);
  fpe_sub(tfpe1, tfpe9, rop->m_x);
  fpe_mul(tfpe2, tfpe1, tfpe6);
  fpe_mul(tfpe3, op1->m_y, tfpe8);
  fpe_sub(rop->m_y, tfpe2, tfpe3);
  fpe_short_coeffred(rop->m_y);
  fpe_mul(rop->m_z, op1->m_z, tfpe5);
}

void curvepoint_fp_double(curvepoint_fp_t rop, const curvepoint_fp_t op)
{
  fpe_t tfpe1, tfpe2, tfpe3, tfpe4;	// Temporary variables needed for intermediary results
  fpe_square(tfpe1, op->m_y);
  fpe_mul(tfpe2, tfpe1, op->m_x);
  fpe_double(tfpe2, tfpe2);
  fpe_double(tfpe2, tfpe2);
  fpe_square(tfpe3, tfpe1);
  fpe_double(tfpe3, tfpe3);
  fpe_double(tfpe3, tfpe3);
  fpe_double(tfpe3, tfpe3);
  fpe_square(tfpe4, op->m_x);
  fpe_triple(tfpe4, tfpe4);
  fpe_short_coeffred(tfpe4);
  fpe_square(rop->m_x, tfpe4);
  fpe_double(tfpe1, tfpe2);
  fpe_sub(rop->m_x, rop->m_x, tfpe1);
  fpe_short_coeffred(rop->m_x);
  fpe_sub(tfpe1, tfpe2, rop->m_x);
  fpe_short_coeffred(tfpe1);
  fpe_mul(rop->m_z, op->m_y, op->m_z);
  fpe_double(rop->m_z, rop->m_z);
  fpe_mul(rop->m_y, tfpe4, tfpe1);
  fpe_sub(rop->m_y, rop->m_y, tfpe3);
  fpe_short_coeffred(rop->m_y);
}

void curvepoint_fp_mul(curvepoint_fp_t rop, const curvepoint_fp_t op, const scalar_t scalar, const unsigned int scalar_bitsize)
{
  size_t i;
  curvepoint_fp_t r;
  curvepoint_fp_set(r, op);
  for (i = scalar_bitsize - 1; i > 0; i--) {
    curvepoint_fp_double(r, r);
    if (scalar_getbit(scalar, i - 1))
      curvepoint_fp_mixadd(r, r, op);
  }
  curvepoint_fp_set(rop, r);
}

// Negate a point, store in rop:
void curvepoint_fp_neg(curvepoint_fp_t rop, const curvepoint_fp_t op)
{
  fpe_neg(curvepoint_dummy_fpe1, op->m_y);
  fpe_set(rop->m_x, op->m_x);
  fpe_set(rop->m_y, curvepoint_dummy_fpe1);
  fpe_set(rop->m_z, op->m_z);
}

// Transform to Affine Coordinates (z=1)
void curvepoint_fp_makeaffine(curvepoint_fp_t point)
{
  if (fpe_isone(point->m_z))
    return;

  fpe_invert(curvepoint_dummy_fpe1, point->m_z);
  fpe_mul(point->m_x, point->m_x, curvepoint_dummy_fpe1);
  fpe_mul(point->m_x, point->m_x, curvepoint_dummy_fpe1);

  fpe_mul(point->m_y, point->m_y, curvepoint_dummy_fpe1);
  fpe_mul(point->m_y, point->m_y, curvepoint_dummy_fpe1);
  fpe_mul(point->m_y, point->m_y, curvepoint_dummy_fpe1);

  fpe_setone(point->m_z);
}

// Print a point:
void curvepoint_fp_print(FILE * outfile, const curvepoint_fp_t point)
{
  fprintf(outfile, "[");
  fpe_print(outfile, point->m_x);
  fprintf(outfile, ", ");
  fpe_print(outfile, point->m_y);
  fprintf(outfile, ", ");
  fpe_print(outfile, point->m_z);
  fprintf(outfile, "]");
}
