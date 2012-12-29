/*
 * File:   dclxvi-20110718/scalar.h
 * Author: Ruben Niederhagen, Peter Schwabe
 * Public Domain
 */

#ifndef SCALAR_H
#define SCALAR_H
#include <stdio.h>

typedef unsigned long long scalar_t[4] ;

void scalar_setrandom(scalar_t rop, const scalar_t bound);

void scalar_set_lluarray(scalar_t rop, unsigned long long v[4]);

int scalar_getbit(const scalar_t s, unsigned int pos);

// Returns the position of the most significant set bit
int scalar_scanb(const scalar_t s);

int scalar_le(const scalar_t a, const scalar_t b);

void scalar_print(FILE *fh, const scalar_t t);

#endif
