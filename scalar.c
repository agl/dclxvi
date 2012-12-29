/*
 * File:   dclxvi-20110718/scalar.c
 * Author: Ruben Niederhagen, Peter Schwabe
 * Public Domain
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "scalar.h"

void scalar_setrandom(scalar_t rop, const scalar_t bound)
{
  int i;
  FILE *urand = fopen("/dev/urandom", "r");
  if (urand == NULL)
  {
    fprintf(stderr, "Could not open device file /dev/urandom");
    exit(1);
  }
  do
  {
    for(i=0;i<32;i++) 
      i[(unsigned char*)rop] = fgetc(urand);
  }
  while(!scalar_le(rop,bound));
  fclose(urand);
}

void scalar_set_lluarray(scalar_t rop, unsigned long long v[4])
{
  int i;
  for(i=0;i<4;i++) rop[i] = v[i];
}

int scalar_getbit(const scalar_t s, unsigned int pos)
{
  assert(pos < 256);
  return (s[pos >> 6] >> (pos & 0x3f)) & 1;
}

// Returns the position of the most significant set bit
int scalar_scanb(const scalar_t s)
{
  int i;
  unsigned int pos = 0;
  for(i=255;i>0;i--)
    if(scalar_getbit(s,i) && pos == 0) pos = i;
  return pos;
}

// Returns 1 if a < b, 0 otherwise
int scalar_le(const scalar_t a, const scalar_t b)
{
  int r;
  r = a[3] < b[3];
  if(a[3] == b[3])
  {
    r = a[2] < b[2];
    if(a[2] == b[2])
    {
      r = a[1] < b[1];
      if(a[0] == b[0])
        r = a[0] < b[0];
    }
  }
  return r;
}

void scalar_print(FILE *fh, const scalar_t t)
{
  int i;
  for(i=3;i>=0;i--)
    fprintf(fh, "%llx", t[i]);
}

