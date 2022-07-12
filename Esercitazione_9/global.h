#ifndef _GL_
#define _GL_

#include <fstream>
#include <iostream>
#include <string>
#include "random.h"

extern Random rnd;

struct city{
  double x;
  double y;
  int label;
};

void initialize_rnd(Random &);
int mod(int n, int m);

#endif