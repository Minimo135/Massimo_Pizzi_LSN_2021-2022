#ifndef _GL_
#define _GL_

#include <fstream>
#include <iostream>
#include <string>
#include "random.h"
#include "mpi.h"

extern Random rnd;

struct city{
  double x;
  double y;
  int label;
  std::string Name;
  std::string State;
};

void initialize_rnd(Random &, int);
int mod(int n, int m);

void Select_pairs(int *to_nodes, int* from_nodes, int rank, int size);

#endif