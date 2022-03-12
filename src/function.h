#ifndef _FUNCTION_H

#include "data_file.h"

class Function {
private:
    //Pointeur de la classe DataFile
    DataFile* _df;

    double _Lx, _Ly;

public:
  // Constructeur
  Function(DataFile* data_file);

  double f(double _x, double _y, double _t);
  double g(double _x, double _y);
  double h(double _x, double _y);
  double sol_exacte(double _x, double _y);
};

#define _FUNCTION_H
#endif
