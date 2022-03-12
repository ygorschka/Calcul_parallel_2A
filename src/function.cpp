#ifndef _FUNCTION_CPP

#include "function.h"
#include <math.h>

Function::Function(DataFile* data_file) :
    _df(data_file)
{
    _Ly = _df->Get_Ly();
    _Lx = _df->Get_Lx();
}

// Condition initiale
double Function::f(double _x, double _y, double _t)
{
    if (_df->Get_scenario() == "cas1") {
        return 2*(_y-pow(_y,2)+_x-pow(_x,2));
    }
    else if (_df->Get_scenario() == "cas2") {
        return sin(_x) + cos(_y);
    }
    else if (_df->Get_scenario() == "cas3")
    {
        return exp(-pow((_x-_Lx/2),2))*exp(-pow((_y-_Ly/2),2))*cos((M_PI/2.)*_t);
    }
    //Cas par défaut cas1
    else {
        return 2*(_y-pow(_y,2)+_x-pow(_x,2));
    }
}

double Function::g(double _x, double _y)
{
    if (_df->Get_scenario() == "cas1") {
        return 0;
    }
    else if (_df->Get_scenario() == "cas2") {
        return sin(_x) + cos(_y);
    }
    else if (_df->Get_scenario() == "cas3")
    {
        return 0;
    }
    //Cas par défaut cas1
    else {
        return 0;
    }
}

double Function::h(double _x, double _y)
{
    if (_df->Get_scenario() == "cas1") {
        return 0;
    }
    else if (_df->Get_scenario() == "cas2") {
        return sin(_x) + cos(_y);
    }
    else if (_df->Get_scenario() == "cas3")
    {
        return 1;
    }
    //Cas par défaut cas1
    else {
        return 0;
    }
}

double Function::sol_exacte(double _x, double _y)
{
    if (_df->Get_scenario() == "cas1") {
        return _x*(1.0-_x)*_y*(1.0-_y);
    }
    else if (_df->Get_scenario() == "cas2") {
        return sin(_x) + cos(_y);
    }
    else {
        return 0.0;
    }
}

#define _FUNCTION_CPP
#endif
