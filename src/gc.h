#ifndef _GC_H

#include "data_file.h"
#include "tensor.h"
#include "function.h"
#include <vector>

class GC
{
private:
    //Pointeur de la classe DataFile
    DataFile* _df;
    //Pointeur de la classe tenseur
    Tensor* _tensor;
    //Pointeur de la classe Function
    Function* _fct;

    std::vector<double> _F, _A;
    int _Nx, _Ny, _rank, _size_vect;

public:
    // Constructeur
    GC(Function* function, DataFile* data_file, Tensor* tensor);

    std::vector<double> gradient(double epsilon, int kmax, double t, std::vector<double> U);

    double produit_sc(std::vector<double> a, std::vector<double> b);
    std::vector<double>  produit_cst(double a, std::vector<double> b);
    std::vector<double> somme_vec(std::vector<double> a, std::vector<double> b);
    std::vector<double> mat_vec(std::vector<double> a, std::vector<double> b);
};

#define _GC_H
#endif
