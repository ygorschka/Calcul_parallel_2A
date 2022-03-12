#ifndef _TENSEUR_H

#include "data_file.h"
#include "function.h"
#include <vector>

using namespace std;

class Tensor{
private:
    //Pointeur de la classe DataFile
    DataFile* _df;
    //Pointeur de la classe function
    Function* _fct;

    std::vector<double> _F, _A, _U_bottom, _U_top;

    double _alpha, _beta, _gamma, _Lx, _Ly, _D, _dt, _dx, _dy;
    int _Nx, _Ny, _size_vect, _iBeg, _iEnd, _dec, _rank, _nproc;

public:

    //Constructeur
    Tensor(DataFile* data_file, Function* function);

    //Construit le second membre F
    void Build_F(const double& t);

    //Renvoie le second membre F
    const std::vector<double>& Get_F() const {return _F;};

    //Construit la matrice A sous forme de vecteur
    void Build_A();

    //Renvoie le matrice A
    const std::vector<double>& Get_A() const {return _A;};

    //Produit matrice vecteur
    std::vector<double> Prd_Mat_Vec(std::vector<double> A, std::vector<double> F);

};

#define _TENSEUR_H
#endif
