#ifndef _DATA_FILE_CPP

#include "data_file.h"
#include <iostream>
#include <fstream>
#include <mpi.h>

using namespace std;

//Constructeur
DataFile::DataFile(std::string file_name):
    _file_name(file_name)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &_nproc);
}

void DataFile::Read_data_file()
{
    //Lecture dans le fichier
    ifstream mon_flux("data.txt");
    mon_flux >> _Nx;
    mon_flux >> _Ny;
    mon_flux >> _Lx;
    mon_flux >> _Ly;
    mon_flux >> _D;
    mon_flux >> _dt;
    mon_flux >> _Tmax;
    mon_flux >> _wich_scenario;
    mon_flux.close();

    //Définit iBeg et iEnd
    charge_c(_rank, _Ny, _nproc, &_iBeg, &_iEnd);

    //Définit la taille des vecteurs
    _size_vect = _Nx*(_iEnd-_iBeg + 1);

    //Définit le décalage dans les vecteurs
    _dec = _iBeg*_Nx;
}

void DataFile::charge_c(int me, int N, int np, int *iBeg, int *iEnd)
{
    int r = N%np;

    if (me < r){
        *iBeg = me*(N/np + 1);
        *iEnd = *iBeg + (N/np + 1) - 1;
    }
    else{
        *iBeg = r + me*(N/np);
        *iEnd = *iBeg + (N/np) - 1;
    }
}

#define _DATA_FILE_CPP
#endif
