#ifndef _TENSEUR_CPP

#include "tensor.h"
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <mpi.h>

using namespace std;

//Constructeur
Tensor::Tensor(DataFile* data_file, Function* function) :
    _df(data_file), _fct(function)
{
    _Nx = _df->Get_Nx();
    _Ny = _df->Get_Ny();
    _Lx = _df->Get_Lx();
    _Ly = _df->Get_Ly();
    _D = _df->Get_D();
    _dt = _df->Get_dt();
    _size_vect = _df->Get_size_vect();
    _iBeg = _df->Get_iBeg();
    _iEnd = _df->Get_iEnd();
    _dec = _df->Get_dec();
    _rank = _df->Get_rank();
    _nproc = _df->Get_nproc();

    _dx = _Lx/(_Nx+1);
    _dy = _Ly/(_Ny+1);
}

//Construit le second membre F
void Tensor::Build_F(const double& t)
{
    _F.resize(_size_vect);

    for (int i=0; i<_Nx; ++i) {
        for (int j=_iBeg; j<=_iEnd; ++j) {
            if ((i == 0) && (j == 0)) {
                _F[j*_Nx + i - _dec] = _dt*(_fct->f(_dx, _dy, t)) + _D*_dt*(_fct->g(_dx,0)/pow(_dy,2) + _fct->h(0,_dy)/pow(_dx,2));
            }
            else if ((i == 0) && (j == _Ny-1)) {
                _F[j*_Nx + i - _dec] = _dt*(_fct->f(_dx, (j+1)*_dy, t)) + _D*_dt*(_fct->g(_dx,_Ly)/pow(_dy,2) + _fct->h(0,(j+1)*_dy)/pow(_dx,2));
            }
            else if ((i == _Nx-1) && (j == 0)) {
                _F[j*_Nx + i - _dec] = _dt*(_fct->f((i+1)*_dx, _dy, t)) + _D*_dt*(_fct->g((i+1)*_dx,0)/pow(_dy,2) + _fct->h(_Lx,_dy)/pow(_dx,2));
            }
            else if ((i == _Nx-1) && (j == _Ny-1)) {
                _F[j*_Nx + i - _dec] = _dt*(_fct->f((i+1)*_dx, (j+1)*_dy, t)) + _D*_dt*(_fct->g((i+1)*_dx,0)/pow(_dy,2) + _fct->h(_Lx,(j+1)*_dy)/pow(_dx,2));
            }
            else if (i == 0) {
                _F[j*_Nx + i - _dec] = _dt*(_fct->f(_dx, (j+1)*_dy, t)) + _D*_dt*(_fct->h(0,(j+1)*_dy)/pow(_dx,2));
            }
            else if (i == _Nx-1) {
                _F[j*_Nx + i - _dec] = _dt*(_fct->f((i+1)*_dx, (j+1)*_dy, t)) + _D*_dt*(_fct->h(_Lx,(j+1)*_dy)/pow(_dx,2));
            }
            else if (j == 0) {
                _F[j*_Nx + i - _dec] = _dt*(_fct->f((i+1)*_dx, _dy, t)) + _D*_dt*(_fct->g((i+1)*_dx,0)/pow(_dy,2));
            }
            else if (j == _Ny-1) {
                _F[j*_Nx + i - _dec] = _dt*(_fct->f((i+1)*_dx, (j+1)*_dy, t)) + _D*_dt*(_fct->g((i+1)*_dx,_Ly)/pow(_dy,2));
            }
            else {
                _F[j*_Nx + i - _dec] = _dt*(_fct->f((i+1)*_dx, (j+1)*_dy, t));
            }
        }
    }
}

void Tensor::Build_A()
{
    _A.resize(3);

    _alpha = 1.0 + _D*_dt*(2/pow(_dx,2) + 2/pow(_dy,2));
    _beta = -_D*_dt/pow(_dx,2);
    _gamma = -_D*_dt/pow(_dy,2);

    _A[0] = _alpha;
    _A[1] = _beta;
    _A[2] = _gamma;
}

std::vector<double> Tensor::Prd_Mat_Vec(std::vector<double> A, std::vector<double> F)
{
    int ierr, send_beg, send_end;

    std::vector<double> Res(_size_vect);

    _U_top.resize(_Nx);
    _U_bottom.resize(_Nx);

    send_beg = _iBeg*_Nx-_dec;
    send_end = _iEnd*_Nx-_dec;

    if (_rank == 0) {
        ierr = MPI_Send(&F[send_end], _Nx, MPI_DOUBLE, _rank+1, 0, MPI_COMM_WORLD);
        ierr = MPI_Recv(&_U_top[0], _Nx, MPI_DOUBLE, _rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else if (_rank == _nproc - 1) {
        ierr = MPI_Send(&F[send_beg], _Nx, MPI_DOUBLE, _rank-1, 0, MPI_COMM_WORLD);
        ierr = MPI_Recv(&_U_bottom[0], _Nx, MPI_DOUBLE, _rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else {
        ierr = MPI_Send(&F[send_beg], _Nx, MPI_DOUBLE, _rank-1, 0, MPI_COMM_WORLD);
        ierr = MPI_Send(&F[send_end], _Nx, MPI_DOUBLE, _rank+1, 0, MPI_COMM_WORLD);
        ierr = MPI_Recv(&_U_bottom[0], _Nx, MPI_DOUBLE, _rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        ierr = MPI_Recv(&_U_top[0], _Nx, MPI_DOUBLE, _rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    for(int i = 0; i<_Nx; ++i) {
        for(int j = _iBeg; j<=_iEnd; ++j) {
            if ((j == 0) && (_iBeg == _iEnd)){
                if (i == 0) {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*F[j*_Nx+i+1-_dec] + A[2]*_U_top[i];
                }
                else if (i == _Nx-1) {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*F[j*_Nx+i-1-_dec] + A[2]*_U_top[i];
                }
                else {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*(F[j*_Nx+i+1-_dec]+F[j*_Nx+i-1-_dec]) + A[2]*_U_top[i];
                }
            }
            else if (j == 0){
                if (i == 0) {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*F[j*_Nx+i+1-_dec] + A[2]*F[(j+1)*_Nx+i-_dec];
                }
                else if (i == _Nx-1) {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*F[j*_Nx+i-1-_dec] + A[2]*F[(j+1)*_Nx+i-_dec];
                }
                else {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*(F[j*_Nx+i+1-_dec]+F[j*_Nx+i-1-_dec]) + A[2]*F[(j+1)*_Nx+i-_dec];
                }
            }
            else if ((j == _Ny-1) && (_iBeg == _iEnd)) {
                if (i == 0) {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*F[j*_Nx+i+1-_dec] + A[2]*_U_bottom[i];
                }
                else if (i == _Nx-1) {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*F[j*_Nx+i-1-_dec] + A[2]*_U_bottom[i];
                }
                else {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*(F[j*_Nx+i+1-_dec]+F[j*_Nx+i-1-_dec]) + A[2]*_U_bottom[i];
                }
            }
            else if (j == _Ny-1) {
                if (i == 0) {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*F[j*_Nx+i+1-_dec] + A[2]*F[(j-1)*_Nx+i-_dec];
                }
                else if (i == _Nx-1) {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*F[j*_Nx+i-1-_dec] + A[2]*F[(j-1)*_Nx+i-_dec];
                }
                else {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*(F[j*_Nx+i+1-_dec]+F[j*_Nx+i-1-_dec]) + A[2]*F[(j-1)*_Nx+i-_dec];
                }
            }
            else if ((j == _iBeg) && (_iBeg == _iEnd)) {
                if (i == 0) {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*F[j*_Nx+i+1-_dec] + A[2]*(_U_top[i] + _U_bottom[i]);
                }
                else if (i == _Nx-1) {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*F[j*_Nx+i-1-_dec] + A[2]*(_U_top[i] + _U_bottom[i]);
                }
                else {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*(F[j*_Nx+i+1-_dec]+F[j*_Nx+i-1-_dec]) + A[2]*(_U_top[i] + _U_bottom[i]);
                }
            }
            else if (j == _iBeg) {
                if (i == 0) {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*F[j*_Nx+i+1-_dec] + A[2]*(F[(j+1)*_Nx+i-_dec] + _U_bottom[i]);
                }
                else if (i == _Nx-1) {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*F[j*_Nx+i-1-_dec] + A[2]*(F[(j+1)*_Nx+i-_dec] + _U_bottom[i]);
                }
                else {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*(F[j*_Nx+i+1-_dec]+F[j*_Nx+i-1-_dec]) + A[2]*(F[(j+1)*_Nx+i-_dec] + _U_bottom[i]);
                }
            }
            else if (j == _iEnd) {
                if (i == 0) {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*F[j*_Nx+i+1-_dec] + A[2]*(_U_top[i] + F[(j-1)*_Nx+i-_dec]);
                }
                else if (i == _Nx-1) {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*F[j*_Nx+i-1-_dec] + A[2]*(_U_top[i] + F[(j-1)*_Nx+i-_dec]);
                }
                else {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*(F[j*_Nx+i+1-_dec]+F[j*_Nx+i-1-_dec]) + A[2]*(_U_top[i] + F[(j-1)*_Nx+i-_dec]);
                }
            }
            else {
                if (i == 0) {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*F[j*_Nx+i+1-_dec] + A[2]*(F[(j+1)*_Nx+i-_dec] + F[(j-1)*_Nx+i-_dec]);
                }
                else if (i == _Nx-1) {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*F[j*_Nx+i-1-_dec] + A[2]*(F[(j+1)*_Nx+i-_dec] + F[(j-1)*_Nx+i-_dec]);
                }
                else {
                    Res[j*_Nx+i-_dec] = A[0]*F[j*_Nx+i-_dec] + A[1]*(F[j*_Nx+i+1-_dec]+F[j*_Nx+i-1-_dec]) + A[2]*(F[(j+1)*_Nx+i-_dec] + F[(j-1)*_Nx+i-_dec]);
                }
            }
        }
    }

    return Res;

}

#define _TENSEUR_CPP
#endif
