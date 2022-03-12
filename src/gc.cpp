#ifndef _GC_CPP

#include "gc.h"
#include <math.h>
#include <cstdio>
#include <mpi.h>

// Constructeur
GC::GC(Function* function, DataFile* data_file, Tensor* tensor) :
_df(data_file), _tensor(tensor), _fct(function)
{
    _Nx = _df->Get_Nx();
    _Ny = _df->Get_Ny();
    _size_vect = _df->Get_size_vect();
    _rank = _df->Get_rank();
    _F.resize(_size_vect);
    _A.resize(3);
}

std::vector<double> GC::gradient(double epsilon, int kmax, double t, std::vector<double> U)
{
    _tensor->Build_F(t);
    _tensor->Build_A();
    _F = _tensor->Get_F();
    _A = _tensor->Get_A();

    int n, ierr;
    n = _size_vect;
    std::vector<double> r(n,0.), p(n,0.), r_next(n,0.), z(n,0.);

    double beta_p, beta, sum_beta;
    int k;

    std::vector<double> x_aux(n, 1.), x_aux2(n, 1.), p_aux(n, 0.), x_aux3(n, 1.);

    x_aux2 = _tensor->Prd_Mat_Vec(_A,U);
    x_aux = produit_cst(-1,x_aux2);
    x_aux3 = somme_vec(_F, U);
    r = somme_vec(x_aux3, x_aux);

    p = r;

    //Calcul du beta pour le vecteur entier
    sum_beta = 0.0;
    beta_p = produit_sc(p,p);
    ierr = MPI_Allreduce(&beta_p, &sum_beta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    beta = sqrt(sum_beta);

    k = 0;

    while ((beta > epsilon)&&(k< kmax))
    {
        double alpha_p_top, alpha_p_bottom, alpha, alpha_top_sum, alpha_bottom_sum;
        double gamma_p_top, gamma, gamma_top_sum, gamma_bottom_sum;

        z = _tensor->Prd_Mat_Vec(_A,p);

        //Calcul de alpha pour le vecteur entier
        alpha_top_sum = 0.0;
        alpha_bottom_sum = 0.0;
        alpha_p_top = produit_sc(r,r);
        alpha_p_bottom = produit_sc(z,p);
        ierr = MPI_Allreduce(&alpha_p_top, &alpha_top_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        ierr = MPI_Allreduce(&alpha_p_bottom, &alpha_bottom_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        alpha = alpha_top_sum/alpha_bottom_sum;

        p_aux = produit_cst(alpha,p);
        U = somme_vec(U,p_aux);
        p_aux = produit_cst(-alpha,z);
        r_next = somme_vec(r,p_aux);

        //Calcul de gamma pour le vecteur entier
        gamma_top_sum = 0.0;
        gamma_p_top = produit_sc(r_next,r_next);
        ierr = MPI_Allreduce(&gamma_p_top, &gamma_top_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        gamma_bottom_sum = alpha_top_sum;
        gamma = gamma_top_sum/gamma_bottom_sum;

        p_aux = produit_cst(gamma,p);
        p = somme_vec(r_next,p_aux);
        //Calcul de beta pour le vecteur entier
        sum_beta = 0.0;
        beta_p = produit_sc(r_next,r_next);
        ierr = MPI_Allreduce(&beta_p, &sum_beta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        beta = sqrt(sum_beta);

        k += 1;
        r = r_next;
    };

    if (k>kmax) {
        printf("Tolérance %f non atteinte \n",beta );
    }

    return U;
};

double GC::produit_sc(std::vector<double> a, std::vector<double> b)
{int n;

    double somme = 0.0;
    n=a.size();

    for (int i = 0; i < n; i++) {
        somme=somme + a[i]*b[i];
    }

    return somme;
};

std::vector<double> GC::produit_cst(double a, std::vector<double> b)
{
    int n;
    n=b.size();

    for (int i = 0; i < n; i++) {
        b[i]=a*b[i];
    }

    return b;
};

std::vector<double> GC::somme_vec(std::vector<double> a, std::vector<double> b)
{
    int n;
    double somme;
    n=a.size();

    for (int i = 0; i < n; i++) {
        a[i]=a[i]+b[i];
    }

    return a;
};

#define _GC_CPP
#endif
