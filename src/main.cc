// Projet realise dans le cadre du cours de 2eme annee de l'enseirb-matmeca AN202 "Calcul haute 
// performance"
// Auteurs : Joseph Bregeat - jbregeat@enseirb-marmeca.fr
//           Yoan Gorschka - ygorschka@enseirb-matmeca.fr

#include <cstdio>
#include <string>
#include <cstdlib>
#include <unistd.h>
#include <vector>
#include "gc.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <chrono>

using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
      printf("Please, enter the name of your data file.\n");
      exit(0);
    }
    const string data_file_name = argv[1];


    // Démarrage du chrono
    auto start = chrono::high_resolution_clock::now();

    MPI_Init(&argc, &argv);

    // ----------------------- Fichier de données --------------------------------
    DataFile* data_file = new DataFile(data_file_name);
    data_file->Read_data_file();

    //----------------------------- Fonction -------------------------------------
    Function* fct = new Function(data_file);

    //----------------------------- Tenseur --------------------------------------
    Tensor* tensor = new Tensor(data_file, fct);

    //-------------------------------- GC ----------------------------------------
    GC* gdt = new GC(fct, data_file, tensor);

    //Variables du problème
    double dt, Tmax, tn, dx, dy, prov, Lx, Ly;
    int nb_iter, Nx, Ny, rank, nproc, iEnd, iBeg, size_vect, dec;

    //Initialisation de certaines variables
    dt = data_file -> Get_dt();
    Tmax = data_file -> Get_Tmax();
    Nx = data_file -> Get_Nx();
    Ny = data_file -> Get_Ny();
    Lx = data_file -> Get_Lx();
    Ly = data_file -> Get_Ly();
    rank = data_file -> Get_rank();
    nproc = data_file -> Get_nproc();
    iBeg = data_file -> Get_iBeg();
    iEnd = data_file -> Get_iEnd();
    size_vect = data_file -> Get_size_vect();
    dec = data_file -> Get_dec();
    dx = Lx/(Nx+1);
    dy = Ly/(Ny+1);

    //Vecteur U
    std::vector<double> U(size_vect, 1.0);

    //Nombre d'itérations
    nb_iter = int(ceil(Tmax/dt));

    //Boucle en temps
    for (int i = 0; i <= nb_iter; ++i) {
        tn = i*dt;
        U = gdt->gradient(0.0001, 10000, tn, U);
    }

    //Pour le nom du fichier
    string name_file = data_file->Get_scenario()+'_'+"rank"+std::to_string(rank)+'_'+std::to_string(nproc)+ ".txt";

    //Sauvegarde de la solution au temps _Tmax
    ofstream mon_flux; // Contruit un objet "ofstream"
    mon_flux.open(name_file, ios::out); // Ouvre un fichier appelé name_file
    if(mon_flux) {// Vérifie que le fichier est bien ouvert
        for (int i = 0; i <= Nx+1; ++i) {
            if (rank == 0) {
                for (int j = iBeg; j <= iEnd+1; ++j) {
                    if ((i == 0) && (j == 0)) {
                        prov = fct->g(i*dx, j*dy);
                        mon_flux << i*dx << " " << j*dy << " " << prov << " " << endl;
                    }
                    else if ((i == Nx+1) && (j == 0)) {
                        prov = fct->g(i*dx, j*dy);
                        mon_flux << i*dx << " " << j*dy << " " << prov << " " << endl;
                    }
                    else if (i == Nx+1) {
                        prov = fct->h(i*dx, j*dy);
                        mon_flux << i*dx << " " << j*dy << " " << prov << " " << endl;
                    }
                    else if (i == 0) {
                        prov = fct->h(i*dx, j*dy);
                        mon_flux << i*dx << " " << j*dy << " " << prov << " " << endl;
                    }
                    else if (j == 0) {
                        prov = fct->g(i*dx, j*dy);
                        mon_flux << i*dx << " " << j*dy << " " << prov << " " << endl;
                    }
                    else {
                        mon_flux << i*dx << " " << j*dy << " " << U[(j-1)*Nx+i-1-dec] << " " << endl;
                    }
                }
            }
            else if (rank == nproc-1) {
                for (int j = iBeg+1; j <= iEnd+2; ++j) {
                    if ((i == 0) && (j == Ny+1)) {
                        prov = fct->g(i*dx, j*dy);
                        mon_flux << i*dx << " " << j*dy << " " << prov << " " << endl;
                    }
                    else if ((i == Nx+1) && (j == Ny+1)) {
                        prov = fct->g(i*dx, j*dy);
                        mon_flux << i*dx << " " << j*dy << " " << prov << " " << endl;
                    }
                    else if (i == Nx+1) {
                        prov = fct->h(i*dx, j*dy);
                        mon_flux << i*dx << " " << j*dy << " " << prov << " " << endl;
                    }
                    else if (i == 0) {
                        prov = fct->h(i*dx, j*dy);
                        mon_flux << i*dx << " " << j*dy << " " << prov << " " << endl;
                    }
                    else if (j == Ny+1) {
                        prov = fct->g(i*dx, j*dy);
                        mon_flux << i*dx << " " << j*dy << " " << prov << " " << endl;
                    }
                    else {
                        mon_flux << i*dx << " " << j*dy << " " << U[(j-1)*Nx+i-1-dec] << " " << endl;
                    }
                }
            }
            else {
                for (int j = iBeg+1; j <= iEnd+1; ++j) {
                    if (i == Nx+1) {
                        prov = fct->h(i*dx, j*dy);
                        mon_flux << i*dx << " " << j*dy << " " << prov << " " << endl;
                    }
                    else if (i == 0) {
                        prov = fct->h(i*dx, j*dy);
                        mon_flux << i*dx << " " << j*dy << " " << prov << " " << endl;
                    }
                    else {
                        mon_flux << i*dx << " " << j*dy << " " << U[(j-1)*Nx+i-1-dec] << " " << endl;
                    }
                }
            }
        }
    }
    else { // Renvoie un message d’erreur si ce n’est pas le cas
        printf("ERREUR: Impossible d’ouvrir le fichier.\n");
    }
    mon_flux.close(); // Ferme le fichier

    MPI_Finalize();

    // Fin du chrono
    auto finish = chrono::high_resolution_clock::now();
    double t = chrono::duration_cast<chrono::seconds>(finish-start).count();

    printf("Cela a pris %f secondes \n", t);

    return 0;
}
