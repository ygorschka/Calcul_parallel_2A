#ifndef _DATA_FILE_H

#include <string>

class DataFile {
private:
    std::string _file_name, _wich_scenario;

    int _Nx, _Ny, _iBeg, _iEnd, _rank, _nproc, _size_vect, _dec;
    double _Lx, _Ly, _D, _dt, _Tmax;

public:
    DataFile(std::string file_name);

    const int Get_Nx() const {return _Nx;};
    const int Get_Ny() const {return _Ny;};
    const double Get_Lx() const {return _Lx;};
    const double Get_Ly() const {return _Ly;};
    const double Get_D() const {return _D;};
    const double Get_dt() const {return _dt;};
    const double Get_Tmax() const {return _Tmax;};
    const std::string Get_scenario() const {return _wich_scenario;};
    const int Get_rank() const {return _rank;};
    const int Get_nproc() const {return _nproc;};
    const int Get_iBeg() const {return _iBeg;};
    const int Get_iEnd() const {return _iEnd;};
    const int Get_dec() const {return _dec;};
    const int Get_size_vect() const {return _size_vect;};

    void Read_data_file();
    void charge_c(int me, int N, int np, int *iBeg, int *iEnd);

};

#define _DATA_FILE_H
#endif
