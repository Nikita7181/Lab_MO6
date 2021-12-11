#pragma once
#include <iostream>
#include <vector>
#include <iomanip>
#include <string>
#include <fstream>
#include "nlohmann/json.hpp"


class Simplex
{
    int sizex, sizey, sizec, sizea, sizeb;
    std::vector <std::vector<double>> A;
    std::vector<double> B;
    std::vector<double> C;
    std::vector <std::vector<double>> origA;
    std::vector<double> origB;
    std::vector<double> origC;
    double funcMax;
    std::vector<std::string> titleRow;
    std::vector<std::string> titleCol;

    bool canSolve;
    bool infinite;



public:
    Simplex();
    Simplex(std::vector <std::vector<double>> a, std::vector<double> b, std::vector<double> c);
    bool optimal();
    int findRow(std::vector<std::vector<double>> tA, std::vector<double> tB, std::vector<double> tC);
    int findColumn(std::vector<double> tB, std::vector<double> tC);
    std::vector<double> calculate();
    void doTransform(int row, int col);
    void endless (int column);
    void getElem( int & row, int & col);
    void print();
    void print2();
    void makeDual();
    double getFMax();
};