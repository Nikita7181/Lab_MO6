#include <iostream>
#include "simplex.h"
#include "nlohmann/json.hpp"
#include <fstream>

using json = nlohmann::json;

double calc_g(double W)
{
    return 1/fabs(W);
}

double minmax (std::vector<std::vector<double>> A)
{
    double max = -999999;
    double min =  999999;

    int sizex = A[0].size();
    int sizey = A.size();

    for (int i = 0; i < sizey; i++)
    {

        for (int j = 0; j < sizex; j++)
        {
            if (A[i][j] < min)
                min = A[i][j];
        }
        if (max < min)
            max = min;
        min = 999999;
    }

    return fabs(max);
}

double maxmin (std::vector<std::vector<double>> A)
{
    double max = -999999;
    double min =  999999;

    int sizex = A[0].size();
    int sizey = A.size();

    for (int j = 0; j < sizex; j++)
    {
        for (int i = 0; i < sizey; i++)
        {
            if (A[i][j] > max)
                max = A[i][j];
        }
        if (min > max)
            min = max;
        max =  -999999;
    }

    return fabs(min);
}

std::vector<double> calcOptimalStrategy(std::vector<double> U, double g)
{
    std::vector<double> result;
    result.resize(U.size(), 0);

    for (int i = 0; i < U.size(); i++)
    {
        result[i] = U[i] * g;
    }
    return result;

}

bool checkSumXEqualToOne(std::vector<double> X)
{
    double sum = 0.0;
    for (int i = 0; i < X.size(); i++)
    {
        sum = sum + X[i];
    }

    double d = sum - 1.0;
    if ( fabs(d) < 0.00001)
        return true;

    return false;

}

void print_v (std::vector<double> v, int start = 0, int end = 0)
{

    if (end == 0)
        end=v.size()-1;
    for (int i = start; i <= end; i++)
        if (v[i] < 0 )
            std::cout << "y" << i + 1 << "=" << std::fixed << std::setprecision(5)  << v[i] << "  ";
        else
            std::cout << "y" << i + 1 << "=" << std::fixed << std::setprecision(5) << " " << v[i] << "  ";
    std::cout << "\n";
}

int main(int argc, char * argv[])
{

    if (argc < 2)
    {
        throw std::runtime_error ("Path is required as parameter!");
    }
    std::string jsonPath1 = argv[1];

    std::ifstream file1{jsonPath1};
    if (!file1 )  {
        throw std::runtime_error{"unable to open json: " + jsonPath1};
    }

    json data1;

    try {
        file1 >> data1;
    }
    catch (std::exception e)
    {
        throw std::runtime_error{"Wrong json format"};
    }


    std::vector<double > c1;
    std::vector<double > b1;
    std::vector<std::vector<double>> a1;

    if (data1.at("C").is_array())
    {
        c1 = data1.at("C").get<std::vector<double>>();
        b1 = data1.at("B").get<std::vector<double>>();
        a1 = data1.at("A").get<std::vector<std::vector<double>>>();
    }

    Simplex s1 (a1, b1, c1);

    std::vector<double> U1 = s1.calculate();

    Simplex s2 (a1, b1, c1);
    s2.makeDual();

    std::vector<double> U2 = s2.calculate();

    double g1 = calc_g(s1.getFMax());
    std::vector<double> X1 = calcOptimalStrategy( U1, g1);
    std::cout << "\n\nOptimal strategy for Player A: \n";
    print_v(X1);

    std::string res1 = (checkSumXEqualToOne(X1) == true) ? "true" : "false";
    std::cout << "sum Xi = 1 : " << res1 <<"\n" ;

    double g2 = calc_g(s2.getFMax());
    std::vector<double> X2 = calcOptimalStrategy( U2, g2);
    std::cout << "Optimal strategy for Player B: \n";
    print_v(X2);

    std::string res2 = (checkSumXEqualToOne(X2) == true) ? "true" : "false";
    std::cout << "sum Xi = 1 : " << res2 << "\n";

    std::cout << "upper coast = " << minmax(a1) << "\n";
    std::cout << "lower coast = " << maxmin(a1) << "\n";
    std::cout << "Game cost =  " << calc_g(s1.getFMax()) << "\n";

    return 0;
}
