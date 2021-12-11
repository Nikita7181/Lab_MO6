#include "simplex.h"
#include <sstream>

template <typename T>
void swap(T  &a, T &b)
{
    T c = b;
    b = a;
    a = c;
}

template <typename T>
std::string toString (const T val)
{
    std::stringstream ss;
    ss << val;
    return ss.str();
}

void swap(std::string &a, std::string &b)
{
    std::string c = b;
    b = a;
    a = c;
}

Simplex::Simplex()
{
    sizex = sizey = sizec = sizea= sizeb = 0;
    funcMax = 0;
    canSolve=true;
    infinite=false;
}

Simplex::Simplex(std::vector<std::vector<double>> a, std::vector<double> b, std::vector<double> c)
{
    sizey = a.size();
    sizeb = b.size();
    sizea = a[0].size();
    sizex = sizea;
    sizec = c.size();
    funcMax = 0;
    canSolve = true;
    infinite = false;
    A.resize(sizey, std::vector<double>(sizex));
    B.resize(b.size());
    C.resize(sizec);
    origA.resize(sizey, std::vector<double>(sizea));

    origB.resize(b.size());
    origC.resize(c.size());

    for (int i = 0; i < sizey; i++)
        for(int j = 0; j < sizea; j++)
            A[i][j] = a[i][j];

    for (int i = 0; i < b.size(); i++)
        B[i] =b[i];
    for (int i = 0; i < c.size(); i++)
    {
        C[i] = c[i];
    }

    titleCol.resize(sizey);
    titleRow.resize(sizex);

    for (int i = 0; i < sizey; i++)
        titleCol[i] = 'u' + toString(i + sizea + 1);
    for (int i = 0; i < sizex; i++)
        titleRow[i]='u'+ toString(i+1);


    for (int i = 0; i < b.size(); i++)
        origB[i] = B[i];
    for (int i = 0; i < c.size(); i++)
        origC[i] = C[i];
    for (int i = 0; i < sizey; i++)
        for(int j = 0; j < sizea; j++)
            origA[i][j] = a[i][j];
}

void Simplex::makeDual()
{
    std::vector <std::vector<double>> newA;
    std::vector<double> newB;
    std::vector<double> newC;

    newB.resize(C.size());
    newC.resize(B.size());
    newA.resize(sizex, std::vector<double>(sizey));

    for (int i = 0; i < newC.size(); i++)
        newC[i] = B[i] + 0.0;

    for (int i = 0; i < newB.size(); i++)
        newB[i] = C[i] + 0.0;

    B.resize(newB.size());
    C.resize(newC.size());
    for (int i = 0; i < newB.size(); i++)
        B[i] = -1 * newB[i];
    for (int i = 0; i < newC.size(); i++)
        C[i] = -1 * newC[i];

    for (int i = 0; i < sizey; i++)
    {
        for (int j = 0; j < sizex; j++)
        {
            newA[j][i] = -1 * A[i][j] + 0.0;
        }
    }

    A = newA;

    swap(sizex,sizey);
    titleCol.resize(sizey);
    titleRow.resize(sizex);

    for (int i = 0; i < sizey; i++)
        titleCol[i] = 'u' + toString(i + sizex + 1);
    for (int i = 0; i < sizex; i++)
        titleRow[i]='u'+ toString(i+1);
}

bool Simplex::optimal() {
    bool b_allPositive = false;
    for (int j = 0; j < B.size(); j++)
    {
        if (B[j] < 0)
            return false;
    }
    b_allPositive = true;

    for (int i = 0; i < sizea; i++)
        if (C[i] > 0) {
            return false;
        }

    return true;
}

int Simplex::findColumn(std::vector<double> tB, std::vector<double> tC) {
    int index = 0;

    double max = tC[0];
    for (int i = 0; i < tC.size(); i++)
        if (tC[i] >= max)
        {
            max = tC[i];
            index = i;
        }


    return index;
}

int Simplex::findRow(std::vector<std::vector<double>> tA, std::vector<double> tB, std::vector<double> tC) {

    int index = 0;
    int col = findColumn(tB, tC);
    double min = tB[index] / tA[index][col];

    int negCount = 0;
    for (int i = 1; i < tB.size(); i++)
    {
        if (tA[i][col] < 0)
            negCount++;
    }
    if (negCount == tB.size())
    {
        canSolve = false;

    }

    for (int i = 1; i < tB.size(); i++)
    {
        if (tB[i] / tA[i][col] < min) //&& A[i][col]!=0)
        {
            min = fabs(tB[i] / tA[i][col]);
            index = i;
        }
    }

    return index;
}

void Simplex::doTransform(int row, int col)
{
    double elem = A[row][col];
    double b = B[row];

    std::vector<double> oldRow; oldRow.resize(sizex);
    std::vector<double> oldCol; oldCol.resize(sizey);
    std::vector<double> newRow; newRow.resize(sizex);
    std::vector<double> newCol; newCol.resize(sizey);

    funcMax = funcMax - (C[col]*(B[row]/elem));

    double  new_elem = 1 / elem;

    for (int i = 0; i < sizex; i++)
        oldRow[i]=A[row][i];

    for (int i = 0; i < sizey; i++)
        oldCol[i]=A[i][col];

    for (int i = 0; i < sizex; i++)
        newRow[i] = oldRow[i]/elem;

    for (int i = 0; i < sizey; i++)
        newCol[i] = -1 * oldCol[i]/elem;

    for (int i = 0; i < sizey; i++)
        for(int j = 0; j < sizex; j++)
            A[i][j] = A[i][j] - oldCol[i]/elem * oldRow[j] + 0.0;

    for (int i = 0; i < sizey; i++)
        A[i][col] = newCol[i] + 0.0;

    for (int i = 0; i < sizex; i++)
        A[row][i] = newRow[i] + 0.0;

    B[row]=B[row]/elem + 0.0;
    for (int i = 0; i < sizey; i++)
        if (i!=row)
            B[i]= B[i] - oldCol[i]* B[row] + 0.0;

    for (int i = 0; i < sizex; i++)
        if (i!=col)
            C[i]= C[i]-C[col]/elem*oldRow[i] + 0.0;

    C[col] = -1 * C[col]/elem + 0.0;

    A[row][col] = new_elem;


}

std::vector<double> Simplex::calculate()
{
    std::vector<double> result;
    int step = 0;
    while ( !optimal()) {
        //endless(findColumn());
        if (infinite )
            break;
        step++;
        print2();
        int row = findRow(A, B, C);
        int col = findColumn(B, C);

        std::cout << "step " << step << "\n";
        //std::cout << row << "  " << col << "  " << A[row][col] << "\n";
        getElem( row, col);
        //std::cout << row << "  " << col << "  " << A[row][col] << "\n";

        swap(titleCol[row], titleRow[col]);
        doTransform(row, col);
        //doTransform(findRow(), findColumn());
    }
    if (canSolve && !infinite) {
        if(funcMax < 0)
        std::cout << "Function maximum value: " << -funcMax << "\n";
        else
            std::cout << "Function maximum value: " << funcMax << "\n";
        print2();

        //std::vector<double> result;
        result.resize(sizex, 0);

        for (int i = 0; i < sizeb; i++)
        {
            if (static_cast<int> (titleCol[i][1]) - 48 <= sizeb) {
                result[static_cast<int> (titleCol[i][1]) - 48 - 1] = B[i];
            }
            //std::cout << titleCol[i] << " = " << B[i] << "  ";
        }

        for (int i = 0; i < sizex; i++) {
            std::cout << "u" << i + 1<< " = " << result[i] << " " ;
        }
        std::cout << "\n";
        return result;
    }
    else if (infinite)
    {
        std::cout << "Function has infinite number of solutions: " <<  "\n";
    }
    else
    {
        std::cout << "Cant solve...";
    }


    return result;
}

double Simplex::getFMax()
{
    return funcMax;
}

void Simplex::print2()
{
    std::cout << std::fixed << std::setprecision(3) << " \t  S0\t";
    for (int i = 0; i < sizex; i++) {
        std::cout << "   " << titleRow[i] << "\t";

    }
    std::cout << "\n";
    for (int i = 0; i < sizey; i++) {
        std::cout << titleCol [i] << "\t";
        B[i] < 0 ?  std::cout  << B[i]: std::cout << " " << B[i];
        std::cout << "\t";
        for (int j = 0; j < sizex; j++)
        {
            A[i][j] < 0 ? std::cout << A[i][j] << "\t" : std::cout << " " << A[i][j] << "\t" ;
        }
        std::cout << "\n";
    }
    std::cout <<"F\t"; funcMax < 0 ?  std::cout  << funcMax : std::cout << " " << funcMax;
    std::cout << "\t";
    for (int i = 0; i < sizex; i++)
    {
        C[i] < 0 ? std::cout << C[i] << "\t" : std::cout << " " << C[i] << "\t" ;
    }
    std::cout << "\n\n";

}

void Simplex::getElem( int & row, int & col)
{
    for (int i = 0; i < B.size(); i++)
        if (B[i] < 0)
        {
            for (int j = 0; j < C.size(); j++)
                if (A[i][j] < 0)
                {
                    row = i;
                    col = j;
                    return;
                }
        }
}
