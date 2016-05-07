#include <iostream>
#include "functions.h"
#include "gnuplot.h"
#include "plots.h"
using namespace std;

int main()
{

    Richardson Rid;
    Plots Plot;

    //инициализация данных
    double a,b,h;
    int N; //чтение колличества ячеек
    double s;//колличесство итераций
    vector<double> f;//функция, правая часть
    vector<double> y;//искомое расспределение темпиратур, вектор неизвестных
    vector<double> deltak;// вектор ошибок
    vector<vector<double>> Matrix;
//    vector<vector<double> >MatrixB;
//    vector<vector<double> >MatrixC;

    //какие графики итерации нужно выводить, kr=0-ни одного, kr=1-(s-1)-графики кратные kr, kr=s -только последнюю итерацию
    int kr=4;
    Plot.setKr(kr);

    //число ячеек
    N=100;
    cout<<"The number of cells \n"<<N<<'\n';

    //настройки а и b по умолчанию
    Rid.setA(0);
    Rid.setB(1);

    //шаг
    h=(b-a)/N;

    //число итераций
    s=64;
    cout<<"Enter the number of iterations\n"<<s<<'\n';
    Rid.setS(s);

    //заполнение матрицы СЛАУ A
    Matrix.resize(N+1);
    for (int i=0; i<=N; i++){
        Matrix[i].resize(N+1);
    }
    for (int i=0; i<=N; i++){
        for (int j=0; j<=N; j++){
            if(j>0)Matrix[j-1][j]=1;
            Matrix[j][j]=-2;
            if(j<N-1)Matrix[j+1][j]=1;
        }
    }

    //нулевое заполнение y
    y.resize(N+1);
    for (int i=0; i<=N; i++){
        y[i]=1;
    }
    y[0]=0;
    y[N]=0;

    //считать значения вектора f
    f.resize(N+1);
    for (int i=0; i<=N; i++){
        f[i]=0;
    }


//    MatrixB.resize(N+1);
//    for (int i=0; i<=N; i++){
//        MatrixB[i].resize(N+1);
//    }
//    Rid.CreateB(Matrix, MatrixB);
//    Rid.ReB(MatrixB);

//    MatrixC.resize(N+1);
//    for (int i=0; i<=N; i++){
//        MatrixC[i].resize(N+1);
//    }


    //итерация
    Plot.setKr(kr);
    Rid.ItartionR(y, Matrix, f,kr);
//    Rid.ItartionRWithGer(y, Matrix, f,kr);
    deltak=Rid.getErrors();


    //вывод у
    cout<<endl;
    for (int i=0; i<=N; i++){
        cout<<y[i]<<"  ";
    }
    cout<<endl;

    for (int i=0; i<=s; i++){
        cout<<deltak[i]<<"  ";
    }
    cout<<endl;

    Plot.YPlot(y);
    Plot.PlotWithE(deltak);
    Plot.AveragePlot(deltak);
    return 0;
}

