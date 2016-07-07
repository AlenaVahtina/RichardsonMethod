#include <iostream>
#include "functions.h"
#include "gnuplot.h"
#include "plots.h"
using namespace std;

int main()
{
    //удаление файлов данных и графиков предыдущей работы программы
    system("rm *.png *.dat");
    Richardson Rid;
    Plots Plot;

    //инициализация данных
    double a=0,b=1,h;
    int N; //чтение колличества ячеек
    double s;//колличесство итераций
    vector<double> f;//функция, правая часть
    vector<double> y;//искомое расспределение темпиратур, вектор неизвестных
    vector<double> deltak;// вектор ошибок
    vector<vector<double>> Matrix;

    //какие графики итерации нужно выводить, kr=0-ни одного, kr=1-(s-1)-графики кратные kr, kr=s -только последнюю итерацию
    int kr=4;
    Plot.setKr(kr);

    //число ячеек
    N=10;
    cout<<"The number of cells \n"<<N<<'\n';

    //настройки а и b по умолчанию
    Rid.setA(a);
    Rid.setB(b);

    //шаг
    h=(b-a)/N;

    //число итераций
    s=16;
    cout<<"Enter the number of iterations\n"<<s<<'\n';
    Rid.setS(s);

    //заполнение матрицы СЛАУ A
    Matrix.resize(N);
    for (int i=0; i<=N; i++){
        Matrix[i].resize(N);
    }
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            if(j>0) Matrix[j-1][j]=1/(h*h);
            Matrix[j][j]=-2/(h*h);
            if(j<N-1) Matrix[j+1][j]=1/(h*h);
        }
    }

    //нулевое заполнение y
    y.resize(N);
    for (int i=0; i<N; i++){
        y[i]=1;
    }
    y[0]=0;
    y[N-1]=0;


    //считать значения вектора f
    f.resize(N);
    for (int i=0; i<N; i++){
        f[i]=0;
    }

     //итерация
     Plot.setKr(kr);
//    Rid.ItartionR(y, Matrix,f,kr);
    Rid.ItartionRFin(y, Matrix,f,kr);
//    Rid.ItartionRWithGer1(y,Matrix,f,kr);

    deltak=Rid.getErrors();


    //вывод у
    cout<<endl;
    for (int i=0; i<N; i++){
        cout<<y[i]<<"  ";
    }
    cout<<endl;

//    for (int i=0; i<=s; i++){
//        cout<<deltak[i]<<"  ";
//    }
//    cout<<endl;

    Plot.YPlot(y);
    //Plot.PlotWithE(deltak);
   // Plot.AveragePlot(deltak);
    return 0;
}

