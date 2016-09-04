#include <iostream>
#include "richardsonslau.h"
#include "functions.h"
#include "gnuplot.h"
#include "plots.h"
#include "normalmatrix.h"
#include "crsmatrix.h"

using namespace std;

int main()
{
    //удаление файлов данных и графиков предыдущей работы программы
    system("rm *.png *.dat");


    RichardsonSLAU *SLAU=new normalmatrix(vector <vector <double> >({
                                                                    {2, 1, 0, 0, 0, 0},
                                                                    {1, 2, 1, 0, 0, 0},
                                                                    {0, 1, 2, 1, 0, 0},
                                                                    {0, 0, 1, 2, 1, 0},
                                                                    {0, 0, 0, 1, 2, 1},
                                                                    {0, 0, 0, 0, 1, 2}
                                                                    }));
    RichardsonSLAU *SLAU1=new crsmatrix(vector <int> ({0,2,5,8,11,14,16}), vector <int> ({0,1,2,1,2,3,2,3,4,3,4,5,4,5}), vector <double> ({2,1,1,2,1,1,2,1,1,2,1,1,2,1,1,2}));
//    SLAU->CreateC()->WriteMatrix();

//    cout<<"Types"<<SLAU->getType()<<"  "<<SLAU1->getType()<<'\n';
//    RichardsonSLAU *testslau2=new normalmatrix(Matrix);
//    testslau2>CreateB()->WriteMatrix();
//    testslau2->CreateC()->WriteMatrix();
    Richardson Rid;
    Plots Plot;

    //инициализация данных
    double a=0,b=1,step;//левый и правый конец, шаг
    int nAmountPoints; //чтение колличества ячеек
    double iterationNomber;//колличесство итераций
    vector<double> f;//функция, правая часть
    vector<double> y;//искомое расспределение темпиратур, вектор неизвестных
    vector<double> deltak;// вектор ошибок
    vector<vector<double>> Matrix;


    //какие графики итерации нужно выводить, fold=0-ни одного, fold=1-(iterationNomber-1)-графики кратные fold, fold=iterationNomber -только последнюю итерацию
    int fold=4;
    Plot.setfold(fold);


    //число ячеек (узлов)
    nAmountPoints=100;
    cout<<"The number of cells \n"<<nAmountPoints<<'\n';

    //настройки а и b по умолчанию
    Rid.setA(a);
    Rid.setB(b);

    //шаг
    step=(b-a)/nAmountPoints;

    //число итераций
    iterationNomber=32;
    cout<<"Enter the number of iterations\n"<<iterationNomber<<'\n';
    Rid.setS(iterationNomber);

    //заполнение матрицы СЛАУ A
    Matrix.resize(nAmountPoints);
    for (int i=0; i<nAmountPoints; i++){
        Matrix[i].resize(nAmountPoints);
    }
    for (int i=0; i<nAmountPoints; i++){
        for (int j=0; j<nAmountPoints; j++){
            if(j>0) Matrix[j-1][j]=1/(step*step);
            Matrix[j][j]=-2/(step*step);
            if(j<nAmountPoints-1) Matrix[j+1][j]=1/(step*step);
        }
    }

    //нулевое заполнение y
    y.resize(nAmountPoints);
    for (int i=0; i<nAmountPoints; i++){
        y[i]=1;
    }
    y[0]=0;
    y[nAmountPoints-1]=0;


    //считать значения вектора f
    f.resize(nAmountPoints);
    for (int i=0; i<nAmountPoints; i++){
        f[i]=0;
    }


     //вычисление у (основное решение задачи)
     RichardsonSLAU *testslau=new normalmatrix(Matrix);
//     Rid.computeResultVectorForE(y, testslau,f,fold);
     Rid.computeResultVectorForEWithRivalProcess(y, testslau,f,fold);
//   Rid.computeResultVectorForNotEWithRivalProcess(y,testslau,f,fold);

     deltak=Rid.getErrors();


    //вывод у
    cout<<endl;
    for (int i=0; i<nAmountPoints; i++){
        cout<<y[i]<<"  ";
    }
    cout<<endl;



    //построение графиков (по умолчанию выведены функции построения графиков для единичной матрицы без конкурирующих процессов)
    Plot.setfold(fold);
    Plot.YPlot(y);
    Plot.PlotWithE(deltak);
    Plot.AveragePlot(deltak);
    return 0;
}

