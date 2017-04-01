#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <iostream>
#include <../matrix/basematrix.h>
#include "../functions.h"
#include "../gnuplot.h"
#include "../plots.h"
#include "../matrix/normalmatrix.h"
#include "../matrix/crsmatrix.h"
#include "../common.h"
#include "../richardsonMethodWithChebyshevOrderedSetOfParameters.h"

#include <QDir>
#include <QDebug>

using namespace std;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    system("rm *.png *.dat");


    RichardsonMethod Rid;
    Plots Plot;


    //инициализация данных
    double a=0,b=1,step;//левый и правый конец, шаг
    int nAmountPoints; //чтение колличества ячеек
    double iterationNomber;//колличесство итераций
    vector<double> f;//функция, правая часть
    vector<double> y;//искомое расспределение темпиратур, вектор неизвестных
    vector<double> deltak;// вектор ошибок
    vector <vector<double> > Matrix;


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
    iterationNomber=16;
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

    //считать параметры q и p
    Rid.setP(0);
    Rid.setQ(0);

     //вычисление у (основное решение задачи)
     BaseMatrix *testslau=new NormalMatrix(Matrix);
//     Rid.computeResultVectorForE(y, testslau,f,fold);
//     Rid.computeResultVectorForC(y, testslau,f,fold);
//     Rid.computeResultVectorForEWithRivalProcess(y, testslau, f, fold);
       Rid.computeResultVectorForNotEWithRivalProcess(y, testslau, f, fold);

     deltak=Rid.getErrors();


    //вывод у
    cout<<endl;
    for (int i=0; i<nAmountPoints; i++){
        cout<<y[i]<<"  ";
    }
    cout<<endl;


    ui->pushButton_4->click();
}

void MainWindow::showImage(QString path){
        QPixmap pixmap(path);
        ui->label->setPixmap(pixmap);
       // ui->label->setMask(pixmap.mask());
        ui->label->repaint();
}

void MainWindow::on_pushButton_4_clicked()
{
    QStringList filters;
    QDir dir(".");
    filters << "*.png" << "*.jpg" << "*.bmp";
    allFiles = dir.entryList(filters, QDir::Files|QDir::NoDotAndDotDot);
    if(allFiles.size()!=0){
        showImage(allFiles.at(0));
        current=0;
    }
}

void MainWindow::on_pushButton_2_clicked()
{
    if(allFiles.size()!=0){
        if(allFiles.size()-1 == current)
            current=0;
        else
            current++;
        showImage(allFiles.at(current));
    }

}

void MainWindow::on_pushButton_3_clicked()
{
    if(allFiles.size()!=0){
        if(0 == current)
            current=allFiles.size()-1;
        else
            current--;
        showImage(allFiles.at(current));
    }
}
