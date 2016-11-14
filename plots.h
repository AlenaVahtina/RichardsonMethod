#ifndef PLOTS_H
#define PLOTS_H

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>

#include "gnuplot.h"

using namespace std;

class Plots
{
public:
    Plots();
    int fold;
    void setfold(int _fold){fold=_fold;}


    //построение графика у
    void yPlot(vector<double> &y);


    //построение графика у на каждой итерации с пресвоением нового имени
    void iteratPlot(vector<double> &y,std::string plotname,std::string filename);


    //построение графика ошибки
    void plotWithError(vector<double> &deltak);


    //построение приведенного графика ошибки
    void averagePlot(vector<double> &deltak);


    //построение графиков для конкурирующих процессов
    void averagePlotDoble(vector<double> &deltak,std::string plotname,std::string filename);


    //построение приведенных графиков для конкурирующих процессов
    void averagePlotDoble2(vector<double> &deltak, vector<double> &deltak2, std::string plotname,std::string filename);

};
#endif // PLOTS_H
