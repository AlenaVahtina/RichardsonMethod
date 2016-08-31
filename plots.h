#ifndef PLOTS_H
#define PLOTS_H

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <fstream>
#include "gnuplot.h"

using namespace std;

class Plots
{
public:
    Plots();
    int fold;
    void setfold(int _fold){fold=_fold;}


    //построение графика у
    void YPlot(vector<double> &y){
        ofstream f;
        f.open("output.dat");
        if (!f.is_open())
          {
            cout << "Error opening file output.dat.\n";
            exit(EXIT_FAILURE);
          }

        for (int i=0; i<y.size(); i++){
            f<<y[i]<<"   "<<i<<endl;
        }
        f.close();
        Gnuplot plot;
        plot("set terminal png size 900,800 enhanced font \"Helvetica,20\"");
        plot("set output 'output.png'");
        plot("plot 'output.dat' using 2:1 with lines title 'y(i)'");
    }


    //построение графика у на каждой итерации с пресвоением нового имени
    void IteratPlot(vector<double> &y,std::string plotname,std::string filename){
        ofstream f;
            f.open(filename);
            if (!f.is_open())
              {
                cout << "Error opening file output.dat.\n";
                exit(EXIT_FAILURE);
              }
            for (int i=0; i<y.size(); i++){
                f<<y[i]<<"   "<<i<<endl;
            }
            f.close();
            Gnuplot plot;
            //f<<"set terminal size 900,800"<<endl<<"set output '"<<i<<".png'"<<endl
            plot("set terminal png size 900,800 enhanced font \"Helvetica,20\"");
            plot("set output '"+plotname+"'");
            plot("plot '"+filename+"' using 2:1 with lines title 'y(i)'");
    }


    //построение графика ошибки
    void PlotWithE(vector<double> &deltak){
        ofstream f;
        f.open("outdelta.dat");
        if (!f.is_open())
          {
            cout << "Error opening file output.dat.\n";
            exit(EXIT_FAILURE);
          }

        for (int i=0; i<deltak.size(); i++){
            f<<deltak[i]<<"   "<<log10(deltak[i])<<"  "<<i<<endl;
        }
        f.close();
        Gnuplot plot;
        plot("set terminal png size 900,800 enhanced font \"Helvetica,20\"");
        plot("set output 'outdelta.png'");
        plot("plot 'outdelta.dat' using 3:2 with linespoints title 'lg(delta k(s))'");
    }


    //построение приведенного графика ошибки
    void AveragePlot(vector<double> &deltak){
        ofstream f;
        f.open("AveragePlot.dat");
        if (!f.is_open())
          {
            cout << "Error opening file output.dat.\n";
            exit(EXIT_FAILURE);
          }
        if (deltak.size()==0)return;
        f<<deltak[0]<<"   "<<log10(deltak[0])<<"  "<<0<<endl;
        double predk=deltak[0];
        for (int i=4; i<deltak.size()-1; i++){
            if ((deltak[i]<predk) && (i%4==0)) {
                f<<deltak[i]<<"   "<<log10(deltak[i])<<"  "<<i<<endl;
                predk=deltak[i];
            }
        }
        f.close();
        Gnuplot plot;
        plot("set terminal png size 900,800 enhanced font \"Helvetica,20\"");
        plot("set output 'AveragePlot.png'");
        plot("set xtics 4");
        plot("set grid xtics mxtics ytics");
        plot("set format x '' ");
        plot("plot 'AveragePlot.dat' using 3:2 with linespoints title 'Average log (delta k(s))'");
    }


    //построение графиков для конкурирующих процессов
    void AveragePlotDoble(vector<double> &deltak,std::string plotname,std::string filename){
        ofstream f;
        f.open(filename);
        if (!f.is_open())
          {
            cout << "Error opening file output.dat.\n";
            exit(EXIT_FAILURE);
          }
        f<<deltak[0]<<"   "<<log10(deltak[0])<<"  "<<0<<endl;
        double predk=deltak[0];
        for (int i=4; i<deltak.size()-1; i++){
            if ((deltak[i]<predk) && (i%4==0)) {
                f<<deltak[i]<<"   "<<log10(deltak[i])<<"  "<<i<<endl;
                predk=deltak[i];
            }
        }
        f.close();
        Gnuplot plot;
        plot("set terminal png size 900,800 enhanced font \"Helvetica,20\"");
        plot("set output '"+plotname+"'");
        plot("set xtics 4");
        plot("set grid xtics mxtics ytics");
        plot("set format x '' ");
        plot("plot '"+filename+"' using 3:2 with linespoints title 'Average log (delta k(s))'");
    }


    //построение приведенных графиков для конкурирующих процессов
    void AveragePlotDoble2(vector<double> &deltak, vector<double> &deltak2, std::string plotname,std::string filename){
        ofstream f;
        f.open(filename);
        if (!f.is_open())
          {
            cout << "Error opening file output.dat.\n";
            exit(EXIT_FAILURE);
          }
        f<<deltak[0]<<"   "<<log10(deltak[0])<<"   "<<deltak2[0]<<"  "<<log10(deltak2[0])<<"   "<<0<<endl;
        double predk=deltak[0];
        for (int i=4; i<deltak.size()-1; i++){
            if ((deltak[i]<predk) && (deltak2[i]<predk) && (i%4==0)) {
                f<<deltak[i]<<"   "<<log10(deltak[i])<<"   "<<deltak2[i]<<"  "<<log10(deltak2[i])<<"   "<<i<<endl;
                predk=deltak[i];
            }
        }
        f.close();
        Gnuplot plot;
        plot("set terminal png size 900,800 enhanced font \"Helvetica,20\"");
        plot("set output '"+plotname+"'");
        plot("set xtics 4");
        plot("set grid xtics mxtics ytics");
        plot("set format x '' ");
        plot("plot \""+filename+"\" using 5:2 with linespoints ,\\");
        plot("  \""+filename+"\" using 5:4 with linespoints");
    }

};
#endif // PLOTS_H
