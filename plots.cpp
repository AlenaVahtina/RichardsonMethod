#include "plots.h"
#include <fstream>
Plots::Plots()
{
}

//построение графика у
void Plots::yPlot(vector<double> &y){
    ofstream f;
    f.open("output.dat");
    if (!f.is_open())
      {
        cout << "Error opening file output.dat.\n";
//        exit(EXIT_FAILURE);
        return;
      }

    for (int i=0; i<y.size(); i++){
        f<<y[i]<<"   "<<i<<endl;
    }
    f.close();
    Gnuplot plot;
    plot("set terminal postscript enhanced");
    plot("set terminal png size 900,800 enhanced font \"Helvetica,20\"");
    plot("set output 'output.png'");
    plot("set xlabel 's' ");
    plot("set ylabel 'y' ");
    plot("plot 'output.dat' using 2:1 with lines title 'y(i)'");
}


//построение графика у на каждой итерации с пресвоением нового имени
void Plots::iteratPlot(vector<double> &y,std::string plotname,std::string filename){
    ofstream f;
        f.open(filename.data());
        if (!f.is_open())
          {
            cout << "Error opening file output.dat.\n";
            //exit(EXIT_FAILURE);
            return;
          }
        for (int i=0; i<y.size(); i++){
            f<<y[i]<<"   "<<i<<endl;
        }
        f.close();
        Gnuplot plot;
        //f<<"set terminal size 900,800"<<endl<<"set output '"<<i<<".png'"<<endl
        plot("set terminal postscript enhanced");
        plot("set terminal png size 900,800 enhanced font \"Helvetica,20\"");
        plot("set xlabel 's' ");
        plot("set ylabel 'y' ");
        plot("set output '"+plotname+"'");
        plot("plot '"+filename+"' using 2:1 with lines title 'y(i)'");
}


//построение графика ошибки
void Plots::plotWithError(vector<double> &deltak){
    ofstream f;
    f.open("outdelta.dat");
    if (!f.is_open())
      {
        cout << "Error opening file output.dat.\n";
        //exit(EXIT_FAILURE);
        return;
      }

    for (int i=0; i<deltak.size(); i++){
        f<<deltak[i]<<"   "<<log10(deltak[i])<<"  "<<i<<endl;
    }
    f.close();
    Gnuplot plot;
    plot("set terminal postscript enhanced");
    plot("set terminal png size 900,800 enhanced font \"Helvetica,20\"");
    plot("set output 'outdelta.png'");
    plot("set xlabel 's' ");
    plot("set ylabel 'y' ");
    plot("plot 'outdelta.dat' using 3:2 with linespoints title 'lg(delta k(s))'");
}


//построение приведенного графика ошибки
void Plots::averagePlot(vector<double> &deltak){
    ofstream f;
    f.open("AveragePlot.dat");
    if (!f.is_open())
      {
        cout << "Error opening file output.dat.\n";
        //exit(EXIT_FAILURE);
        return;
      }
    if (deltak.size()==0)return;
    f<<deltak[0]<<"   "<<log10(deltak[0])<<"  "<<0<<endl;
    double predk=deltak[0];
    for (int i=3; i<deltak.size()-1; i+=4){
        if ((deltak[i]<predk)) {
            f<<deltak[i]<<"   "<<log10(deltak[i])<<"  "<<i<<endl;
            predk=deltak[i];
        }
    }
    f.close();
    Gnuplot plot;
    plot("set terminal postscript enhanced");
    plot("set terminal png size 900,800 enhanced font \"Helvetica,20\"");
    plot("set output 'AveragePlot.png'");
    plot("set xtics 4");
    plot("set grid xtics mxtics ytics");
    plot("set format x '' ");
    plot("set xlabel 's' ");
    plot("set ylabel 'y' ");
    plot("plot 'AveragePlot.dat' using 3:2 with linespoints title 'Average log (delta k(s))'");
}


//построение графиков для конкурирующих процессов
void Plots::averagePlotDoble(vector<double> &deltak,std::string plotname,std::string filename){
    ofstream f;
    f.open(filename.data());
    if (!f.is_open())
      {
        cout << "Error opening file output.dat.\n";
        //exit(EXIT_FAILURE);
        return;
      }
    f<<deltak[0]<<"   "<<log10(deltak[0])<<"  "<<0<<endl;
    double predk=deltak[0];
    for (int i=3; i<deltak.size()-1; i+=4){
        if ((deltak[i]<predk) ) {
            f<<deltak[i]<<"   "<<log10(deltak[i])<<"  "<<i<<endl;
            predk=deltak[i];
        }
    }
    f.close();
    Gnuplot plot;
    plot("set terminal postscript enhanced");
    plot("set terminal png size 900,800 enhanced font \"Helvetica,20\"");
    plot("set output '"+plotname+"'");
    plot("set xtics 4");
    plot("set grid xtics mxtics ytics");
    plot("set format x '' ");
    plot("set xlabel 's' ");
    plot("set ylabel 'y' ");
    plot("plot '"+filename+"' using 3:2 with linespoints title 'Average log (delta k(s))'");
}


//построение приведенных графиков для конкурирующих процессов
void Plots::averagePlotDoble2(vector<double> &deltak, vector<double> &deltak2, std::string plotname,std::string filename){
    ofstream f;
    f.open(filename.data());
    if (!f.is_open())
      {
        cout << "Error opening file output.dat.\n";
        //exit(EXIT_FAILURE);
        return;
      }
    f<<deltak[0]<<"   "<<log10(deltak[0])<<"   "<<deltak2[0]<<"  "<<log10(deltak2[0])<<"   "<<0<<endl;
    double predk=deltak[0];
    for (int i=3; i<deltak.size()-1; i+=4){
        if ((deltak[i]<predk) && (deltak2[i]<predk) ) {
            f<<deltak[i]<<"   "<<log10(deltak[i])<<"   "<<deltak2[i]<<"  "<<log10(deltak2[i])<<"   "<<i<<endl;
            predk=deltak[i];
        }
    }
    f.close();
    Gnuplot plot;
    plot("set terminal postscript enhanced");
    plot("set terminal png size 900,800 enhanced font \"Helvetica,20\"");
    plot("set output '"+plotname+"'");
    plot("set xtics 4");
    plot("set grid xtics mxtics ytics");
    plot("set format x '' ");
    plot("set xlabel 's' ");
    plot("set ylabel 'y' ");
    plot("plot \""+filename+"\" using 5:2 with linespoints ,\\");
    plot("  \""+filename+"\" using 5:4 with linespoints");
}
