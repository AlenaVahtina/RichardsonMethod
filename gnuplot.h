// gnuplot.h
#ifndef GNUPLOT_H_
#define GNUPLOT_H_

#include <cstdio>
#include <string>
#include <iostream>

#ifdef WIN32
    #define GNUPLOT_NAME "pgnuplot -persist"
#else
    #define GNUPLOT_NAME "gnuplot -persist"
#endif
static int ggg=0;
using std::string;
using std::cerr;

class Gnuplot
{
public:
    Gnuplot() {
        #ifdef WIN32
            gnuplotpipe = _popen(GNUPLOT_NAME, "w");
        #else
            gnuplotpipe  = popen(GNUPLOT_NAME, "w");
        #endif

        if (!gnuplotpipe)
        {
            cerr << ("Gnuplot not found !");
        }
    }
    ~Gnuplot(){
        fprintf(gnuplotpipe,"exit\n");

        #ifdef WIN32
           _pclose(gnuplotpipe);
        #else
            pclose(gnuplotpipe);
        #endif
 //       std::cout<<ggg<<std::endl;
 //       ggg++;
    }
    void operator ()(const string & command)
    {
        fprintf(gnuplotpipe,"%s\n",command.c_str());
        fflush(gnuplotpipe); //без fflush ничего рисоваться не будет
    }

protected:
    FILE *gnuplotpipe;
};

#endif // #ifndef _GNUPLOT_H_
