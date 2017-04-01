#-------------------------------------------------
#
# Project created by QtCreator 2017-04-01T21:36:24
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = gui
TEMPLATE = app

CONFIG += c++11
QMAKE_CXXFLAGS += -std=c++11

SOURCES += main.cpp\
        mainwindow.cpp\
        ../plots.cpp\
        ../matrix/basematrix.cpp\
        ../matrix/crsmatrix.cpp\
        ../matrix/normalmatrix.cpp\
        ../richardsonMethodWithChebyshevOrderedSetOfParameters.cpp\
        ../common.cpp

HEADERS  += mainwindow.h\
        ../plots.h\
        ../matrix/basematrix.h\
        ../matrix/crsmatrix.h\
        ../matrix/normalmatrix.h\
        ../richardsonMethodWithChebyshevOrderedSetOfParameters.h\
        ../common.h\
        ../gnuplot.h

FORMS    += mainwindow.ui
