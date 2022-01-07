TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += ../OptimTemplates \
                               D:/CppLibs/eigen-3.4.0

SOURCES += \
        imp_TestingGenetic.cpp \
        imp_TestingMatrix.cpp \
        main.cpp

HEADERS += \
    def_TestingGenetic.h \
    def_TestingMatrix.h
