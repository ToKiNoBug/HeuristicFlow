TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += ../OptimTemplates \
                               D:/CppLibs/eigen-3.4.0

SOURCES += \
        imp_TestingFuns.cpp \
        main.cpp

HEADERS += \
    def_TestingFuns.h
