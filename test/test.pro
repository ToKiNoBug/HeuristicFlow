TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += ../ \
                               D:/CppLibs/eigen-3.4.0

SOURCES += \
        imp_TestingGenetic.cpp \
        imp_TestingMatrix.cpp \
        main.cpp \
        testNsga2.cpp

HEADERS += \
    def_TestingGenetic.h \
    def_TestingMatrix.h \
    testNsga2.h
