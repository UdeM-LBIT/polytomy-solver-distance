#-------------------------------------------------
#
# Project created by QtCreator 2013-09-09T14:20:36
#
#-------------------------------------------------

INCLUDEPATH += src
DEPENDPATH += src

INCLUDEPATH += /usr/include/python2.7

unix:LIBS += -lboost_python -lpython2.7

TARGET = PolytomySolver
#CONFIG   += console
#CONFIG   -= app_bundle

TEMPLATE = app

DEFINES -= UNICODE
DEFINES += _MBCS
QMAKE_CXXFLAGS -= DUNICODE
QMAKE_CXXFLAGS += -std=c++0x -fPIC

SOURCES += main.cpp \
    src/trees/treeiterator.cpp \
    src/trees/treeinfo.cpp \
    src/trees/polysolver.cpp \
    src/trees/paralogycorrector.cpp \
    src/trees/node.cpp \
    src/trees/newicklex.cpp \
    src/trees/polysolver_distance.cpp \
    src/trees/polysolver_nad.cpp \
    src/trees/genespeciestreeutil.cpp \
    main_tests.cpp


HEADERS += \
    src/trees/treeiterator.h \
    src/trees/treeinfo.h \
    src/trees/polysolver.h \
    src/trees/paralogycorrector.h \
    src/trees/node.h \
    src/trees/newicklex.h \
    src/div/util.h \
    src/div/define.h \
    src/trees/polysolver_distance.h \
    src/trees/polysolver_nad.h \
    src/trees/genespeciestreeutil.h \
    src/div/tinydir.h
