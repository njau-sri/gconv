TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -std=c++11
QMAKE_LFLAGS += -static

SOURCES += \
        main.cpp \
    cmdline.cpp \
    geno.cpp \
    hmp.cpp \
    ped.cpp \
    util.cpp \
    vcf.cpp \
    gconv.cpp

HEADERS += \
    cmdline.h \
    geno.h \
    hmp.h \
    ped.h \
    split.h \
    util.h \
    vcf.h
