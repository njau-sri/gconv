QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = gconv-gui
TEMPLATE = app

DEFINES += QT_NO_CAST_FROM_ASCII QT_NO_CAST_TO_ASCII
DEFINES += QT_DEPRECATED_WARNINGS

SOURCES += main.cpp \
        dialog.cpp

HEADERS += dialog.h

FORMS += dialog.ui
