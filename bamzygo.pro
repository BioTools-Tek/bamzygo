#-------------------------------------------------
#
# Project created by QtCreator 2013-05-24T11:52:28
#
#-------------------------------------------------

QT       += core
QT       -= gui

TARGET = bamzygo
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

LIBS += -L -lsamtools -lbam -lz

SOURCES += main.cpp

HEADERS += \
    Mpileup.h \
    ParseArgs.h \
    Zygosity.h
