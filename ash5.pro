TEMPLATE = app

CONFIG += console c++11

CONFIG -= app_bundle

CONFIG -= qt

SOURCES += \
    ash5.cpp \
    ash5suite.cpp

HEADERS += \
    ash5.h \
    dyn_utils.h

INCLUDEPATH += C:\cpp_environment\Eigen3\
INCLUDEPATH += C:\cpp_environment\boost_1_64_0\


DISTFILES += \
    suite-config-2.ini
