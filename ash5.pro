TEMPLATE = app

CONFIG += console c++11

CONFIG -= app_bundle

CONFIG -= qt

SOURCES += \
    ash5.cpp \
#    ash5suite.cpp \
#    configuration.cpp \
#    run_info.cpp \
    ash5tests.cpp

HEADERS += \
    ash5.h \
    dyn_utils.h \
#    configuration.h \
    newt_sol.h \
    num_jac.h \
#    run_info.h
    misc_defs.h

INCLUDEPATH += C:\cpp-libs\Eigen3\
INCLUDEPATH += C:\cpp-libs\boost_1_65_0\

DISTFILES += \
#    suite-config-2.ini \
    make-tests.py
