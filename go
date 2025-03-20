g++ -o eIDcut_DCfiducial eIDcut_DCfiducial.cpp \
        -I$HIPO/include -I$CLAS12ROOT/Clas12Root -I$CLAS12ROOT/Clas12Banks \
        -L$HIPO/lib -lhipo4 -L$CLAS12ROOT/lib -lClas12Banks -lClas12Root \
        `root-config --cflags --glibs`
