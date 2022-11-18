#!/bin/bash

ln -sf lmodel_relxill.dat lmodel.dat

sed -i "s,#define RELXILL_TABLE_PATH.*,#define RELXILL_TABLE_PATH "'"'`pwd`'"'"," relbase.h
# sed -i "s,#define RELXILL_TABLE_PATH.*,#define RELXILL_TABLE_PATH "'"'`$(RELXILL_TABLE_PATH)`'"'"," relbase.h
# sed -i.ori "s,#define RELXILL_TABLE_PATH.*,#define RELXILL_TABLE_PATH "'"'"$RELXILL_TABLE_PATH"'"'"," relbase.h

echo "initpackage relxill_nk lmodel_relxill.dat `pwd` \nexit" | xspec
echo "lmod relxill_nk . \nexit" | xspec

rm -f *~ *.o *FunctionMap.* lpack_* *.mod Makefile
