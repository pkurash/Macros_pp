#!/bin/bash

#names="LHC16g LHC16h  LHC16i  LHC16j  LHC16k  LHC16l  LHC16o  LHC16p"
#names="LHC17c LHC17e LHC17f LHC17h LHC17i LHC17k  LHC17l  LHC17m LHC17o LHC17r"
names="LHC18b LHC18d LHC18e LHC18f LHC18g LHC18h LHC18i LHC18j LHC18k LHC18l LHC18m LHC18n LHC18o LHC18p"

for name in $names
  do
  root -b -q 'clusterQA.C("'$name'", 1)'
  echo '"'$name'"'
done
