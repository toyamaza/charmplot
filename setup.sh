#!/bin/bash

LHAPDF_PATH=$(command -v lhapdf)
LHAPDF_PATH=${LHAPDF_PATH//\/bin\/lhapdf/}
export LHAPDF_DATA_PATH=/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current/:/cvmfs/sft.cern.ch/lcg/views/LCG_96bpython3/x86_64-centos7-gcc8-opt/share/LHAPDF
