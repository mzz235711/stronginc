#!/bin/bash -e

echo "make sure you are in the build dir."

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/BIGDATA1/buaa_wffan_1/local/lib64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/BIGDATA1/buaa_wffan_1/local/lib

CXX=g++ CC=gcc cmake .. -DCMAKE_BUILD_TYPE=Release -DTIANHE=true -DGRAPE_EPOD=true -DGRAPE_VPOD=true -DBOOST_ROOT="/BIGDATA1/buaa_wffan_1/local" -DFOLLY_ROOT_DIR="/BIGDATA1/buaa_wffan_1/local" -Dgflags_DIR=/BIGDATA1/buaa_wffan_1/local -Dglog_DIR=/BIGDATA1/buaa_wffan_1/local -DGFLAGS_INCLUDE_DIR=/BIGDATA1/buaa_wffan_1/local/include -DGFLAGS_LIBRARY=/BIGDATA1/buaa_wffan_1/local/lib/libgflags.so -DGLOG_INCLUDE_DIR=/BIGDATA1/buaa_wffan_1/local/include -DGLOG_LIBRARY=/BIGDATA1/buaa_wffan_1/local/lib64/libglog.so -DFOLLY_INCLUDE_DIR=/BIGDATA1/buaa_wffan_1/local/include
make -j24

