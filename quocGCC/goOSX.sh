#!/bin/bash
cmake -GXcode -DCMAKE_BUILD_TYPE=Release -DUSE_PNG=1 -DPARSE_GCC_ERRORS=0 -DUSE_C++11=1 -DUSE_TIFF=1 ../quocmesh
