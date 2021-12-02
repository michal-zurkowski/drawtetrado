#!/bin/bash

GCC_VERSION=$(gcc --version | grep ^gcc | sed 's/^.* //g')
GCC_MAJOR=$(echo $GCC_VERSION | cut -d. -f1)

echo "Detected GCC ${GCC_VERSION} with Major ${GCC_MAJOR}"

if [[ $GCC_MAJOR == 8 || $GCC_MAJOR == 9 ]]; then

  CXXFLAGS="-O3 -std=c++2a"

elif [[ $GCC_MAJOR -ge 10 ]]; then

  CXXFLAGS="-O3 -std=c++20" 

else
  echo "GCC version >= 8 is required for compilation"
  exit
fi

g++ -o svg_optimizer ${CXXFLAGS} svg_optimizer.cpp
