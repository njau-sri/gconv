#!/bin/bash

rm -rf gconv-$1
mkdir gconv-$1

if [ $1 == "lnx64" ]; then

    g++ *.cpp -o gconv-$1/gconv -s -O2 -std=c++11 -static

elif [ $1 == "win32" ]; then

    i686-w64-mingw32-g++ *.cpp -o gconv-$1/gconv.exe -s -O2 -std=c++11 -static

elif [ $1 == "win64" ]; then

    x86_64-w64-mingw32-g++ *.cpp -o gconv-$1/gconv.exe -s -O2 -std=c++11 -static

fi

tar zcf gconv-$1.tar.gz gconv-$1
