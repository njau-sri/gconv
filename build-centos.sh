#!/bin/bash

rm -rf gconv-$1
mkdir gconv-$1
make distclean

if [ $1 == "lnx64" ]; then

    g++ *.cpp -o gconv-$1/gconv -s -O2 -std=c++11 -static
    qmake-qt4 gui
    make
    strip -s gconv-gui
    mv gconv-gui gconv-$1/

elif [ $1 == "win32" ]; then

    i686-w64-mingw32-g++ *.cpp -o gconv-$1/gconv.exe -s -O2 -std=c++11 -static
    i686-w64-mingw32-qmake-qt4 gui
    make
    i686-w64-mingw32-strip -s release/gconv-gui.exe
    mv release/gconv-gui.exe gconv-$1/

elif [ $1 == "win64" ]; then

    x86_64-w64-mingw32-g++ *.cpp -o gconv-$1/gconv.exe -s -O2 -std=c++11 -static
    x86_64-w64-mingw32-qmake-qt4 gui
    make
    x86_64-w64-mingw32-strip -s release/gconv-gui.exe
    mv release/gconv-gui.exe gconv-$1/

fi

tar zcf gconv-$1.tar.gz gconv-$1
