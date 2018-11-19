#!/bin/bash

VER=2019.0.dev

rm -rf gconv-$1
mkdir gconv-$1
make distclean

if [ $1 == "lnx64" ]; then

    g++ *.cpp -o gconv-$1/gconv -s -O2 -std=c++11 -static
    qmake-qt4 gui
    make
    strip gconv-gui
    mv gconv-gui gconv-$1/

elif [ $1 == "win32" ]; then

    i686-w64-mingw32-g++ *.cpp -o gconv-$1/gconv.exe -s -O2 -std=c++11 -static
    i686-w64-mingw32-qmake-qt4 gui
    make release
    i686-w64-mingw32-strip release/gconv-gui.exe
    mv release/gconv-gui.exe gconv-$1/

elif [ $1 == "win64" ]; then

    x86_64-w64-mingw32-g++ *.cpp -o gconv-$1/gconv.exe -s -O2 -std=c++11 -static
    x86_64-w64-mingw32-qmake-qt4 gui
    make release
    x86_64-w64-mingw32-strip release/gconv-gui.exe
    mv release/gconv-gui.exe gconv-$1/

elif [ $1 == "macos" ]; then

    export PATH="/usr/local/opt/qt/bin:$PATH"
    export LDFLAGS="-L/usr/local/opt/qt/lib"
    export CPPFLAGS="-I/usr/local/opt/qt/include"

    g++ *.cpp -o gconv -O2 -std=c++11
    qmake gui
    make
    macdeployqt gconv-gui.app
    mv gconv gconv-gui.app/Contents/MacOS/
    mv gconv-gui.app gconv-$1/

fi

tar zcf gconv-${VER}-${1}.tar.gz gconv-$1
