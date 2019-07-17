#!/bin/bash

VER=v1.0.dev

PKG=gconv-$VER-$1

rm -rf $PKG
mkdir $PKG

make distclean

if [ $1 == "glnx64" ]; then

    g++ src/*.cpp -o $PKG/gconv -s -O2 -std=c++11 -static
    qmake-qt5 src/gui
    make
    strip gconv-gui
    mv gconv-gui $PKG/

elif [ $1 == "win32" ]; then

    i686-w64-mingw32-g++ src/*.cpp -o $PKG/gconv.exe -s -O2 -std=c++11 -static
    i686-w64-mingw32-qmake-qt4 src/gui
    make release
    i686-w64-mingw32-strip release/gconv-gui.exe
    mv release/gconv-gui.exe $PKG/

elif [ $1 == "win64" ]; then

    x86_64-w64-mingw32-g++ src/*.cpp -o $PKG/gconv.exe -s -O2 -std=c++11 -static
    x86_64-w64-mingw32-qmake-qt4 src/gui
    make release
    x86_64-w64-mingw32-strip release/gconv-gui.exe
    mv release/gconv-gui.exe $PKG/

elif [ $1 == "macos" ]; then

    export PATH="/usr/local/opt/qt/bin:$PATH"
    export LDFLAGS="-L/usr/local/opt/qt/lib"
    export CPPFLAGS="-I/usr/local/opt/qt/include"

    g++ src/*.cpp -o $PKG/gconv -O2 -std=c++11
    qmake src/gui
    make
    macdeployqt gconv-gui.app
    mv gconv-gui.app $PKG/

fi

if [[ $1 == win* ]]; then
    zip -qr $PKG.zip $PKG
else
    tar zcf $PKG.tar.gz $PKG
fi
