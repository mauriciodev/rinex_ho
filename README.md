# rinex_ho

## Requirements
git, cmake, g++/mingw

### On Ubuntu
  sudo apt-get install g++ cmake build-essential libncurses-dev wget unzip zip gzip

### On windows, install these packages: 

  https://cmake.org/download/
  
  It's recommended to install a C++ IDE with a compiler. Ex.: [Qt Creator](https://www.qt.io/download-qt-installer)
  
  It's recommended to install a Git client. Ex.: [GitCola](https://git-cola.github.io/downloads.html), [GitHub Desktop](https://desktop.github.com/)

## Download from GitHub
git clone git@github.com:mauriciodev/rinex_ho.git

## Or download and unzip: 
https://github.com/mauriciodev/rinex_ho/archive/refs/heads/main.zip

## Installation from source
cd rinex_ho

mkdir bin

cd bin

cmake ..

make -j4

## Testing (on the bin directory)
./rinex_ho



