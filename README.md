# rinex_ho

# Google Colaboratory example:
[Open notebook on Google Colab](https://colab.research.google.com/drive/1RpyO_5xuy44qs30U3xIRyYTYpcSjSvWu?usp=sharing)

## Requirements
git, cmake, g++/mingw

### On Ubuntu
  sudo apt-get install g++ cmake build-essential libncurses-dev wget unzip zip gzip

### On windows, install these packages: 
  https://cmake.org/download/\
  It's recommended to install a C++ IDE with a compiler. Ex.: [Qt Creator](https://www.qt.io/download-qt-installer)\
  It's recommended to install a Git client. Ex.: [GitCola](https://git-cola.github.io/downloads.html), [GitHub Desktop](https://desktop.github.com/)

## Download from GitHub
git clone git@github.com:mauriciodev/rinex_ho.git

## Or download and unzip: 
https://github.com/mauriciodev/rinex_ho/archive/refs/heads/main.zip

## Installation from source
cd rinex_ho\
mkdir bin\
cd bin\
cmake ..\
make -j4

### On windows (using Git-cola and QtCreator):
Download the source via Git-cola. \
Open QtCreator and open the CMakeLists.txt from Rinex_Ho source folder. \
Choose "Desktop Qt MinGW" on the Projets tab, under Build & Run. \
Go to the "Edit" tab, right click on Rinex_Ho project, press "Run CMake". This should prepare the Makefile. \
On the same menu, press "Build". \
To debug, use the Debug menu.



## Testing (on the bin directory)
./rinex_ho



