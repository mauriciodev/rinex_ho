# rinex_ho

## Requirements
git, cmake, g++/mingw

On Ubuntu: sudo apt-get install g++ cmake build-essential libncurses-dev wget unzip zip gzip

## Download from GitHub
git clone git@github.com:mauriciodev/rt_ppp.git

## Or download and unzip: 
https://github.com/mauriciodev/rt_ppp/archive/refs/heads/main.zip

unzip Binary.zip on rt_ppp folder

## Installation from source
cd rt_ppp

mkdir bin

cd bin

cmake ..

make -j4

## Testing (on the bin directory)
./rt_ppp OBS_DATA/2021/SJRP/10/sjrp0101.21o Config_file/RTPPP_Config_SJRP_2021.inp



## Installation with Docker 
### Install docker
https://docs.docker.com/engine/install/ubuntu/

### Building the container (run on rt_ppp folder)
docker build . -t rt_ppp

### Testing with Docker
docker run -v "$PWD"/OBS_DATA/2021/SJRP/10/sjrp0101.21o:/opt/input.xxo -v "$PWD"/OUT_DIR:/opt/rt_ppp/OUT_DIR/ -it rt_ppp ./rt_ppp /opt/input.xxo Config_file/RTPPP_Config_SJRP_2021.inp
