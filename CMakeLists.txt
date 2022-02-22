cmake_minimum_required(VERSION 3.5)
project(Rinex_ho2 LANGUAGES CXX C)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_EXE_LINKER_FLAGS " -static")

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

# set everything up for c++ 17 features
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


#find_library(GFORTRAN_LIBRARY NAMES libgfortran.so)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/rinex_ha.inp DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/rinex_ha_param.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
#file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/CGLALO.DAT DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/BRFT DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/IGRF_COEF DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/RECF DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/ASC1 DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)



#find_package(Qt5 COMPONENTS Widgets LinguistTools REQUIRED)

#set(TS_FILES Rinex_ho2/Rinex_ho2_pt_BR.ts)
#set(UI_SOURCES Rinex_ho2/main.cpp
#     Rinex_ho2/mainwindow.cpp
#     Rinex_ho2/mainwindow.h
#     Rinex_ho2/mainwindow.ui)
# 
# add_executable(Rinex_ho2
#     ${UI_SOURCES}
#     ${TS_FILES}
#   )

# target_link_libraries(Rinex_ho2 PRIVATE Qt5::Widgets)



set(RINEXHO_HEADERS
    Comp_Bias.h
    Class_Iono.h
    Class_TEC.h
    CycleSlip.h
    Filter_Code.h
    Class_IGRF11.h
    rinex.h
    datetime.h
    rinex_ho.h
    Variable.h
)

set(RINEXHO_SOURCES
    Comp_Bias.cpp
    Class_Iono.cpp
    Class_TEC.cpp
    CycleSlip.cpp
    Filter_Code.cpp
    Class_IGRF11.cpp
    rinex.cpp
    datetime.cpp
    rinex_ho.cpp
)

add_executable(Rinex_ho
    ${RINEXHO_SOURCES}
    ${RINEXHO_HEADERS}
  )
#target_link_libraries(Rinex_ho stdc++fs)
#target_link_libraries(Rinex_ho cgm_util)
#qt5_create_translation(QM_FILES ${CMAKE_SOURCE_DIR} ${TS_FILES})
