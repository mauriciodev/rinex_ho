cmake_minimum_required(VERSION 3.5)
project(Rinex_ho2 LANGUAGES CXX C Fortran)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


#find_library(GFORTRAN_LIBRARY NAMES libgfortran.so)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/rinex_ha.inp DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/rinex_ha_param.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/CGLALO.DAT DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/BRFT DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/IGRF11 DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)


find_package(Qt5 COMPONENTS Widgets LinguistTools REQUIRED)

set(TS_FILES Rinex_ho2/Rinex_ho2_pt_BR.ts)
set(UI_SOURCES Rinex_ho2/main.cpp
    Rinex_ho2/mainwindow.cpp
    Rinex_ho2/mainwindow.h
    Rinex_ho2/mainwindow.ui)

add_executable(Rinex_ho2
    ${UI_SOURCES}
    ${TS_FILES}
  )

target_link_libraries(Rinex_ho2 PRIVATE Qt5::Widgets)



set(FORTRAN_SOURCES
    CGM_DLL_Lahey_Fujtsu/cgm_util.f
    CGM_DLL_Lahey_Fujtsu/precip.f
    )
add_library(cgm_util STATIC ${FORTRAN_SOURCES})

set(RINEXHO_HEADERS
    Class_Iono/Class_Iono.h
    Class_TEC/Class_TEC.h
    Cycle_Slip/CycleSlip.h
    FilterCode/Filter_Code.h
    IGRF11/Class_IGRF11.h
    Rinex_Class/rinex.h
    Rinex_Class/datetime.h
    rinex_ho.h
)

set(RINEXHO_SOURCES
    Class_Iono/Class_Iono.cpp
    Class_TEC/Class_TEC.cpp
    Cycle_Slip/CycleSlip.cpp
    FilterCode/Filter_Code.cpp
    IGRF11/Class_IGRF11.cpp
    Rinex_Class/rinex.cpp
    Rinex_Class/datetime.cpp
    rinex_ho.cpp
)

add_executable(Rinex_ho
    ${RINEXHO_SOURCES}
    ${RINEXHO_HEADERS}
  )
target_link_libraries(Rinex_ho cgm_util)
qt5_create_translation(QM_FILES ${CMAKE_SOURCE_DIR} ${TS_FILES})
