cmake_minimum_required(VERSION 3.26)

project(Curves)

set(CMAKE_EXPORT_COMPILE_COMMANDS True)

add_compile_options(-O3)

set(SRC 
        src/Merson.cpp 
        src/Euler.cpp 
        src/CurveFlow.cpp 
        src/CurveFlowForce.cpp
        src/DeTurck.cpp 
        src/functions.cpp 
        src/RedistributionInSpace.cpp 
        src/RedisForce.cpp
        src/AreaPreserving.cpp
)

add_executable(curves ${SRC} src/curve.cpp)

include_directories(
    "./include/"
)
