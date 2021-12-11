cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include("cmake/HunterGate.cmake")

HunterGate(
        URL "https://github.com/cpp-pm/hunter/archive/v0.23.297.tar.gz"
        SHA1 "3319fe6a3b08090df7df98dee75134d68e2ef5a3"
)

project(Game_theory)

hunter_add_package(nlohmann_json)
find_package(nlohmann_json 3.8.0 REQUIRED)

add_executable(${PROJECT_NAME}
        sources/main.cpp
        sources/simplex.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
        "$<INSTALL_INTERFACE:include>"
        )

target_link_libraries(${PROJECT_NAME} nlohmann_json::nlohmann_json)
