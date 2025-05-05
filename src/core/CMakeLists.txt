project(core)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(OpenMP)

if(OPENMP_FOUND)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#include_directories ("${PROJECT_SOURCE_DIR}/../../dependencies/OpenMesh/include" "${PROJECT_SOURCE_DIR}/../../dependencies/Eigen/include" "${PROJECT_SOURCE_DIR}/../../dependencies/cplex/include")

#link_directories("${PROJECT_SOURCE_DIR}/../../dependencies/OpenMesh/lib" "${PROJECT_SOURCE_DIR}/../../dependencies/cplex/lib")

set(CORE_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include "${PROJECT_SOURCE_DIR}/../../dependencies/OpenMesh/include" "${PROJECT_SOURCE_DIR}/../../dependencies/Eigen/include" "${PROJECT_SOURCE_DIR}/../../dependencies/cplex/include")

set(CORE_LINK_DIRS "${PROJECT_SOURCE_DIR}/../../dependencies/OpenMesh/lib" "${PROJECT_SOURCE_DIR}/../../dependencies/cplex/lib")

link_directories(${CORE_LINK_DIRS})

add_definitions(-D_USE_MATH_DEFINES)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY AUTOGEN_TARGETS_FOLDER AutoMoc)

file(GLOB_RECURSE CORE_SOURCES *.cpp *.h)

add_library(${PROJECT_NAME} ${CORE_SOURCES})

target_include_directories( ${PROJECT_NAME}
    PUBLIC ${CORE_INCLUDE_DIRS}
)


target_link_libraries(
    ${PROJECT_NAME}
	debug OpenMeshCored
	optimized OpenMeshCore
	optimized cplex1280
	debug cplex1280d
)

#set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT main)


