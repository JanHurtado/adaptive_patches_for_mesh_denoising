project(AdaptivePatchesForMeshDenoising)
cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

cmake_policy(SET CMP0053 NEW)
cmake_policy(SET CMP0020 NEW)

set(CMAKE_CXX_STANDARD 14)

add_definitions(-D_USE_MATH_DEFINES)

add_subdirectory(core)
add_subdirectory(visualization)
add_subdirectory(gui)
#add_subdirectory(application)

link_directories("${PROJECT_SOURCE_DIR}/../dependencies/OpenMesh/lib" "${PROJECT_SOURCE_DIR}/../dependencies/cplex/lib" "${PROJECT_SOURCE_DIR}/../dependencies/glew/lib")

add_executable(main application/main.cpp application/test.h)

target_link_libraries(
    main
	core
	debug OpenMeshCored
	optimized OpenMeshCore
	optimized cplex1280
	debug cplex1280d
	visualization
	gui
)

install( TARGETS main DESTINATION bin )

file( GLOB QT5_DLLS ../dependencies/Qt5/bin/*.dll )
file( GLOB QT5_PLUGINS_DLLS ../dependencies/Qt5/plugins/* )
file( GLOB OPENMESH_DLLS ../dependencies/OpenMesh/bin/*.dll )
file( GLOB GLEW_DLLS ../dependencies/glew/bin/*.dll )
file( GLOB SHADERS ../shaders/*.glsl )

install( FILES ${QT5_DLLS} DESTINATION bin )
install( DIRECTORY ${QT5_PLUGINS_DLLS} DESTINATION bin )
install( FILES ${OPENMESH_DLLS} DESTINATION bin )
install( FILES ${GLEW_DLLS} DESTINATION bin )
install( FILES ${SHADERS} DESTINATION bin )