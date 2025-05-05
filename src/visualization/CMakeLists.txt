project(visualization)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(OpenGL REQUIRED)

set(VISUALIZATION_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include "${PROJECT_SOURCE_DIR}/../../dependencies/glew/include" "${PROJECT_SOURCE_DIR}/../../dependencies/glm" "${PROJECT_SOURCE_DIR}/../../dependencies/OpenMesh/include")

set(VISUALIZATION_LINK_DIRS "${PROJECT_SOURCE_DIR}/../../dependencies/glew/lib" "${PROJECT_SOURCE_DIR}/../../dependencies/OpenMesh/lib")

link_directories(${VISUALIZATION_LINK_DIRS})

add_definitions(-D_USE_MATH_DEFINES)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY AUTOGEN_TARGETS_FOLDER AutoMoc)

file(GLOB_RECURSE VISUALIZATION_SOURCES *.cpp *.h)

add_library(${PROJECT_NAME} ${VISUALIZATION_SOURCES})

target_include_directories( ${PROJECT_NAME}
    PUBLIC ${VISUALIZATION_INCLUDE_DIRS}
)

target_link_libraries(
    ${PROJECT_NAME}
	debug OpenMeshCored
	optimized OpenMeshCore
	${OPENGL_gl_LIBRARY}
	glew32
	glew32s
	core
)

