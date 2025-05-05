project(gui)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

#find_package(OpenGL REQUIRED)
find_package(Qt5Core)
find_package(Qt5Widgets)
find_package(Qt5Gui)
find_package(Qt5OpenGL)

set(GUI_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include "${PROJECT_SOURCE_DIR}/../../dependencies/glew/include" "${PROJECT_SOURCE_DIR}/../../dependencies/glm" "${PROJECT_SOURCE_DIR}/../../dependencies/OpenMesh/include")

set(GUI_LINK_DIRS "${PROJECT_SOURCE_DIR}/../../dependencies/glew/lib" "${PROJECT_SOURCE_DIR}/../../dependencies/OpenMesh/lib")

link_directories(${GUI_LINK_DIRS})

add_definitions(-D_USE_MATH_DEFINES)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY AUTOGEN_TARGETS_FOLDER AutoMoc)

file(GLOB_RECURSE GUI_SOURCES *.cpp *.h)

add_library(${PROJECT_NAME} ${GUI_SOURCES})

target_include_directories( ${PROJECT_NAME}
    PUBLIC ${GUI_INCLUDE_DIRS}
)

target_link_libraries(
    ${PROJECT_NAME}
	Qt5::Widgets
	core
	visualization
)

