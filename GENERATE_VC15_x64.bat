set VTK_DIR=""
set ITK_DIR=""
set PATH=%CD%\dependencies\OpenMesh\bin;%CD%\dependencies\cplex\bin;%CD%\dependencies\glew\bin;%CD%\dependencies\Qt5\bin;%PATH%
set CMAKE_PREFIX_PATH=%CD%\dependencies\cplex;%CD%\dependencies\Qt5;
mkdir build
chdir build
cmake -G "Visual Studio 15 2017 Win64" -DCMAKE_INSTALL_PREFIX=install ../

@pause