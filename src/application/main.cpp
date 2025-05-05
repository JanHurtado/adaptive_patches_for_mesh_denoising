
#include <gui/myMainWindow.h>
#include "test.h"

#include <string>
#include <fstream>
#include <streambuf>


/** @defgroup mesh_processing Mesh Processing
*   Module containing mesh processing tools useful for the application
*/

/** @defgroup visualization Visualization based on OpenGL
*   Module containing OpenGL abstractions for visualization of triangular meshes
*/

int main(int argc, char **argv)
{
	QApplication app(argc, argv);
	myMainWindow main_window;
	main_window.show();
	/*myGLWindow myWindow;
	myWindow.show();*/
	printf("OpenGL version supported by this platform (%s): \n",glGetString(GL_VERSION));
	//initializeTestCases();
	//runTests();
	//comparison1();
	//comparison3();
	//latexCodeGenerator("results_devil.txt", "results_devil_new.txt");
	//latexCodeGenerator("results_block.txt", "results_block_new.txt");
	//latexCodeGenerator("results_fandisk.txt", "results_fandisk_new.txt");
	//latexCodeGenerator("results_joint.txt", "results_joint_new.txt");
	//latexCodeGenerator("results_sharpSphere.txt", "results_sharpSphere_new.txt");

	return app.exec();
}
