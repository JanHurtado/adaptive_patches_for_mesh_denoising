#include "gui/myGLWindow.h"

myGLWindow::myGLWindow(QWidget *parent)
: QOpenGLWidget(parent)
{
}


myGLWindow::~myGLWindow()
{
	delete renderer;
}

void myGLWindow::sendDataToOpenGL()
{
	makeCurrent();
	renderer->sendDataSingleBuffer();
	doneCurrent();
}

void myGLWindow::initializeGL()
{
	renderer = new myRenderer;
	setMouseTracking(true);
	renderer->addShader(GL_VERTEX_SHADER, "VertexShaderCodePhong.glsl");
	renderer->addShader(GL_FRAGMENT_SHADER, "FragmentShaderCodePhong.glsl");
	renderer->initialize();
}

void myGLWindow::paintGL()
{
    renderer->setWidth(width());
    renderer->setHeight(height());
	renderer->draw();
}

bool myGLWindow::event(QEvent *event)
{
	if (event->type() == QEvent::MouseMove) {
		QMouseEvent * e = static_cast<QMouseEvent *>(event);
		if (e->buttons() == Qt::LeftButton)
			renderer->rotateObjects(glm::vec2(e->x(), e->y()));
		if (e->buttons() == Qt::RightButton)
			renderer->translateCamera(glm::vec2(e->x(), e->y()));
		repaint();
		return true;
	}
	else if (event->type() == QEvent::Wheel)
	{
		QWheelEvent * e = static_cast<QWheelEvent *>(event);
		renderer->zoom(static_cast<float>(e->delta()));
		repaint();
		return true;
	}
	return QOpenGLWidget::event(event);
}

void myGLWindow::setVisualizationMode(int mode)
{ 
	if (renderer->getNumberOfShapes() > 0)
	{
		renderer->setShapeDrawMode(0, (myDrawFlags)mode);
		repaint();
	}
}

void myGLWindow::setShape(ShapeData * _shape)
{
	renderer->clearShapes();
	renderer->addShape(_shape);
	makeCurrent();
	renderer->resendDataSingleBuffer();
	doneCurrent();
	renderer->setDefaultValues();
	repaint();
}

void myGLWindow::removeSelection()
{
	if (renderer->getNumberOfShapes() > 1)
	{
		renderer->removeShape(1);
		makeCurrent();
		renderer->resendDataSingleBuffer();
		doneCurrent();
		repaint();
	}
	else
	{
		repaint();
	}
}

void myGLWindow::setSelection(ShapeData * _selection)
{
	if (renderer->getNumberOfShapes() <= 1)
	{
		renderer->addShape(_selection, e_draw_selection);
	}
	else
	{
		renderer->removeShape(1);
		renderer->addShape(_selection,e_draw_selection);
	}
	makeCurrent();
	renderer->resendDataSingleBuffer();
	doneCurrent();
	repaint();
}

void myGLWindow::updateMesh()
{
	if (renderer->getNumberOfShapes() > 0)
	{
		makeCurrent();
		renderer->updateVertexBuffer(0);
		doneCurrent();
		repaint();
	}
}

void myGLWindow::addShader(GLenum _shaderType, const string & _fileName)
{
	makeCurrent();
	renderer->addShader(_shaderType,_fileName);
	doneCurrent();
}

void myGLWindow::clearAndDeleteShaders()
{
	makeCurrent();
	renderer->clearAndDeleteShaders();
	doneCurrent();
}

void myGLWindow::installShaders()
{
	makeCurrent();
	renderer->installShaders();
	doneCurrent();
}

glm::vec2 myGLWindow::getCurrentMousePosition()
{ 
	QPoint p = mapFromGlobal(QCursor::pos()); return glm::vec2(p.x(), p.y());
}

glm::vec3 myGLWindow::getRayDirection(glm::vec2 & pos)
{ 
	return renderer->getRayDirection(pos);
}

myCamera * myGLWindow::getCamera()
{ 
    return renderer->getCamera();
}

glm::mat4 myGLWindow::getModelToWorldMatrix()
{ 
    return renderer->getModelToWorldMatrix();
}

void myGLWindow::screenShot()
{
	int scale = 2;
	resize(width()*scale, height()*scale);
	makeCurrent();
	glReadBuffer(GL_BACK);
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	int w = renderer->getWidth();
	int h = renderer->getHeight();
	int nSize = w*h * 3;
	// First let's create our buffer, 3 channels per Pixel
	char* dataBuffer = (char*)malloc(nSize*sizeof(char));

	if (!dataBuffer) return;

	// Let's fetch them from the backbuffer	
	// We request the pixels in GL_BGR format, thanks to Berzeger for the tip
	glReadPixels((GLint)0, (GLint)0, (GLint)w, (GLint)h, GL_RGB, GL_UNSIGNED_BYTE, dataBuffer);

	CImg<unsigned char> image(w, h, 1, 3);


	int pos = 0;
	for (int r = h - 1; r >= 0; r--)
	{
		for (int c = 0; c < w; c++)
		{

			image(c, r, 0, 0) = dataBuffer[pos];
			pos++;
			image(c, r, 0, 1) = dataBuffer[pos];
			pos++;
			image(c, r, 0, 2) = dataBuffer[pos];
			pos++;
		}
	}


	free(dataBuffer);

	image.save("ss.bmp");
	doneCurrent();
	resize(width()/scale, height()/scale);
}