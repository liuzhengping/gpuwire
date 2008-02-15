#ifndef GLWIDGET_H_
#define GLWIDGET_H_

//#include<QWidget>

#include<GL/glut.h>
	

#include<QGLWidget>
#include<QWheelEvent>
#include "Window.h"

#define ORIGINALTEX 0
#define SOBELTEX 1

//extern void APIENTRY glBlendColor (GLclampf, GLclampf, GLclampf, GLclampf);

//Isolated definition
typedef struct{
    unsigned char x, y, z, w;
} uchar4;

////////////////////////////////////////////////////////////////////////////////
// Small BMP loading utility
////////////////////////////////////////////////////////////////////////////////
extern "C" void LoadBMPFile(uchar4 **, int *, int *, const char *);

//CUDA kernel wrapper
void 
runTestGraph( int argc, char** argv,int iw, int ih, int startNode, float* gradient);
void printPathSource(int destination,int* path);

void empty();

class GLWidget : public QGLWidget
{
	Q_OBJECT

public:
	GLWidget(QWidget* parent=0);
	virtual ~GLWidget();
	QSize minimumSizeHint() const;
	QSize maximumSizeHint() const;
    QSize sizeHint() const;
    
protected:
     void initializeGL();
     void paintGL();	
     void resizeGL(int width, int height);
     void mousePressEvent(QMouseEvent *event);
     void mouseMoveEvent(QMouseEvent *event);
     void wheelEvent ( QWheelEvent * event );
     void keyPressEvent(QKeyEvent* event);  
     float zoomFactor,zoomFactorX,zoomFactorY;
     float ry;
     float dx,dy;
     int mouseX, mouseY;
     //float wWidth,wHeight;
     QPoint lastPos;
     Window* myWindow;
     GLuint tex, mask, sobel;
     uchar4* h_Src;
     int imageW, imageH;
     float * gradient;
     int chosenTexture;
     int drawTeapot;
     float approach;
     
     
signals: 
	void coordChanged(int x,int y);    
     
private:
	bool calculatedGradient;
	
};


#endif /*GLWIDGET_H_*/
