
#include <QtGui>
#include <QtOpenGL>
#include <math.h>
#include <stdlib.h>




#include "GLWidget.h"


void filter(uchar4 ** h_Src,int width, int height, float * gradient){		
	uchar4 temp[width*height];
	for(int i=0;i<width*height;i++){
		temp[i].x = (*h_Src)[i].x;
		temp[i].y = (*h_Src)[i].y;
		temp[i].z = (*h_Src)[i].z;
		temp[i].w = (*h_Src)[i].w;
	}
	for(int i=width+1;i<height*width-width-1;i++){
		if((i%width==0)||(i%width==width-1)) continue;
		
		short ulx = temp[i+width-1].x;
		short upx = temp[i+width].x;
		short urx = temp[i+width+1].x;
		
		short clx = temp[i-1].x;
		short cex = temp[i].x;
		short crx = temp[i+1].x;
		
		short dlx = temp[i-width-1].x;
		short downx = temp[i-width].x;
		short drx = temp[i-width+1].x;
		
		short Horzx = urx + 2*crx + drx - ulx - 2*clx - dlx;
    	short Vertx = ulx + 2*upx + urx - dlx - 2*downx - drx;
    	Horzx /= 8;
    	Vertx /= 8;
    	short Sumx = (short) sqrt( Horzx * Horzx + Vertx * Vertx );
    	unsigned char chanX;    	
	    if ( Sumx < 0 ) chanX= 0; 
	    else if ( Sumx > 0xff ) chanX= 0xff;
    	else chanX = (unsigned char) Sumx;
		
		short uly = temp[i+width-1].y;
		short upy = temp[i+width].y;
		short ury = temp[i+width+1].y;
		
		short cly = temp[i-1].y;
		short cey = temp[i].y;
		short cry = temp[i+1].y;
		
		short dly = temp[i-width-1].y;
		short downy = temp[i-width].y;
		short dry = temp[i-width+1].y;
		
		short Horzy = ury + 2*cry + dry - uly - 2*cly - dly;
    	short Verty = uly + 2*upy + ury - dly - 2*downy - dry;
    	Horzy /= 8;
    	Verty /= 8;
    	short Sumy = (short) sqrt( Horzy*Horzy + Verty*Verty);
    	unsigned char chanY;    	
	    if ( Sumy < 0 ) chanY= 0; 
	    else if ( Sumy > 0xff ) chanY= 0xff;
    	else chanY = (unsigned char) Sumy;
    	
    	short ulz = temp[i+width-1].z;
		short upz = temp[i+width].z;
		short urz = temp[i+width+1].z;
		
		short clz = temp[i-1].z;
		short cez = temp[i].z;
		short crz = temp[i+1].z;
		
		short dlz = temp[i-width-1].z;
		short downz = temp[i-width].z;
		short drz = temp[i-width+1].z;
		
		short Horzz = urz + 2*crz + drz - ulz - 2*clz - dlz;
    	short Vertz = ulz + 2*upz + urz - dlz - 2*downz - drz;
    	Horzz /= 8;
    	Vertz /= 8;
    	short Sumz = (short) sqrt ( Horzz*Horzz + Vertz*Vertz );
    	unsigned char chanZ;    	
	    if ( Sumz < 0 ) chanZ= 0; 
	    else if ( Sumz > 0xff ) chanZ= 0xff;
    	else chanZ = (unsigned char) Sumz;
		
		
		(*h_Src)[i].x = chanX;
		(*h_Src)[i].y = chanY;
		(*h_Src)[i].z = chanZ;
		gradient[i] = 0.3 * chanX/255.0 + 0.59 * chanY/255.0 + 0.11 * chanZ/255.0;
		
		//(*h_Src)[i].y = std::max(temp[i].y - 50,0);
		//(*h_Src)[i].z = std::max(temp[i].z - 50,0);
		//(*h_Src)[i].w = temp[i].w ;
	}

}

GLWidget::GLWidget(QWidget* parent) : QGLWidget(parent)
{
	zoomFactor = 1.0f;
	zoomFactorX = 1.0f;
	zoomFactorY = 1.0f;
	ry=0.0f;	
	dx=0.0f;
	setFocusPolicy (Qt::StrongFocus);
	setMouseTracking(true);
	myWindow = (Window*)parent;
	calculatedGradient=false;
	chosenTexture = ORIGINALTEX;
	drawTeapot = 0;
	approach=0.0f;
		
}

GLWidget::~GLWidget()
{
}

 QSize GLWidget::minimumSizeHint() const
 {
     return QSize(50, 50);
 }

 QSize GLWidget::sizeHint() const
 {
     return QSize(512, 384);
 }
  QSize GLWidget::maximumSizeHint() const
 {
     return QSize(2000, 2000);
 }
 void GLWidget::initializeGL()
 {
 	int argc=0;
 	char** test;
 	glutInit(&argc,test); 	
 	QColor trolltechPurple = QColor::fromCmykF(0.39, 0.39, 0.0, 0.0);
 	
 	glMatrixMode(GL_PROJECTION);
 	glLoadIdentity();
 	
 	//gluPerspective(45.0f, (GLfloat)400.0/(GLfloat)400.0,0.1f,100.0f);
 	//gluPerspective(zoomFactor, (GLfloat)400.0/(GLfloat)400.0,0.1f,100.0f);
 	
     qglClearColor(trolltechPurple.dark());

     
    //loading the texture
    LoadBMPFile(&h_Src, &imageW, &imageH, "plane.bmp");
    
    glGenTextures(1,&tex);
    glBindTexture(GL_TEXTURE_2D, tex);    
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);    
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, imageW, imageH, 0, GL_RGBA, GL_UNSIGNED_BYTE, h_Src);
    
    for(int i=0;i<imageW*imageH;i++){
    	h_Src[i].w=255;
    }
    
    
    gradient = (float*) malloc( imageW*imageH*sizeof(float));
    
    
    
    
    filter(&h_Src,imageW,imageH,gradient);
    
    glGenTextures(1,&sobel);
    glBindTexture(GL_TEXTURE_2D, sobel);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, imageW, imageH, 0, GL_RGBA, GL_UNSIGNED_BYTE, h_Src);
    
    /*for(int i=0;i < imageH;i++){
    	for(int j=0; j<imageW ; j++){
    		printf("%.3f ",gradient[i*imageW+j]);
    	}
    	printf("\n");
    }*/
    //float* dummy;
    
        
    printf("My width %d my height %d\n",imageW,imageH);
    printf("byte0 %d %d %d %d \n",h_Src[320].x,h_Src[320].y,h_Src[320].z,h_Src[320].w);
    printf("byte1 %d %d %d %d \n",h_Src[321].x,h_Src[321].y,h_Src[321].z,h_Src[321].w);    
    	
	glGenTextures(1,&mask);
	glBindTexture(GL_TEXTURE_2D, mask);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
	for(int i=0;i<imageW*imageH;i++){
		h_Src[i].x=0;
		h_Src[i].y=0;
		h_Src[i].z=0;
		h_Src[i].w=255;
	}
	
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, imageW, imageH, 0, GL_RGBA, GL_UNSIGNED_BYTE, h_Src);
	   
	
	 //object = makeObject();
	 glEnable(GL_TEXTURE_2D);
     glShadeModel(GL_SMOOTH);
     glEnable(GL_DEPTH_TEST);
     
     //glEnable(GL_CULL_FACE);
     //glEnable(GL_BACKFACE);
     glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
     glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
     
     glEnable(GL_BLEND);
     //glBlendColor(1.0f,1.0f,1.0f,1.0f);      
     glBlendFunc(GL_ONE,  GL_ONE);
     //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//     glBlendFunc(GL_ONE, GL_ONE);
     
 }
 
  void GLWidget::resizeGL(int width, int height)
 {
     int side = qMax(width, height);
     //glViewport((width - side) / 2, (height - side) / 2, side, side);
     glViewport( 0, 0, width, height);
     printf("My width %d height %d\n",width,height);

     glMatrixMode(GL_PROJECTION);
     glLoadIdentity();
     //glOrtho(-4.5, +4.5, +4.5, -4.5, +1.0, 350.0);
     glMatrixMode(GL_MODELVIEW);
 }

 void GLWidget::paintGL()
 {
 	glMatrixMode(GL_PROJECTION);
 	glLoadIdentity();
 	
 	
 	if(drawTeapot){
 		glEnable(GL_DEPTH_TEST); 		
 		gluPerspective(45.0f, 1.0,0.1f,100.0f);
 	}
 	else{
 		glDisable(GL_DEPTH_TEST);
 		glOrtho(0,1/zoomFactorX,0,1/zoomFactorY,-100,100);
 	}
 	//glOrtho(-1*zoomFactor,+1*zoomFactor,-1*zoomFactor,+1*zoomFactor,0,100);
 	//glOrtho(-zoomFactorX/2,zoomFactorX/2,-zoomFactorY/2,zoomFactorY/2,-100,100);
 	
 	
 	
 	//gluOrtho2D(-1*zoomFactor, 1*zoomFactor, -1*zoomFactor, 1*zoomFactor);  
	//gluPerspective(45.0f+distance, (GLfloat)400.0/(GLfloat)400.0,0.1f,100.0f);
	
 	glMatrixMode(GL_MODELVIEW);
     glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
     glLoadIdentity();
     
     //glScalef(1+.1*distance,1+.1*distance,1+.1*distance);
     if(drawTeapot){
     	glTranslated(dx/width()-0.5f, -dy/height()-0.5f, approach);//0.1*distance);
     }
     else{
     	glTranslated(dx/width(), -dy/height(), 0.0);//0.1*distance);
     }
          
     //glRotated(10*distance, 1.0, 0.0, 0.0);
     glRotated(ry, 0.0, 1.0, 0.0);
     glRotated(0.0, 0.0, 0.0, 1.0);
     
          

     
     //glEnable(GL_TEXTURE_2D);
     
glDisable(GL_BLEND);

	if(chosenTexture==ORIGINALTEX)     
		glBindTexture(GL_TEXTURE_2D,tex);
	else
		glBindTexture(GL_TEXTURE_2D,sobel);
 
 
glBegin( GL_QUADS );
glTexCoord2d(0.0,0.0); glVertex2d(0.0,0.0);
glTexCoord2d(1.0,0.0); glVertex2d(1.0,0.0);
glTexCoord2d(1.0,1.0); glVertex2d(1.0,1.0);
glTexCoord2d(0.0,1.0); glVertex2d(0.0,1.0);
glEnd();

glEnable(GL_BLEND);

     glBindTexture(GL_TEXTURE_2D,mask);
 
 
glBegin( GL_QUADS );
glTexCoord2d(0.0,0.0); glVertex2d(0.0,0.0);
glTexCoord2d(1.0,0.0); glVertex2d(1.0,0.0);
glTexCoord2d(1.0,1.0); glVertex2d(1.0,1.0);
glTexCoord2d(0.0,1.0); glVertex2d(0.0,1.0);
glEnd();
     




	
     
     
     if(drawTeapot){
     	glDisable(GL_BLEND);

	if(chosenTexture==SOBELTEX)     
		glBindTexture(GL_TEXTURE_2D,tex);
	else
		glBindTexture(GL_TEXTURE_2D,sobel);

		glBegin( GL_QUADS );
		glTexCoord2d(0.0,0.0); glVertex3d(0.0,0.0,-1.0);
		glTexCoord2d(1.0,0.0); glVertex3d(1.0,0.0,-1.0);
		glTexCoord2d(1.0,1.0); glVertex3d(1.0,1.0,-1.0);
		glTexCoord2d(0.0,1.0); glVertex3d(0.0,1.0,-1.0);
		glEnd();    	
     	
     	glBindTexture(GL_TEXTURE_2D,0);
     	glTranslated(0.0f,1.0f,-.5f);
     	glutWireTeapot(0.5f);
     }
 
     // Draw the top face
 /*    glBegin(GL_QUADS);
      //glTexCoord2f(1.0f,1.0f); 
     //glBegin(GL_TRIANGLES);
    	 
     	//glColor3f(0.0f,0.0f,1.0f);
     	glTexCoord2f(0.0f,0.0f);    
		glVertex3f(0.0f,0.0f,0.0f);		 
		
		glTexCoord2f(1.0f,0.0f);		
		glVertex3f(1.0f,0.0f,0.0f);		
		//glColor3f(0.0f,0.0f,0.0f);
		
		glTexCoord2f(1.0f,1.0f);
		glVertex3f(1.0f,1.0f,0.0f);
		
		
		glTexCoord2f(0.0f,1.0f);
		glVertex3f(0.0f,1.0f,0.0f);
     glEnd();
     */
     //printf("ZoomfactorX %f\n",zoomFactorX);
          
     
     
     //glCallList(object);
 }
 void GLWidget::wheelEvent ( QWheelEvent * event ){ 	
 	if(zoomFactorX+event->delta()/480.0>0)
 		zoomFactorX+=event->delta()/480.0;
 	if(zoomFactorY+event->delta()/480.0>0)
 		zoomFactorY+=event->delta()/480.0; 	
 	updateGL();
 }  
 
 void GLWidget::mousePressEvent(QMouseEvent *event){
 	int argc=0;
 	char** test;
 	lastPos = event->pos();
 	if(event->button() == Qt::LeftButton){
 		calculatedGradient = true;
 		int xpos = (int) ( -dx + event->x()/zoomFactorX + 0.5);
 		int ypos = (int) ( height() - dy  - (  (height()-event->y())/zoomFactorY )+ 0.5);
 		int startNode = (imageH-ypos)*imageW + xpos;
// 		int startNode = event->y()*imageW + event->x(); 	
 		runTestGraph( argc, test,imageW,imageH,startNode,gradient);
 	}
 	updateGL(); 	
 }
 void GLWidget::mouseMoveEvent(QMouseEvent *event){
 	int mydx = event->x() - lastPos.x();
 	int mydy = event->y() - lastPos.y();
 	int xpos = (int) ( -dx + event->x()/zoomFactorX + 0.5);
 	int ypos = (int) ( height() - dy  - (  (height()-event->y())/zoomFactorY )+ 0.5);
 	//printf("x %d y %d dy %d\n",xpos, ypos,dy);
 	myWindow->setX(xpos);
 	myWindow->setY(ypos);
 	emit coordChanged(xpos,ypos);
 	
 	mouseX = xpos;
 	mouseY = ypos;
 	
 	if(event->buttons() && Qt::LeftButton){ 
 		//printf("Delta x %f\n",dx);
 		dx += mydx / zoomFactorX;
 		dy += mydy / zoomFactorY;
 		//dx += zoomFactorX*mydx/width();
 		//dy += zoomFactorY*mydy/height();
 	}
 	else{ 		
 		if(calculatedGradient){
 			if( xpos>=0 && xpos < imageW && ypos >=0 && ypos < imageH){ 
		 		int node = (imageH-ypos)*imageW + xpos;
		 		int path[5000];
		 		for(int i=0;i<5000;i++)
		 			path[i]=-1;

		 		printPathSource(node,path);		 		
		 		glBindTexture(GL_TEXTURE_2D, mask);
				glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
				for(int i=0;i<imageW*imageH;i++){
					h_Src[i].x=0;
					h_Src[i].y=0;
					h_Src[i].z=0;
					h_Src[i].w=255;
				}			
	 		
		 		for(int i=0;i<1000;i++){
		 			if(path[i]==-1) break;
		 			h_Src[path[i]].x=255;
		 			h_Src[path[i]].y=255;
		 			h_Src[path[i]].z=0;
		 			printf("Received %d\n",path[i]);
		 		}
	 			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, imageW, imageH, 0, GL_RGBA, GL_UNSIGNED_BYTE, h_Src);
 			}
 		}
 	}
 	
 	lastPos = event->pos();
 	updateGL();
 	
 }
 
 void GLWidget::keyPressEvent(QKeyEvent* event){
 	if(event->key()==Qt::Key_R){
 		printf("Ry %f\n",ry);
 		ry+=10.0f;
 	}
 	else if(event->key()==Qt::Key_T){
 		printf("Ry %f\n",ry);
 		ry-=10.0f;
 	}
 	else if(event->key()==Qt::Key_O){
 		if(chosenTexture==ORIGINALTEX)
 			chosenTexture = SOBELTEX;
 		else
 			chosenTexture = ORIGINALTEX;
 	}
 	else if(event->key()==Qt::Key_P){
 		drawTeapot = !drawTeapot;
 	}
 	else if(event->key()==Qt::Key_W){
 		approach -= 0.1f;
 	}
 	else if(event->key()==Qt::Key_S){
 		approach += 0.1f;
 	}
 	else{
 		printf("Ignoring\n");
 		event->ignore();
 	} 
 	updateGL();
 		
 		
 }
