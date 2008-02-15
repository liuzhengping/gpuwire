#include <QtGui>
#include "Window.h"
#include "GLWidget.h"


Window::Window()
{
	glWidget = new GLWidget(this);
	QVBoxLayout *mainLayout = new QVBoxLayout;
	mainLayout->addWidget(glWidget);		
	
	
	
	labelX = new QLabel("X:");
	labelY = new QLabel("Y:");	
		
	//mainLayout->addWidget(labelX);
	//mainLayout->addWidget(labelY);
	
	coordWidget = new CoordWidget();
	
	mainLayout->addWidget(coordWidget);
	
	setLayout(mainLayout);
	setWindowTitle("GPU LiveWire");
	setMinimumSize(QSize(0, 0));
	
	QObject::connect(glWidget,SIGNAL(coordChanged(int,int)), coordWidget, SLOT(setCoord(int,int)));
}

void Window::setX(int x){
	labelX->setText(tr("X: %1").arg(x));
}
void Window::setY(int y){
	labelY->setText(tr("Y: %1").arg(y));
}

Window::~Window()
{
}
