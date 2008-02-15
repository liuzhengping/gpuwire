#include "CoordWidget.h"

CoordWidget::CoordWidget()
{
	coordLabel = new QLabel("CoordLabel");
	QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(coordLabel);
    setLayout(layout);
}



CoordWidget::~CoordWidget()
{
}

void CoordWidget::setCoord(int x, int y){
	char temp[200];
	sprintf(temp,"%d,%d",x,y);
	QString str = temp;
	coordLabel->setText("Coord " + str);
	
}