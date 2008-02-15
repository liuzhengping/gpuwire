#ifndef WINDOW_H_
#define WINDOW_H_

#include <QWidget>
#include "CoordWidget.h"

class GLWidget;
class QLabel;

class Window : public QWidget
{
	Q_OBJECT
	
public:
	Window();
	virtual ~Window();
	void setX(int x);
	void setY(int y);
	
	
private:
	GLWidget *glWidget;
	QLabel *labelX;
	QLabel *labelY;
	CoordWidget *coordWidget;
	
};

#endif /*WINDOW_H_*/
