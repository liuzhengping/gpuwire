#ifndef COORDWIDGET_H_
#define COORDWIDGET_H_

#include <QLabel>
#include <QVBoxLayout>
#include <string>

class CoordWidget : public QWidget
{
	Q_OBJECT
public:
	CoordWidget();
	virtual ~CoordWidget();
private:
	QLabel* coordLabel;
public slots:
	void setCoord(int x, int y);
	
};

#endif /*COORDWIDGET_H_*/
