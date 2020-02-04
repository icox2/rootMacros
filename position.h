#include <cmath>

pair <double, double> position(double xa, double xb, double ya, double yb){

	double xpos = 25 * ((xa+xb)-(ya+yb))/(xa+xb+ya+yb);
	double ypos = 25 * ((xa+yb)-(xb+ya))/(xa+xb+ya+yb);

	return std::make_pair (xpos,ypos);
}
