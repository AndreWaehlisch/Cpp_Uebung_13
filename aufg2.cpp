#include <iostream>
#include "helfer.h" 
#include "stepperdopr853.h"
#include "stepperbs.h"

// Aus Aufgabe1:
class Lorentz{

public: 
  VecDoub B,ystart; 

  Lorentz(double x1, double x2, double x3, double v1, double v2, double v3, double b1, double b2, double b3) 
  {
    B.resize(3); 
    ystart.resize(6); 

      //Anfangspositionen 
      ystart[0] = x1; 
      ystart[1] = x2; 
      ystart[2] = x3; 
      //Anfangsgeschwindigkeiten 
      ystart[3] = v1; 
      ystart[4] = v2; 
      ystart[5] = v3; 
      B[0] = b1; 
      B[1] = b2; 
      B[2] = b3; 
    }

  void operator () ( const double t, const VecDoub& y, VecDoub& dydt) 
  {
    dydt[0] = y[3]; 
    dydt[1] = y[4]; 
    dydt[2] = y[5]; 
    dydt[3] = y[4]*B[2] - y[5]*B[1]; 
    dydt[4] = y[5]*B[0] - y[3]*B[2]; 
    dydt[5] = y[3]*B[1] - y[4]*B[0]; 
  }
};

int main()
{
	Lorentz l(0,0,0,1,0,0,1,3,-1);
	Output out(1000);
	Odeint<StepperDopr853 <Lorentz> > dgl(l.ystart, 0., 100, 1e-10, 1e-10, 0.01, 0, out, l);
	
	cout << out.steps[1] << endl;
}
