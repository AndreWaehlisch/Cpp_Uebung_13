#include <cstdlib>
#include <limits>
#include <iostream>
#include "helfer.h" 

using namespace std; 

VecDoub a(const VecDoub& a, const VecDoub& b)
{
  int n = a.size(); 

  VecDoub sum(n); 
  
  for (int i = 0 ; i < n ; i++) 
    sum[i] = a[i] + b[i] ; 

  return sum; 
}

VecDoub m(const double& f,const VecDoub& a) 
{
  int n = a.size(); 
  VecDoub p(n); 
  for (int i = 0; i < n ; i++)
    p[i] = f*a[i]; 
  return p; 
}

template <class T>
void solver 
(VecDoub &y, const double t, const double h, T& ode)
{
  VecDoub k1(6),k2(6),k3(6),k4(6),aux1(6),aux2(6); 
  
  ode(t,y,k1); 
  ode(t+.5*h,a(y,m(.5*h,k1)),k2); 
  ode(t+.5*h,a(y,m(.5*h,k2)),k3); 
  ode(t+h,a(y,m(h,k3)),k4); 
  
  aux1 = m(1/3.,a(k2,k3)); 
  aux2 = m(1/6.,a(k1,k4)); 
  
  y = a(y,m(h,a(aux1,aux2))); 
}

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
  double t = 0; 
  double h = .01; 
  double a[] = {0.,0.,0.,1.,0.,0.,1.,3.,-1.}; 
  VecDoub y(6,a); 

  for(t = 0; t < 10; t+=h)
    {
      for (int i = 0; i < 6; i++) 
	cout << l.ystart[i] << " " ; 
      
      cout << endl; 
      solver<Lorentz>(l.ystart,t,h,l); 
    }   
}
