#include <cstdlib>
#include <limits>
#include <iostream>
#include <fstream>
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
  
  void reset (double x1, double x2, double x3, double v1, double v2, double v3, double b1, double b2, double b3) 
  {
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
  double tmax[] = {10,100,1e4};
  double h; 
  int N; 
  double a[] = {0.,0.,0.,1.,0.,0.,1.,3.,-1.}; 
  double vbetrag, veps, delta; 
  VecDoub y(6,a); 
 
  //ofstream file1("lorentz_1.txt",ios::trunc);  
  //ofstream file2("lorentz_2.txt",ios::trunc);  
  //ofstream file3("lorentz_3.txt",ios::trunc);  
  
  vbetrag = sqrt(l.ystart[3]*l.ystart[3] +  l.ystart[4]*l.ystart[4] +  l.ystart[5]*l.ystart[5]); 
  veps = .01 * vbetrag; 

 
  for ( int i = 0 ; i < 3; i++) 
    {
      N = 1; 
      delta = veps + 1;       
      while ( delta > veps ) 
	{
	  delta = veps + 1; 
	  h = tmax[i] / N; 
	  //	  cout << i << " "<< h << endl;
	  for(t = 0; t <= tmax[i]; t+=h)
	    solver<Lorentz>(l.ystart,t,h,l); 
	  
	  delta = abs(vbetrag - sqrt(l.ystart[3]*l.ystart[3] +  l.ystart[4]*l.ystart[4] +  l.ystart[5]*l.ystart[5]));
	  l.reset(0,0,0,1,0,0,1,3,-1); 
	  N++; 
	}
      cout << N << endl; 
    }
}

/*
  
      if ( t <= 100 )
	{
	  for (int i = 0; i < 6; i++) 
	    file2 << l.ystart[i] << " " ; 
	  file2 << endl; 
	}
      
      for (int i = 0; i < 6; i++) 
	file3 << l.ystart[i] << " " ; 
      file3 << endl; 
      solver<Lorentz>(l.ystart,t,h,l); 
    } 

  file1.close(); 
  file2.close(); 
  file3.close(); 
*/

