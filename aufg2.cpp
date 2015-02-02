#include <iostream>
#include <fstream>
#include <omp.h>
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
	// Setup OpenMP
	#ifdef _OPENMP
		omp_set_num_threads( omp_get_num_procs() );
		cout << "OMP loaded with #" << omp_get_num_procs() << " threads" << endl;
	#endif

	const double tol = 1e-10;
	const double time = 1e5;
	ofstream file1("stepperDopr853.txt", ios::trunc);
	ofstream file2("stepperBS.txt", ios::trunc);

	Lorentz l1(0,0,0,1,0,0,1,3,-1);
	Output out1(10000);
	Odeint<StepperDopr853 <Lorentz> > dgl1(l1.ystart, 0., time, tol, tol, 0.01, 0, out1, l1);

	Lorentz l2(0,0,0,1,0,0,1,3,-1);
	Output out2(10000);
	Odeint<StepperBS <Lorentz> > dgl2(l2.ystart, 0., time, tol, tol, 0.01, 0, out2, l2);

	#pragma omp parallel sections
	{
			#pragma omp section
				dgl1.integrate();

			#pragma omp section
				dgl2.integrate();
	}

	cout << "Akzeptiert (Dopr853): " << out1.steps[0] << endl;
	cout << "Zurückgewiesen (Dopr853): " << out1.steps[1] << endl;

	cout << "Akzeptiert (BS): " << out2.steps[0] << endl;
	cout << "Zurückgewiesen (BS): " << out2.steps[1] << endl;

	for (int i=0; i < out1.count; i++)
	{
		file1 << out1.xsave[i] << '\t' << out1.ysave[0][i] << '\t' << out1.ysave[1][i] << '\t' << out1.ysave[2][i] << '\t' << out1.ysave[3][i] << '\t' << out1.ysave[4][i] << '\t' << out1.ysave[5][i] << endl;
		file2 << out2.xsave[i] << '\t' << out2.ysave[0][i] << '\t' << out2.ysave[1][i] << '\t' << out2.ysave[2][i] << '\t' << out2.ysave[3][i] << '\t' << out2.ysave[4][i] << '\t' << out2.ysave[5][i] << endl;
	}

	file1.close();
	file2.close();
}
