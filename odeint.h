#ifndef _ODEINT_H_
#define _ODEINT_H_

#include <climits>
#include "helfer.h"

// **********
class StepperBase
{
	public:
	double &x;
	double xold;
	VecDoub &y, &dydx;
	double atol, rtol;
	bool dense;
	double hdid;
	double hnext;
	double EPS;
	int n, neqn;
	VecDoub yout, yerr;

// 	Parent class for stepper classes StepperBS and StepperDopr853
	StepperBase(VecDoub &yy, VecDoub &dydxx, double &xx, const double atoll, const double rtoll, bool dens) :
		x(xx), y(yy), dydx(dydxx), atol(atoll), rtol(rtoll), dense(dens), n(y.size()), neqn(n), yout(n), yerr(n)
	{
	}
};

// **********
struct Output
{
	int kmax;
	int nvar;
	int nsave;
	bool dense;
	int count;
	double x1,x2,xout,dxout;

	inline void init(const int neqn, const double xlo, const double xhi);
	inline void resize(void);

	template <class Stepper>
	inline void save_dense(Stepper &s, const double xout, const double h);
	inline void save(const double x, VecDoub &y);

	template <class Stepper>
	inline void out(const int nstp, const double x, VecDoub &y, Stepper &s, const double h);

	public:
	int steps[2];

	Output(void);
	Output(const int nsavee);

	VecDoub xsave;
	MatDoub ysave;
};

// **********
template <class Stepper>
class Odeint
{
	friend class Output;

	static const int MAXSTP=1000000000;
	double EPS;
	int nvar;
	double x1,x2,hmin;
	bool dense;
	VecDoub y,dydx;
	VecDoub &ystart;
	Output &out;
	typename Stepper::Dtype &derivs;
	Stepper s;
	int nstp;
	double x,h;

	public:
	int nok;
	int nbad;

	// Constructor
	Odeint(VecDoub &ystartt, const double xx1, const double xx2, const double atol, const double rtol, \
	const double h1, const double hminn, Output& outt, typename Stepper::Dtype &derivss);

	// Integration routine
	void integrate(void);
};

// **********
template <class Stepper>
Odeint<Stepper>::Odeint(VecDoub &ystartt, const double xx1, const double xx2, const double atol, \
	const double rtol, const double h1, const double hminn, Output& outt, typename Stepper::Dtype &derivss) :

	nvar(ystartt.size()), y(nvar), dydx(nvar), ystart(ystartt), x(xx1), nok(0), nbad(0), \
	x1(xx1), x2(xx2), hmin(hminn), dense(outt.dense), out(outt), derivs(derivss), \
	s(y,dydx,x,atol,rtol,dense)
{
	EPS = std::numeric_limits<double>::epsilon();
	h = SIGN(h1,x2-x1);
	for(int i=0; i<nvar; i++)
		y[i]=ystart[i];
	out.init(s.neqn,x1,x2);
}

// **********
template <class Stepper>
void Odeint<Stepper>::integrate()
{
	derivs(x,y,dydx);

	// Store initial values
	if (dense)
		out.out(-1,x,y,s,h);
	else
		out.save(x,y);
	for (nstp=0; nstp<MAXSTP; nstp++)
	{

		// If stepsize can overshoot, decrease
		if ((x+h*1.0001-x2)*(x2-x1)>0.0)
			h = x2-x;

		// Take a step
		s.step(h,derivs);
		if (s.hdid == h)
			++nok;
		else
			++nbad;
		if (dense)
			out.out(nstp,x,y,s,s.hdid);
		else
			out.save(x,y);
		if ((x-x2)*(x2-x1) >= 0.0)
		{
			for (int i=0; i<nvar; i++)
				ystart[i]=y[i];
			if (out.kmax > 0 && abs(out.xsave[out.count-1]-x2) > 100.0*abs(x2)*EPS)
				out.save(x,y);
			out.steps[0] = nok;
			out.steps[1] = nbad;
			return;
		}
		if (abs(s.hnext) <= hmin)
			throw ("Step size too small in Odeint");
		h=s.hnext;
	}
	throw ("Too many steps in routine Odeint");
}

// **********
Output::Output() : kmax(-1), dense(false), count(0)
{
}

// **********
Output::Output(const int nsavee) : kmax(500), nsave(nsavee), count(0), xsave(kmax)
{
	dense = nsave > 0 ? true : false;
}

// **********
void Output::init(const int neqn, const double xlo, const double xhi)
{
	nvar=neqn;
	if(kmax == -1)
		return;
	ysave.resize(nvar,kmax);
	if (dense)
	{
		x1=xlo;
		x2=xhi;
		xout=x1;
		dxout=(x2-x1)/nsave;
	}
}

// **********
void Output::resize()
{
	int kold=kmax;
	kmax *= 2;
	VecDoub tempvec(xsave);
	xsave.resize(kmax);
	for (int k=0; k<kold; k++)
		xsave[k]=tempvec[k];
	MatDoub tempmat(ysave);
	ysave.resize(nvar,kmax);
	for (int i=0; i<nvar; i++)
		for (int k=0; k<kold; k++)
			ysave[i][k]=tempmat[i][k];
}

// **********
template <class Stepper>
void Output::save_dense(Stepper &s, const double xout, const double h)
{
	if (count == kmax)
		resize();
	for (int i=0; i<nvar; i++)
		ysave[i][count]=s.dense_out(i,xout,h);
	xsave[count++]=xout;
}

// **********
void Output::save(const double x, VecDoub &y)
{
	if (kmax <= 0)
		return;
	if (count == kmax)
		resize();
	for (int i=0; i<nvar; i++)
		ysave[i][count]=y[i];
	xsave[count++]=x;
}

// **********
template <class Stepper>
void Output::out(const int nstp, const double x, VecDoub &y, Stepper &s, const double h)
{
	if (!dense)
		throw ("Dense output not set in Output!");
	if (nstp == -1)
	{
		save(x,y);
		xout += dxout;
	}
	else
	{
		while ((x-xout)*(x2-x1) > 0.0)
		{
			save_dense(s,xout,h);
			xout += dxout;
		}
	}
}

#endif
