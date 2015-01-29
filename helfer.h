#ifndef _HELFER_H_
#define _HELFER_H_

#include <ctime>
#include <cmath>
#include <limits>
#include <cstdlib>
using namespace std;

// Auxiliary function: square
// **********
template<class T>
inline const T& SQR(const T& a)
{
	return a*a;
}

// Auxiliary function: maximum
// **********
template<class T>
inline const T& MAX(const T& a, const T& b)
{
	return b>a ? b : a;
}

// Auxiliary function: minimum
template<class T>
inline const T& MIN(const T& a, const T& b)
{
	return b<a ? b : a;
}

// Auxiliary function: swapping
// **********
template<class T>
inline void SWAP(T& a, T& b)
{
	const T dummy = a;
	a = b;
	b = dummy;
}

// Auxiliary function: sign
// **********
template<class T>
inline T SIGN(const T &a, const T &b)
{
	return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

// Vector class of arbitrary length
// **********
template <class T>
class Vector
{
	int nn;
	T *v;

	public:
	Vector(void);
	explicit Vector(int n);
	Vector(int n, const T& a);
	Vector(int n, const T* a);
	Vector(const Vector& rhs);
	Vector& operator= (const Vector& rhs);
	typedef T value_type;
	inline T& operator[] (const int i);
	inline const T& operator[] (const int i) const;
	inline int size(void) const;
	void resize(int newn);
	void assign(int newn, const T& a);
	void assign(int newn, const T* a);
	~Vector(void);
};

// **********
template <class T>
Vector<T>::Vector() : nn(0), v(NULL) {}

template <class T>
Vector<T>::Vector(int n) : nn(n), v(n>0 ? new T[n] : NULL) {}

template <class T>
Vector<T>::Vector(int n, const T& a) : nn(n), v(n>0 ? new T[n] : NULL)
{
	for(int i=0; i<n; i++)
		v[i] = a;
}

template <class T>
Vector<T>::Vector(int n, const T* a) : nn(n), v(n>0 ? new T[n] : NULL)
{
	for(int i=0; i<n; i++)
		v[i] = *a++;
}

template <class T>
Vector<T>::Vector(const Vector<T>& rhs) : nn(rhs.nn), v(nn>0 ? new T[nn] : NULL)
{
	for(int i=0; i<nn; i++)
		v[i] = rhs[i];
}

template <class T>
Vector<T>& Vector<T>::operator= (const Vector<T>& rhs)
{
	if (this != &rhs)
	{
		if (nn != rhs.nn)
		{
			if (v != NULL)
				delete[] v;
			nn=rhs.nn;
			v= nn>0 ? new T[nn] : NULL;
		}
		for (int i=0; i<nn; i++)
			v[i]=rhs[i];
	}
	return *this;
}

template <class T>
T& Vector<T>::operator[] (const int i)
{
	return v[i];
}

template <class T>
const T& Vector<T>::operator[] (const int i) const
{
	return v[i];
}

template <class T>
int Vector<T>::size() const
{
	return nn;
}

template <class T>
void Vector<T>::resize(int newn)
{
	if (newn != nn)
	{
		if (v != NULL)
			delete[] v;
		nn = newn;
		v = nn > 0 ? new T[nn] : NULL;
	}
}

template <class T>
void Vector<T>::assign(int newn, const T& a)
{
	if (newn != nn)
	{
		if (v != NULL)
			delete[] v;
		nn = newn;
		v = nn > 0 ? new T[nn] : NULL;
	}
	for (int i=0; i<nn; i++)
		v[i] = a;
}

template <class T>
void Vector<T>::assign(int newn, const T* a)
{
	if (newn != nn)
	{
		if (v != NULL)
			delete[] v;
		nn = newn;
		v = nn > 0 ? new T[nn] : NULL;
	}
	for(int i=0; i<nn; i++)
		v[i] = *a++;
}

template <class T>
Vector<T>::~Vector()
{
	if (v != NULL)
		delete[] v;
}

// Matrix class of arbitrary length
// **********
template <class T>
class Matrix
{
	int nn;
	int mm;
	T **v;

	public:
	Matrix(const Matrix& rhs);
	Matrix(void);
	Matrix(int n, int m);
	Matrix(int n, int m, const T& a);
	Matrix(int n, int m, const T* a);
	Matrix& operator= (const Matrix& rhs);
	typedef T value_type;
	inline T* operator[] (const int i);
	inline const T* operator[] (const int i) const;
	inline int nrows(void) const;
	inline int ncols(void) const;
	void resize(int newn, int newm);
	void assign(int newn, int newm, const T& a);
	~Matrix(void);
};

// **********
template <class T>
Matrix<T>::Matrix() : nn(0), mm(0), v(NULL) {}

template <class T>
Matrix<T>::Matrix(int n, int m) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int nel = m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (int i=1; i<n; i++)
		v[i] = v[i-1] + m;
}

template <class T>
Matrix<T>::Matrix(int n, int m, const T& a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int nel = m*n;
	if (v) v[0] = nel>0 ?
		new T[nel] : NULL;
	for (int i=1; i< n; i++)
		v[i] = v[i-1] + m;
	for (int i=0; i< n; i++)
		for (int j=0; j<m; j++)
			v[i][j] = a;
}

template <class T>
Matrix<T>::Matrix(int n, int m, const T* a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int nel = m*n;
	if (v)
		v[0] = nel>0 ? new T[nel] : NULL;
	for (int i=1; i< n; i++)
		v[i] = v[i-1] + m;
	for (int i=0; i< n; i++)
		for (int j=0; j<m; j++)
			v[i][j] = *a++;
}

template <class T>
Matrix<T>::Matrix(const Matrix& rhs) : nn(rhs.nn), mm(rhs.mm), v(nn>0 ? new T*[nn] : NULL)
{
	int nel=mm*nn;
	if (v)
		v[0] = nel>0 ? new T[nel] : NULL;
	for (int i=1; i< nn; i++)
		v[i] = v[i-1] + mm;
	for (int i=0; i< nn; i++)
		for (int j=0; j<mm; j++)
			v[i][j] = rhs[i][j];
}

template <class T>
Matrix<T> & Matrix<T>::operator= (const Matrix<T>& rhs)
{
	if (this != &rhs)
	{
		if (nn != rhs.nn || mm != rhs.mm)
		{
			if (v != NULL)
			{
				delete[] (v[0]);
				delete[] (v);
			}
			nn=rhs.nn;
			mm=rhs.mm;
			v = nn>0 ? new T*[nn] : NULL;
			int nel = mm*nn;
			if (v)
				v[0] = nel>0 ? new T[nel] : NULL;
			for (int i=1; i< nn; i++)
				v[i] = v[i-1] + mm;
		}
		for (int i=0; i< nn; i++)
			for (int j=0; j<mm; j++)
				v[i][j] = rhs[i][j];
	}
	return *this;
}

template <class T>
T* Matrix<T>::operator[] (const int i)
{
	return v[i];
}

template <class T>
const T* Matrix<T>::operator[] (const int i) const
{
	return v[i];
}

template <class T>
int Matrix<T>::nrows() const
{
	return nn;
}

template <class T>
int Matrix<T>::ncols() const
{
	return mm;
}

template <class T>
void Matrix<T>::resize(int newn, int newm)
{
	if (newn != nn || newm != mm)
	{
		if (v != NULL)
		{
			delete[] v[0];
			delete[] v;
		}
		nn = newn;
		mm = newm;
		v = nn>0 ? new T*[nn] : NULL;
		int nel = mm*nn;
		if (v)
			v[0] = nel>0 ? new T[nel] : NULL;
		for (int i=1; i<nn; i++)
			v[i] = v[i-1] + mm;
	}
}

template <class T>
void Matrix<T>::assign(int newn, int newm, const T& a)
{
	if (newn != nn || newm != mm)
	{
		if (v != NULL)
		{
			delete[] v[0];
			delete[] v;
		}
		nn = newn;
		mm = newm;
		v = nn>0 ? new T*[nn] : NULL;
		int nel = mm*nn;
		if (v)
			v[0] = nel>0 ? new T[nel] : NULL;
		for (int i=1; i< nn; i++)
			v[i] = v[i-1] + mm;
	}
	for (int i=0; i<nn; i++)
		for (int j=0; j<mm; j++)
			v[i][j] = a;
}

template <class T>
Matrix<T>::~Matrix()
{
	if (v != NULL)
	{
		delete[] v[0];
		delete[] v;
	}
}

// Definitions of Vector/Matrix with integer/double entries
// **********
typedef Vector<double> VecDoub;
typedef Vector<int> VecInt;
typedef Matrix<double> MatDoub;
typedef Matrix<int> MatInt;

#endif
