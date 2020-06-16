#ifndef _PEPFRAG_ION_H
#define _PEPFRAG_ION_H

#include <Python.h>
#include <string>
#include <vector>

enum class IonType {
	precursor = 1,
	immonium = 2,
	b = 3,
	y = 4,
	a = 5,
	c = 6,
	z = 7,
	x = 8
};

// Forward declaration
struct Ion;

using Ions = std::vector<Ion>;
		 
struct Ion {
	double mass;
	std::string label;
	long position;
	
	Ion(double _mass, const std::string& _label, long _position)
		: mass(_mass), label(_label), position(_position) {}
		
	explicit operator PyObject*() const {
		PyObject* pMass = PyFloat_FromDouble(mass);
		PyObject* pLabel = PyUnicode_FromString(label.c_str());
		PyObject* pPosition = PyLong_FromLong(position);

		PyObject* tuple = PyTuple_Pack(3, pMass, pLabel, pPosition);

		Py_DECREF(pMass);
		Py_DECREF(pLabel);
		Py_DECREF(pPosition);

		return tuple;
	}
};

#endif // _PEPFRAG_ION_H
