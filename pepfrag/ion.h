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
	z = 7
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
		return PyTuple_Pack(
			3,
			PyFloat_FromDouble(mass),
			PyUnicode_FromString(label.c_str()),
			PyLong_FromLong(position)
		);
	}
};

#endif // _PEPFRAG_ION_H
