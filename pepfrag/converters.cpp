#include <Python.h>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "converters.h"
#include "ion.h"

/*
* Define type checks since these are macros in the CPython source
*/

bool checkFloat(PyObject* obj) {
	return PyFloat_Check(obj);
}

bool checkString(PyObject* obj) {
	return PyUnicode_Check(obj);
}

std::string unicodeToString(PyObject* obj) {
	return std::string(PyUnicode_AsUTF8(obj));
}

/* Python to C++ */

template<class T>
std::vector<T> listToVector(PyObject* source, bool(*check)(PyObject*), T(*convert)(PyObject*)) {
	std::vector<T> data;
	if (PyList_Check(source)) {
		long size = (long) PyList_Size(source);
		data.reserve(size);
		for (Py_ssize_t ii = 0; ii < size; ii++) {
			PyObject* value = PyList_GetItem(source, ii);
			if (check(value)) {
				data.push_back(convert(value));
			}
			else {
				throw std::logic_error("Contained PyObject pointer was not expected type");
			}
		}
	}
	else {
		throw std::logic_error("PyObject pointer was not a list");
	}
	
	return data;
}

std::vector<double> listToDoubleVector(PyObject* source) {
	return listToVector<double>(source, &checkFloat, &PyFloat_AsDouble);
}

std::vector<std::string> listToStringVector(PyObject* source) {
	return listToVector<std::string>(source, &checkString, &unicodeToString);
}

std::unordered_map<IonType, std::vector<std::string>> dictToIonTypeMap(PyObject* source) {
	std::unordered_map<IonType, std::vector<std::string>> types;
	if (PyDict_Check(source)) {
		PyObject *key, *value;
		Py_ssize_t pos = 0;
		while (PyDict_Next(source, &pos, &key, &value)) {
			types.emplace(
				static_cast<IonType>( PyLong_AsLong( key ) ),
				listToStringVector(value));
		}
	}
	else {
		throw std::logic_error("PyObject pointer was not a dict");
	}
	return types;
}

/* C++ to Python */

PyObject* ionVectorToList(const std::vector<Ion>& ions) {
	PyObject* listObj = PyList_New(ions.size());

	for (int ii = 0; ii < (int) ions.size(); ii++) {
		PyList_SetItem(listObj, ii, (PyObject*) ions[ii]);
	}
	
	return listObj;
}
