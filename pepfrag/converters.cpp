#include <Python.h>
#include <stdexcept>
#include <string>
#include <utility>
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
	if (!PySequence_Check(source)) {
		throw std::logic_error("PyObject pointer was not a sequence");
	}
	
	long size = (long) PySequence_Size(source);
	std::vector<T> data;
	data.reserve(size);
	for (Py_ssize_t ii = 0; ii < size; ii++) {
		PyObject* value = PySequence_GetItem(source, ii);
		if (check(value)) {
			data.push_back(convert(value));
		}
		else {
			Py_DECREF(value);
			throw std::logic_error("Contained PyObject pointer was not expected type");
		}
		Py_DECREF(value);
	}
	return data;
}

std::vector<double> listToDoubleVector(PyObject* source) {
	return listToVector<double>(source, &checkFloat, &PyFloat_AsDouble);
}

std::vector<std::string> listToStringVector(PyObject* source) {
	return listToVector<std::string>(source, &checkString, &unicodeToString);
}

std::vector<std::pair<IonType, std::vector<std::string>>> dictToIonTypeMap(PyObject* source) {
	if (!PyDict_Check(source)) {
		throw std::logic_error("PyObject pointer was not a dict");
	}
	
	std::vector<std::pair<IonType, std::vector<std::string>>> types;
	types.reserve((long) PyDict_Size(source));
	PyObject *key, *value;
	Py_ssize_t pos = 0;
	while (PyDict_Next(source, &pos, &key, &value)) {
		types.emplace_back(
			static_cast<IonType>( PyLong_AsLong( key ) ),
                        listToStringVector(value)
		);
	}
	
	return types;
}

std::vector<ModMassSite> modSiteListToVector(PyObject* source, size_t seqLen) {
	if (!PySequence_Check(source)) {
		throw std::logic_error("PyObject pointer was not a sequence");
	}

	Py_ssize_t size = PySequence_Size(source);

	std::vector<ModMassSite> modSites;
	modSites.reserve(size);

	for (Py_ssize_t ii = 0; ii < size; ii++) {
		PyObject* tuple = PySequence_GetItem(source, ii);

		if (!PyTuple_Check(tuple)) {
			throw std::logic_error("PyObject pointer was not a tuple");
		}

		PyObject* site = PyTuple_GET_ITEM(tuple, 1);

		long siteIdx;
		if (PyLong_Check(site)) {
			siteIdx = PyLong_AsLong(site);
		}
		else {
			std::string siteStr = PyUnicode_AsUTF8(site);
			siteIdx = (siteStr == "N-term" || siteStr == "nterm") ? 0 : seqLen + 1;
		}

		modSites.emplace_back(
			siteIdx,
			PyFloat_AsDouble(PyTuple_GET_ITEM(tuple, 0))
		);

		Py_DECREF(tuple);
	}

	return modSites;
}
