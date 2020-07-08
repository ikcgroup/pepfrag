#ifndef _PEPFRAG_CONVERTERS_H
#define _PEPFRAG_CONVERTERS_H

#include <Python.h>
#include <string>
#include <utility>
#include <vector>
#include <map>

#include "ion.h"
#include "mass.h"

std::vector<double> listToDoubleVector(PyObject* source);

std::vector<std::string> listToStringVector(PyObject* source);

using IonTypeMap = std::vector<std::pair<IonType, std::vector<std::pair<std::string, double>>>>;

IonTypeMap dictToIonTypeMap(PyObject* source);

std::map<long, double> modSiteListToMap(PyObject* source, size_t seqLen);

/* C++ to Python */

template<class T>
PyObject* vectorToList(const std::vector<T>& data, PyObject*(*convert)(T)) {
        long size = (long) data.size();
        PyObject* listObj = PyList_New(size);
        for (long ii = 0; ii < size; ii++) {
                PyList_SET_ITEM(listObj, ii, convert(data[ii]));
        }
        return listObj;
}

template<class T>
PyObject* vectorToList(const std::vector<T>& data, PyObject*(*convert)(const T&)) {
        long size = (long) data.size();
        PyObject* listObj = PyList_New(size);
        for (long ii = 0; ii < size; ii++) {
                PyList_SET_ITEM(listObj, ii, convert(data[ii]));
        }
        return listObj;
}

template<class T>
PyObject* vectorToList(const std::vector<T>& data) {
        long size = (long) data.size();
        PyObject* listObj = PyList_New(size);
        for (long ii = 0; ii < size; ii++) {
                PyList_SET_ITEM(listObj, ii, (PyObject*) data[ii]);
        }
        return listObj;
}

#endif // _PEPFRAG_CONVERTERS_H
