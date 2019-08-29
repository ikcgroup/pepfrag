#ifndef _PEPFRAG_CONVERTERS_H
#define _PEPFRAG_CONVERTERS_H

#include <Python.h>
#include <string>
#include <unordered_map>
#include <vector>

#include "ion.h"

std::vector<double> listToDoubleVector(PyObject* source);

std::vector<std::string> listToStringVector(PyObject* source);

std::unordered_map<IonType, std::vector<std::string>> dictToIonTypeMap(PyObject* source);

PyObject* ionVectorToList(const std::vector<Ion>& ions);

#endif // _PEPFRAG_CONVERTERS_H