#ifndef _PEPFRAG_CONVERTERS_H
#define _PEPFRAG_CONVERTERS_H

#include <Python.h>
#include <string>
#include <unordered_map>
#include <vector>

#include "ion.h"
#include "mass.h"

std::vector<double> listToDoubleVector(PyObject* source);

std::vector<std::string> listToStringVector(PyObject* source);

std::unordered_map<IonType, std::vector<std::string>> dictToIonTypeMap(PyObject* source);

std::vector<ModMassSite> modSiteListToVector(PyObject* source, size_t seqLen);

PyObject* ionVectorToList(const std::vector<Ion>& ions);

PyObject* doubleVectorToList(const std::vector<double>& data);

#endif // _PEPFRAG_CONVERTERS_H
