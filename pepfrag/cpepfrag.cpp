#include <Python.h>
#include <string>
#include <unordered_map>
#include <vector>

#include "cpepfrag.h"
#include "converters.h"
#include "iongenerator.h"
#include "ion.h"
#include "mass.h"

std::vector<Ion> generateIons(
	IonType type,
	const std::vector<double>& masses,
	long charge,
	const std::vector<std::string>& neutralLosses,
	bool radical,
	const std::string& sequence)
{
	IonGeneratorPtr generator = IonGenerator::create(type);
	return generator->generate(masses, charge, neutralLosses, radical, sequence);
}

PyObject* python_generateIons(PyObject* module, PyObject* args) {
	PyObject *ionTypes = NULL, *precMass = NULL, *pySeqMasses = NULL, *bMassList = NULL,
	         *yMassList = NULL, *pyCharge = NULL, *pyRadical = NULL, *pySequence = NULL;
	
	try {
		if (!PyArg_UnpackTuple(args, "cpython_generateIons", 8, 8,
							   &ionTypes, &precMass, &pySeqMasses, &bMassList,
							   &yMassList, &pyCharge,
							   &pyRadical, &pySequence)) return NULL;
	} catch (const std::exception& ex) {
		PyErr_SetString(PyExc_RuntimeError, ex.what());
	}
	
	std::unordered_map<IonType, std::vector<std::string>> ionConfigs = dictToIonTypeMap(ionTypes);
	
	long charge = PyLong_AsLong(pyCharge);
	bool radical = PyObject_IsTrue(pyRadical);
	std::string sequence = PyUnicode_AsUTF8(pySequence);
	
	std::vector<double> seqMasses = listToDoubleVector(pySeqMasses);
	std::vector<double> bMasses = listToDoubleVector(bMassList);
	std::vector<double> yMasses = listToDoubleVector(yMassList);
	std::vector<double> precMasses = std::vector<double>{ PyFloat_AsDouble(precMass) };
	
	Ions ions;
	ions.reserve(1000);
	for (const auto& pair : ionConfigs) {
		const std::vector<double>* massList;
		switch (pair.first) {
			case IonType::b:
			case IonType::a:
			case IonType::c:
				massList = &bMasses;
				break;
			case IonType::y:
			case IonType::z:
				massList = &yMasses;
				break;
			case IonType::immonium:
				massList = &seqMasses;
				break;
			case IonType::precursor:
				massList = &precMasses;
				break;
		}
		
		mergeIonVectors(ions, generateIons(pair.first, *massList, charge, pair.second,
			                               radical, sequence));
	}
	
	return ionVectorToList(ions);
}

PyObject* python_calculateMass(PyObject* module, PyObject* args) {
	PyObject* sequence = NULL;
	PyObject* modSites = NULL;

	if (!PyArg_UnpackTuple(args, "cpython_calculateMass", 2, 2, &sequence, &modSites)) return NULL;

	std::string seq = PyUnicode_AsUTF8(sequence);

	return doubleVectorToList(calculateMass(
		seq,
		modSiteListToVector(modSites, seq.size())
	));
};

// Boilerplate code for C++ extension

static PyMethodDef cpepfrag_methods[] = {
	{"generate_ions", python_generateIons, METH_VARARGS, "Fragment ion generation."},
	{"calculate_mass", python_calculateMass, METH_VARARGS, "Peptide mass calculation"},
	{NULL, NULL, 0, NULL} /* SENTINEL */
};

PyDoc_STRVAR(cpepfragModuleDoc, "CPython functions for pepfrag");

static struct PyModuleDef cpepfragExt = {
	PyModuleDef_HEAD_INIT,
	"cpepfrag",
	cpepfragModuleDoc,
	-1,
	cpepfrag_methods
};

PyMODINIT_FUNC PyInit_cpepfrag(void) {
	return PyModule_Create(&cpepfragExt);
}
