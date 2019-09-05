#include <Python.h>
#include <string>
#include <unordered_map>
#include <vector>

#include "cpepfrag.h"
#include "converters.h"
#include "iongenerator.h"
#include "ion.h"
#include "mass.h"

Ions generateIons(
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
	PyObject *ionTypes, *pySeqMasses, *bMassList, *yMassList, *pySequence;
	double precMass;
	long charge;
	int radical;
	
	try {
		if (!PyArg_ParseTuple(args, "OdO&OOliO", &ionTypes, &precMass, &pySeqMasses, &bMassList,
				      &yMassList, &charge, &radical, &pySequence)) return NULL;
	} catch (const std::exception& ex) {
		PyErr_SetString(PyExc_RuntimeError, ex.what());
	}
	
	std::unordered_map<IonType, std::vector<std::string>> ionConfigs = dictToIonTypeMap(ionTypes);
	
	std::string sequence = PyUnicode_AsUTF8(pySequence);
	
	std::vector<double> seqMasses = listToDoubleVector(pySeqMasses);
	std::vector<double> bMasses = listToDoubleVector(bMassList);
	std::vector<double> yMasses = listToDoubleVector(yMassList);
	std::vector<double> precMasses = std::vector<double>{ precMass };
	
	Ions ions;
	ions.reserve(1000);
	const std::vector<double>* massList;
	for (const auto& pair : ionConfigs) {
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
			default:
				PyErr_SetString(PyExc_RuntimeError, "Invalid ion type specified");
				return NULL;
		}
		
		mergeIonVectors(ions, generateIons(pair.first, *massList, charge, pair.second,
			                               (bool) radical, sequence));
	}
	
	return vectorToList<Ion>(ions);
}

PyObject* python_calculateMass(PyObject* module, PyObject* args) {
	PyObject *sequence, *modSites;

	if (!PyArg_UnpackTuple(args, "calculateMass", 2, 2, &sequence, &modSites)) return NULL;

	std::string seq = PyUnicode_AsUTF8(sequence);

	return vectorToList(calculateMass(
		seq,
		modSiteListToVector(modSites, seq.size())
	), &PyFloat_FromDouble);
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
