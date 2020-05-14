#ifndef _PEPFRAG_MASS_H
#define _PEPFRAG_MASS_H

#include <string>
#include <vector>

struct ModMassSite {
	long site;
	double mass;
	
	ModMassSite(long _site, double _mass);
};

std::vector<double> calculateMass(
    const std::string& sequence,
    const std::vector<ModMassSite>& modSites,
    long massType
);

#endif // _PEPFRAG_MASS_H
