#ifndef _PEPFRAG_MASS_H
#define _PEPFRAG_MASS_H

#include <string>
#include <utility>
#include <vector>

using ModMassSite = std::pair<long, double>;

std::vector<double> calculateMass(const std::string& sequence, const std::vector<ModMassSite>& modSites);

#endif // _PEPFRAG_MASS_H
