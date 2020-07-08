#include <functional>
#include <map>
#include <unordered_map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "mass.h"

const std::unordered_map< char, std::pair<double, double> > AA_MASSES{
	{'G', std::pair<double, double> {57.02146372069, 57.051402191402}},
	{'A', std::pair<double, double> {71.03711378515, 71.078019596249}},
	{'S', std::pair<double, double> {87.03202840472, 87.077424520567}},
	{'P', std::pair<double, double> {97.05276384961, 97.115372897831}},
	{'V', std::pair<double, double> {99.06841391407, 99.131254405943}},
	{'T', std::pair<double, double> {101.04767846918, 101.104041925414}},
	{'C', std::pair<double, double> {103.00918495955, 103.142807002376}},
	{'I', std::pair<double, double> {113.08406397853, 113.157871810790}},
	{'L', std::pair<double, double> {113.08406397853, 113.157871810790}},
	{'N', std::pair<double, double> {114.04292744138, 114.102804382804}},
	{'D', std::pair<double, double> {115.02694302429, 115.087565341620}},
	{'Q', std::pair<double, double> {128.05857750584, 128.129421787651}},
	{'K', std::pair<double, double> {128.09496301519, 128.172515776292}},
	{'E', std::pair<double, double> {129.04259308875, 129.114182746467}},
	{'M', std::pair<double, double> {131.04048508847, 131.19604181207}},
	{'H', std::pair<double, double> {137.05891185847, 137.139515217458}},
	{'F', std::pair<double, double> {147.06841391407, 147.174197992883}},
	{'R', std::pair<double, double> {156.10111102405, 156.185922199184}},
	{'Y', std::pair<double, double> {163.06332853364, 163.173602917201}},
	{'W', std::pair<double, double> {186.07931295073, 186.210313751855}}
};

const std::vector< std::function<double(char)> > GET_MASS_FUNCTIONS{
    [](char r) { return AA_MASSES.at(r).first; },
    [](char r) { return AA_MASSES.at(r).second; }
};

ModMassSite::ModMassSite(long _site, double _mass)
	: site(_site), mass(_mass) {}

std::vector<double> calculateMass(
    const std::string& sequence,
    const std::map<long, double>& modSiteMasses,
    long massType
) {
    auto massGetter = GET_MASS_FUNCTIONS[massType];

	size_t seqLen = sequence.size();

	std::vector<double> seqMasses(seqLen + 2);
	// Position 0 is the N-term mass, Position sequence.size() + 1 is the C-term mass
	auto it = modSiteMasses.find(0);
	seqMasses[0] = it != modSiteMasses.end() ? it->second : 0;
    it = modSiteMasses.find(seqLen + 1);
    seqMasses[seqLen + 1] = it != modSiteMasses.end() ? it->second : 0;
	for (size_t ii = 0; ii < seqLen; ii++) {
	    double residueMass;
	    try {
	        residueMass = massGetter(sequence[ii]);
	    }
	    catch (const std::out_of_range& ex) {
	        throw std::out_of_range("Invalid residue detected: " + std::string(1, sequence[ii]));
	    }
		seqMasses[ii + 1] = residueMass;

		it = modSiteMasses.find(ii + 1);
		if (it != modSiteMasses.end()) {
		    seqMasses[ii + 1] += it->second;
		}
	}

	return seqMasses;
}
