#include <unordered_map>
#include <string>
#include <utility>
#include <vector>

#include "mass.h"

const std::unordered_map<char, double> AA_MASSES{
	{'G', 57.02146372069},
	{'A', 71.03711378515},
	{'S', 87.03202840472},
	{'P', 97.05276384961},
	{'V', 99.06841391407},
	{'T', 101.04767846918},
	{'C', 103.00918495955},
	{'I', 113.08406397853},
	{'L', 113.08406397853},
	{'N', 114.04292744138},
	{'D', 115.02694302429},
	{'Q', 128.05857750584},
	{'K', 128.09496301519},
	{'E', 129.04259308875},
	{'M', 131.04048508847},
	{'H', 137.05891185847},
	{'F', 147.06841391407},
	{'R', 156.10111102405},
	{'Y', 163.06332853364},
	{'W', 186.07931295073}
};

// TODO: in just one loop
std::vector<double> calculateMass(const std::string& sequence, const std::vector<ModMassSite>& modSites) {
	// Position 0 is the N-term mass, Position sequence.size() + 1 is the C-term mass
	size_t seqLen = sequence.size();
	std::vector<double> siteModMasses(seqLen + 2, 0);
	for (const auto& modSite : modSites) {
		siteModMasses[modSite.first] += modSite.second;
	}

	std::vector<double> seqMasses(seqLen + 2);
	seqMasses[0] = siteModMasses[0];
	seqMasses[seqLen + 1] = siteModMasses[seqLen + 1];
	for (size_t ii = 0; ii < seqLen; ii++) {
		seqMasses[ii + 1] = AA_MASSES.at(sequence[ii]) + siteModMasses[ii + 1];
	}

	return seqMasses;
}
