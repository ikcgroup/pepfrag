#include <algorithm>
#include <iterator>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "iongenerator.h"
#include "ion.h"

// Defined for use with std::upper_bound, called by std::inplace_merge
bool operator<(const Ion& left, const Ion& right)
{
	return left.position < right.position;
}

/* StringCache */

class StringCache {
    static std::unordered_map<long, std::string> cache;

    public:

        static std::string get(long n) {
            auto it = cache.find(n);
            if (it != cache.end()) {
                return it->second;
            }
            std::string entry = std::to_string(n);
            cache[n] = entry;
            return entry;
        }
};

std::unordered_map<long, std::string> StringCache::cache = std::unordered_map<long, std::string>();

const std::string RADICAL = "•";

/* IonGenerator */

IonGenerator::IonGenerator(const std::string& label) : ionLabel(label) {}

IonGeneratorPtr IonGenerator::create(IonType type) {
	switch (type) {
		case IonType::b:
			return std::make_shared<BIonGenerator>(BIonGenerator());
		case IonType::y:
			return std::make_shared<YIonGenerator>(YIonGenerator());
		case IonType::a:
			return std::make_shared<AIonGenerator>(AIonGenerator());
		case IonType::c:
			return std::make_shared<CIonGenerator>(CIonGenerator());
		case IonType::z:
			return std::make_shared<ZIonGenerator>(ZIonGenerator());
		case IonType::x:
		    return std::make_shared<XIonGenerator>(XIonGenerator());
		case IonType::precursor:
			return std::make_shared<PrecursorIonGenerator>(PrecursorIonGenerator());
		case IonType::immonium:
			return std::make_shared<ImmoniumIonGenerator>(ImmoniumIonGenerator());
	}
	return NULL;
}

/* SimpleIonGenerator */

SimpleIonGenerator::SimpleIonGenerator(const std::string& label) : IonGenerator(label) {};

Ions SimpleIonGenerator::generate(
	const std::vector<double>& masses,
	long charge,
	const std::vector<NeutralLossPair>& neutralLosses,
	bool radical,
	const std::string& sequence) const
{
	std::pair<int, int> massIndices = preProcessMasses(masses);
	
	Ions ions;
	ions.reserve(masses.size() * 10);
	
	for (int ii = massIndices.first; ii < massIndices.second; ii++) {
		double ion_mass = fixMass(masses[ii]);
		
		ions.push_back(generateBaseIon(ion_mass, ii, sequence));
		
		if (radical) {
			generateRadicalIons(ions, ion_mass, ii);
		}
		
		if (!neutralLosses.empty()) {
			generateNeutralLosses(ions, ion_mass, ii, neutralLosses);
		}
	}
	
	// Perform charging based on the singly charged ions in "ions" and update
	// "chargedIons" so that doubly charged ions are not re-submitted to chargeIons
	// since this results in ions like [23+]
	Ions chargedIons;
	chargedIons.reserve(ions.size() * charge);
	chargedIons = ions;
	for (long cs = 1; cs < charge; cs++) {
		chargeIons(ions, chargedIons, cs + 1);
	}

	return chargedIons;
}

std::pair<int, int> SimpleIonGenerator::preProcessMasses(const std::vector<double>& masses) const {
	return std::make_pair(0, masses.size() - 1);
}

Ion SimpleIonGenerator::generateBaseIon(double mass, long position, const std::string& /*sequence*/) const {
	return {mass, ionLabel + StringCache::get(position + 1) + "[+]", position + 1};
}

void SimpleIonGenerator::generateRadicalIons(Ions& ions, double mass, long position) const {
}

void SimpleIonGenerator::generateNeutralLosses(
	Ions& ions,
	double mass,
	long position,
	const std::vector<NeutralLossPair>& neutralLosses) const
{
	for (NeutralLossPair neutralLoss : neutralLosses) {
		ions.push_back(generateNeutralLossIon(ionLabel, neutralLoss, mass, position));
	}
}

double SimpleIonGenerator::fixMass(double mass) const {
	return mass;
}

/* BIonGenerator */

BIonGenerator::BIonGenerator() : SimpleIonGenerator("b") {}

void BIonGenerator::generateRadicalIons(Ions& ions, double mass, long position) const {
	ions.emplace_back(
		mass,
		"[" + ionLabel + StringCache::get(position + 1) + "-H[" + RADICAL + "+]",
		position + 1
	);
}

double BIonGenerator::fixMass(double mass) const {
	return mass + PROTON_MASS;
}

/* YIonGenerator */

YIonGenerator::YIonGenerator() : SimpleIonGenerator("y") {}

double YIonGenerator::fixMass(double mass) const {
	return mass + PROTON_MASS;
}

/* AIonGenerator */

AIonGenerator::AIonGenerator() : SimpleIonGenerator("a") {}

void AIonGenerator::generateRadicalIons(Ions& ions, double mass, long position) const {
	ions.emplace_back(
		mass - PROTON_MASS,
		"[" + ionLabel + StringCache::get(position + 1) + "-H][" + RADICAL + "+]",
		position + 1
	);
	ions.emplace_back(
		mass + PROTON_MASS,
		"[" + ionLabel + StringCache::get(position + 1) + "+H][" + RADICAL + "+]",
		position + 1
	);
}

double AIonGenerator::fixMass(double mass) const {
	return mass + PROTON_MASS - FIXED_MASSES.at("CO");
}

/* CIonGenerator */

CIonGenerator::CIonGenerator() : SimpleIonGenerator("c") {}

void CIonGenerator::generateRadicalIons(Ions& ions, double mass, long position) const {
	ions.emplace_back(
		mass + 2 * PROTON_MASS,
		"[" + ionLabel + StringCache::get(position + 1) + "+2H][" + RADICAL + "+]",
		position + 1
	);
}

double CIonGenerator::fixMass(double mass) const {
	return mass + 4 * PROTON_MASS + FIXED_MASSES.at("N");
}

/* ZIonGenerator */

ZIonGenerator::ZIonGenerator() : SimpleIonGenerator("z") {}

void ZIonGenerator::generateRadicalIons(Ions& ions, double mass, long position) const {
	ions.emplace_back(
		mass - PROTON_MASS,
		"[" + ionLabel + StringCache::get(position + 1) + "-H][" + RADICAL + "+]",
		position + 1
	);
}

double ZIonGenerator::fixMass(double mass) const {
	return mass - FIXED_MASSES.at("N") - 1 * PROTON_MASS;
}

/* XIonGenerator */

XIonGenerator::XIonGenerator() : SimpleIonGenerator("x") {}

void XIonGenerator::generateRadicalIons(Ions& ions, double mass, long position) const {
    ions.emplace_back(
		mass,
		"[" + ionLabel + StringCache::get(position + 1) + "-H[" + RADICAL + "+]",
		position + 1
	);
}

double XIonGenerator::fixMass(double mass) const {
    return mass + FIXED_MASSES.at("CO") - PROTON_MASS;
}

/* ImmoniumIonGenerator */

ImmoniumIonGenerator::ImmoniumIonGenerator() : SimpleIonGenerator("imm") {}

std::pair<int, int> ImmoniumIonGenerator::preProcessMasses(const std::vector<double>& masses) const {
	return std::make_pair(0, masses.size());
}

Ion ImmoniumIonGenerator::generateBaseIon(double mass, long position, const std::string& sequence) const {
	return {
		mass,
		ionLabel + "(" + sequence[position] + ")",
		0
	};
}
	
double ImmoniumIonGenerator::fixMass(double mass) const {
	return mass - FIXED_MASSES.at("CO") + PROTON_MASS;
}

/* PrecursorIonGenerator */

PrecursorIonGenerator::PrecursorIonGenerator() : IonGenerator("M") {}

Ions PrecursorIonGenerator::generate(
	const std::vector<double>& masses,
	long charge,
	const std::vector<NeutralLossPair>& neutralLosses,
	bool radical,
	const std::string& sequence) const
{
	Ions ions;
	ions.reserve(20);
	
	// Only use one mass - if multiple masses are passed to the PrecursorIonGenerator,
	// an exception needs to be thrown
	double mass = masses[0];
	
	long seqLen = (long) sequence.size();
	
	for (long cs = 1; cs < charge + 1; cs++) {
		std::string chargeSymbol = (radical ? RADICAL : "") + (cs > 1 ? StringCache::get(cs) : "") + "+";
		
		ions.emplace_back(
			(mass / (double) cs) + PROTON_MASS,
			"[" + ionLabel + "+H][" + chargeSymbol + "]",
			seqLen
		);

		if (radical) {
			ions.emplace_back(
				mass / (double) cs,
				ionLabel + "[" + chargeSymbol + "]",
				seqLen
			);
		}
		
		for (NeutralLossPair neutralLoss : neutralLosses) {
			ions.emplace_back(
				(mass - neutralLoss.second) / (double) cs + PROTON_MASS,
				"[" + ionLabel + "-" + neutralLoss.first + "][" + chargeSymbol + "]",
				seqLen);
		}
	}
	
	return ions;
}

/* Utility functions */

void chargeIons(const Ions& sourceIons, Ions& target, long chargeState) {
	double hMass = PROTON_MASS * (chargeState - 1);
	long minPos = 2 * chargeState - 1;
	std::string chargeStr = StringCache::get(chargeState) + "+";
	for (const Ion& ion : sourceIons) {
		if (ion.position >= minPos) {
			std::string label = ion.label;
			target.emplace_back(
				(ion.mass + hMass) / (double) chargeState,
				label.replace(ion.label.find('+'), 1, chargeStr),
				ion.position);
		}
	}
}

void mergeIonVectors(Ions& target, const Ions& source) {
	size_t n = target.size();
	target.insert(target.end(), std::make_move_iterator(source.begin()), std::make_move_iterator(source.end()));
	std::inplace_merge(target.begin(), target.begin() + n, target.end());
}
