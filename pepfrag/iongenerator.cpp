#include <algorithm>
#include <iterator>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "iongenerator.h"
#include "ion.h"

// Defined for use with std::upper_bound, called by std::inplace_merge
bool operator<(const Ion& left, const Ion& right)
{
	return left.position < right.position;
}

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
		case IonType::precursor:
			return std::make_shared<PrecursorIonGenerator>(PrecursorIonGenerator());
		case IonType::immonium:
			return std::make_shared<ImmoniumIonGenerator>(ImmoniumIonGenerator());
	}
	return NULL;
}

Ions IonGenerator::generate(
	const std::vector<double>& masses,
	long charge,
	const std::vector<std::string>& neutralLosses,
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

std::pair<int, int> IonGenerator::preProcessMasses(const std::vector<double>& masses) const {
	return std::make_pair(0, masses.size() - 1);
}

Ion IonGenerator::generateBaseIon(double mass, long position, const std::string& /*sequence*/) const {
	return {mass, ionLabel + std::to_string(position + 1) + "[+]", position + 1};
}

void IonGenerator::generateRadicalIons(Ions& ions, double mass, long position) const {
}

void IonGenerator::generateNeutralLosses(
	Ions& ions,
	double mass,
	long position,
	const std::vector<std::string>& neutralLosses) const
{
	for (const std::string& neutralLoss : neutralLosses) {
		ions.push_back(generateNeutralLossIon(ionLabel, neutralLoss, mass, position));
	}
}

double IonGenerator::fixMass(double mass) const {
	return mass;
}

/* BIonGenerator */

BIonGenerator::BIonGenerator() : IonGenerator("b") {}

void BIonGenerator::generateRadicalIons(Ions& ions, double mass, long position) const {
	ions.emplace_back(
		mass,
		"[" + ionLabel + std::to_string(position + 1) + "-H[•+]",
		position + 1
	);
}

double BIonGenerator::fixMass(double mass) const {
	return mass + PROTON_MASS;
}

/* YIonGenerator */

YIonGenerator::YIonGenerator() : IonGenerator("y") {}

double YIonGenerator::fixMass(double mass) const {
	return mass + PROTON_MASS;
}

/* AIonGenerator */

AIonGenerator::AIonGenerator() : IonGenerator("a") {}

void AIonGenerator::generateRadicalIons(Ions& ions, double mass, long position) const {
	ions.emplace_back(
		mass - PROTON_MASS,
		"[" + ionLabel + std::to_string(position + 1) + "-H][•+]",
		position + 1
	);
	ions.emplace_back(
		mass + PROTON_MASS,
		"[" + ionLabel + std::to_string(position + 1) + "+H][•+]",
		position + 1
	);
}

double AIonGenerator::fixMass(double mass) const {
	return mass + PROTON_MASS - FIXED_MASSES.at("CO");
}

/* CIonGenerator */

CIonGenerator::CIonGenerator() : IonGenerator("c") {}

void CIonGenerator::generateRadicalIons(Ions& ions, double mass, long position) const {
	ions.emplace_back(
		mass + 2 * PROTON_MASS,
		"[" + ionLabel + std::to_string(position + 1) + "+2H][•+]",
		position + 1
	);
}

double CIonGenerator::fixMass(double mass) const {
	return mass + 3 * PROTON_MASS + FIXED_MASSES.at("N");
}

/* ZIonGenerator */

ZIonGenerator::ZIonGenerator() : IonGenerator("z") {}

void ZIonGenerator::generateRadicalIons(Ions& ions, double mass, long position) const {
	ions.emplace_back(
		mass - PROTON_MASS,
		"[" + ionLabel + std::to_string(position + 1) + "-H][•+]",
		position + 1
	);
}

double ZIonGenerator::fixMass(double mass) const {
	return mass - FIXED_MASSES.at("N") - 3 * PROTON_MASS;
}

/* ImmoniumIonGenerator */

ImmoniumIonGenerator::ImmoniumIonGenerator() : IonGenerator("imm") {}

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
	const std::vector<std::string>& neutralLosses,
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
		std::string chargeSymbol = (radical ? "•" : "") + (cs > 1 ? std::to_string(cs) : "") + "+";
		
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
		
		for (const std::string& neutralLoss : neutralLosses) {
			ions.emplace_back(
				(mass - FIXED_MASSES.at(neutralLoss)) / (double) cs + PROTON_MASS,
				"[" + ionLabel + "-" + neutralLoss + "][" + chargeSymbol + "]",
				seqLen);
		}
	}
	
	return ions;
}

std::pair<int, int> PrecursorIonGenerator::preProcessMasses(
	const std::vector<double>& /*masses*/) const
{
	// This method intentionally does nothing and should not be called
	throw NotImplementedException();
}
	
Ion PrecursorIonGenerator::generateBaseIon(
	double /*mass*/,
	long /*position*/,
	const std::string& /*sequence*/) const
{
	// This method intentionally does nothing and should not be called
	throw NotImplementedException();
}
	
void PrecursorIonGenerator::generateRadicalIons(
	Ions& /*ions*/,
	double /*mass*/,
	long /*position*/) const 
{
	// This method intentionally does nothing and should not be called
	throw NotImplementedException();
}

void PrecursorIonGenerator::generateNeutralLosses(
	Ions& /*ions*/,
	double /*mass*/,
	long /*position*/,
	const std::vector<std::string>& /*neutralLosses*/) const
{
	// This method intentionally does nothing and should not be called
	throw NotImplementedException();
}

double PrecursorIonGenerator::fixMass(double /*mass*/) const {
	// This method intentionally does nothing and should not be called
	throw NotImplementedException();
}

/* Utility functions */

void chargeIons(const Ions& sourceIons, Ions& target, long chargeState) {
	double hMass = PROTON_MASS * (chargeState - 1);
	long minPos = 2 * chargeState - 1;
	std::string chargeStr = std::to_string(chargeState) + "+";
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
