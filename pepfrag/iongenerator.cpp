#include <algorithm>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

#include "iongenerator.h"
#include "ion.h"

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

std::vector<Ion> IonGenerator::generate(
	const std::vector<double>& _masses,
	long charge,
	const std::vector<std::string>& neutralLosses,
	bool radical,
	const std::string& sequence) const
{
	std::vector<double> masses = preProcessMasses(_masses);
	
	std::vector<Ion> ions;
	ions.reserve(masses.size() * 10);
	
	for (int ii = 0; ii < masses.size(); ii++) {
		double ion_mass = fixMass(masses[ii]);
		
		ions.push_back(generateBaseIon(ion_mass, ii, sequence));
		
		if (radical) {
			mergeIonVectors(ions, generateRadicalIons(ion_mass, ii));
		}
		
		if (!neutralLosses.empty()) {
			mergeIonVectors(ions, generateNeutralLosses(ion_mass, ii, neutralLosses));
		}
	}
	
	// Perform charging based on the singly charged ions in "ions" and update
	// "chargedIons" so that doubly charged ions are not re-submitted to chargeIons
	// since this results in ions like [23+]
	std::vector<Ion> chargedIons;
	chargedIons.reserve(ions.size() * charge);
	chargedIons = ions;
	for (long cs = 1; cs < charge; cs++) {
		mergeIonVectors(chargedIons, chargeIons(ions, cs + 1));
	}

	return chargedIons;
}

std::vector<double> IonGenerator::preProcessMasses(const std::vector<double>& masses) const {
	std::vector<double> _masses = masses;
	_masses.pop_back();
	return _masses;
}

Ion IonGenerator::generateBaseIon(double mass, long position, const std::string& /*sequence*/) const {
	return {mass, ionLabel + std::to_string(position + 1) + "[+]", position + 1};
}

std::vector<Ion> IonGenerator::generateRadicalIons(double mass, long position) const {
	return std::vector<Ion>();
}

std::vector<Ion> IonGenerator::generateNeutralLosses(
	double mass,
	long position,
	const std::vector<std::string>& neutralLosses) const
{
	std::vector<Ion> ions;
	ions.reserve(neutralLosses.size());
	for (const std::string& neutralLoss : neutralLosses) {
		ions.push_back(generateNeutralLossIon(ionLabel, neutralLoss, mass, position));
	}

	return ions;
}

double IonGenerator::fixMass(double mass) const {
	return mass;
}

/* BIonGenerator */

BIonGenerator::BIonGenerator() : IonGenerator("b") {}

std::vector<Ion> BIonGenerator::generateRadicalIons(double mass, long position) const {
	return std::vector<Ion>{
		Ion(mass,
			"[" + ionLabel + std::to_string(position + 1) + "-H[•+]",
			position + 1)
	};
}

double BIonGenerator::fixMass(double mass) const {
	return mass + PROTON_MASS;
}

/* YIonGenerator */

YIonGenerator::YIonGenerator() : IonGenerator("y") {}

std::vector<Ion> YIonGenerator::generateRadicalIons(double mass, long position) const {
	return std::vector<Ion>();
}

double YIonGenerator::fixMass(double mass) const {
	return mass + PROTON_MASS;
}

/* AIonGenerator */

AIonGenerator::AIonGenerator() : IonGenerator("a") {}

std::vector<Ion> AIonGenerator::generateRadicalIons(double mass, long position) const {
	return std::vector<Ion>{
		Ion(mass - PROTON_MASS, "[" + ionLabel + std::to_string(position + 1) + "-H][•+]", position + 1),
		Ion(mass + PROTON_MASS, "[" + ionLabel + std::to_string(position + 1) + "+H][•+]", position + 1)
	};
}

double AIonGenerator::fixMass(double mass) const {
	return mass + PROTON_MASS - FIXED_MASSES.at("CO");
}

/* CIonGenerator */

CIonGenerator::CIonGenerator() : IonGenerator("c") {}

std::vector<Ion> CIonGenerator::generateRadicalIons(double mass, long position) const {
	return std::vector<Ion>{
		Ion(mass + 2 * PROTON_MASS, "[" + ionLabel + std::to_string(position + 1) + "+2H][•+]", position + 1)
	};
}

double CIonGenerator::fixMass(double mass) const {
	return mass + 3 * PROTON_MASS + FIXED_MASSES.at("N");
}

/* ZIonGenerator */

ZIonGenerator::ZIonGenerator() : IonGenerator("z") {}

std::vector<Ion> ZIonGenerator::generateRadicalIons(double mass, long position) const {
	return std::vector<Ion>{
		Ion(mass - PROTON_MASS, "[" + ionLabel + std::to_string(position + 1) + "-H][•+]", position + 1)
	};
}

double ZIonGenerator::fixMass(double mass) const {
	return mass - FIXED_MASSES.at("N") - 3 * PROTON_MASS;
}

/* ImmoniumIonGenerator */

ImmoniumIonGenerator::ImmoniumIonGenerator() : IonGenerator("imm") {}

std::vector<double> ImmoniumIonGenerator::preProcessMasses(const std::vector<double>& masses) const {
	return masses;
}

Ion ImmoniumIonGenerator::generateBaseIon(double mass, long position, const std::string& sequence) const {
	return {
		mass,
		ionLabel + "(" + sequence.at(position) + ")",
		0
	};
}
	
double ImmoniumIonGenerator::fixMass(double mass) const {
	return mass - FIXED_MASSES.at("CO") + PROTON_MASS;
}

/* PrecursorIonGenerator */

PrecursorIonGenerator::PrecursorIonGenerator() : IonGenerator("M") {}

std::vector<Ion> PrecursorIonGenerator::generate(
	const std::vector<double>& masses,
	long charge,
	const std::vector<std::string>& neutralLosses,
	bool radical,
	const std::string& sequence) const
{
	std::vector<Ion> ions;
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

std::vector<double> PrecursorIonGenerator::preProcessMasses(
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
	
std::vector<Ion> PrecursorIonGenerator::generateRadicalIons(
	double /*mass*/,
	long /*position*/) const 
{
	// This method intentionally does nothing and should not be called
	throw NotImplementedException();
}
			
std::vector<Ion> PrecursorIonGenerator::generateNeutralLosses(
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

std::vector<Ion> chargeIons(const std::vector<Ion>& ions, long chargeState) {
	double hMass = PROTON_MASS * (chargeState - 1);
	std::vector<Ion> chargedIons;
	chargedIons.reserve(ions.size());
	long minPos = 2 * chargeState - 1;
	for (const Ion& ion : ions) {
		if (ion.position >= minPos) {
			std::string label = ion.label;
			chargedIons.emplace_back(
				(ion.mass + hMass) / (double) chargeState,
				label.replace(ion.label.find('+'), 1, std::to_string(chargeState) + "+"),
				ion.position);
		}
	}
	return chargedIons;
}

void mergeIonVectors(std::vector<Ion>& target, const std::vector<Ion>& source) {
	size_t n = target.size();
	target.insert(target.end(), std::make_move_iterator(source.begin()), std::make_move_iterator(source.end()));
	std::inplace_merge(target.begin(), target.begin() + n, target.end());
}