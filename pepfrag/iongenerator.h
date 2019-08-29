#ifndef _PEPFRAG_IONGENERATOR_H
#define _PEPFRAG_IONGENERATOR_H

#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "ion.h"

const std::map<std::string, double> FIXED_MASSES = {
	{"H", 1.007276466879},
	{"tag", 304.20536},
	{"H2O", 18.01056468403},
	{"CO", 27.99491461957},
	{"NH3", 17.02654910112},
	{"cys_c", 57.021464},
	{"CO2", 43.989830},
	{"N", 14.003074}
};

const double PROTON_MASS = FIXED_MASSES.at("H");

struct NotImplementedException : public std::exception
{
	const char* what () const throw ()
    {
    	return "Method not implemented";
    }
};

// Forward declaration
class IonGenerator;

using IonGeneratorPtr = std::shared_ptr<IonGenerator>;

class IonGenerator {
	protected:
		const std::string ionLabel;
	
	public:
		IonGenerator(const std::string& label);

		virtual ~IonGenerator() = 0;
		
		static IonGeneratorPtr create(IonType type);
	
		virtual std::vector<Ion> generate(
			const std::vector<double>& masses,
			long charge,
			const std::vector<std::string>& neutralLosses,
			bool radical,
			const std::string& sequence) const;
		
	private:
		virtual std::vector<double> preProcessMasses(
			const std::vector<double>& masses) const;
	
		virtual Ion generateBaseIon(
			double mass,
			long position,
			const std::string& sequence) const;
			
		virtual std::vector<Ion> generateRadicalIons(
			double mass,
			long position) const;
			
		virtual std::vector<Ion> generateNeutralLosses(
			double mass,
			long position,
			const std::vector<std::string>& neutralLosses) const;
			
		virtual double fixMass(double mass) const;
};

class BIonGenerator : public IonGenerator
{
	public:
		BIonGenerator();

		~BIonGenerator() override {};
		
	private:
		std::vector<Ion> generateRadicalIons(
			double mass,
			long position) const override;
			
		double fixMass(double mass) const override;
};

class YIonGenerator : public IonGenerator
{
	public:
		YIonGenerator();

		~YIonGenerator() override {};
		
	private:
		std::vector<Ion> generateRadicalIons(
			double mass,
			long position) const override;
			
		double fixMass(double mass) const override;
};

class AIonGenerator : public IonGenerator
{
	public:
		AIonGenerator();

		~AIonGenerator() override {};
		
	private:
		std::vector<Ion> generateRadicalIons(
			double mass,
			long position) const override;
			
		double fixMass(double mass) const override;
};

class CIonGenerator : public IonGenerator
{
	public:
		CIonGenerator();

		~CIonGenerator() override {};
		
	private:
		std::vector<Ion> generateRadicalIons(
			double mass,
			long position) const override;
			
		double fixMass(double mass) const override;
};

class ZIonGenerator : public IonGenerator
{
	public:
		ZIonGenerator();

		~ZIonGenerator() override {};
		
	private:
		std::vector<Ion> generateRadicalIons(
			double mass,
			long position) const override;
			
		double fixMass(double mass) const override;
};

class ImmoniumIonGenerator : public IonGenerator
{
	public:
		ImmoniumIonGenerator();

		~ImmoniumIonGenerator() override {};
		
	private:
		std::vector<double> preProcessMasses(
			const std::vector<double>& masses) const override;
	
		Ion generateBaseIon(
			double mass,
			long position,
			const std::string& sequence) const override;
			
		double fixMass(double mass) const override;
};

class PrecursorIonGenerator : public IonGenerator
{
	public:
		PrecursorIonGenerator();

		~PrecursorIonGenerator() override {};
		
		std::vector<Ion> generate(
			const std::vector<double>& masses,
			long charge,
			const std::vector<std::string>& neutralLosses,
			bool radical,
			const std::string& sequence) const override;
		
	private:
		std::vector<double> preProcessMasses(
			const std::vector<double>& masses) const override;
	
		Ion generateBaseIon(
			double mass,
			long position,
			const std::string& sequence) const override;
			
		std::vector<Ion> generateRadicalIons(
			double mass,
			long position) const override;
			
		std::vector<Ion> generateNeutralLosses(
			double mass,
			long position,
			const std::vector<std::string>& neutralLosses) const override;
			
		double fixMass(double mass) const override;
};

std::vector<Ion> chargeIons(const std::vector<Ion>& ions, long chargeState);

inline Ion generateNeutralLossIon(
	const std::string& typeChar,
	const std::string& neutralLoss,
	double mass,
	long position)
{
	return {
		mass - FIXED_MASSES.at(neutralLoss),
		"[" + typeChar + std::to_string(position + 1) + "-" + neutralLoss + "][+]",
		position + 1
	};
}

void mergeIonVectors(std::vector<Ion>& target, const std::vector<Ion>& source);

#endif // _PEPFRAG_IONGENERATOR_H
