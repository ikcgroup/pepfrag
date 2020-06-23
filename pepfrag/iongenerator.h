#ifndef _PEPFRAG_IONGENERATOR_H
#define _PEPFRAG_IONGENERATOR_H

#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "ion.h"

using NeutralLossPair = std::pair<std::string, double>;

const std::unordered_map<std::string, double> FIXED_MASSES = {
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

// Forward declaration
class IonGenerator;

using IonGeneratorPtr = std::shared_ptr<IonGenerator>;

class IonGenerator {
    protected:
        const std::string ionLabel;

    public:
        explicit IonGenerator(const std::string& label);

        virtual ~IonGenerator() {};

        static IonGeneratorPtr create(IonType type);

        virtual Ions generate(
            const std::vector<double>& masses,
            long charge,
            const std::vector<NeutralLossPair>& neutralLosses,
            bool radical,
            const std::string& sequence) const = 0;
};

/*
 * This class implements a version of the generate method with behaviour
 * customized by the defined private methods for each subclass.
 */
class SimpleIonGenerator : public IonGenerator {
    public:
        explicit SimpleIonGenerator(const std::string& label);

        virtual Ions generate(
            const std::vector<double>& masses,
            long charge,
            const std::vector<NeutralLossPair>& neutralLosses,
            bool radical,
            const std::string& sequence) const override;

	private:
		virtual std::pair<int, int> preProcessMasses(
			const std::vector<double>& masses) const;
	
		virtual Ion generateBaseIon(
			double mass,
			long position,
			const std::string& sequence) const;
			
		virtual void generateRadicalIons(
			Ions& ions,
			double mass,
			long position) const;
			
		virtual void generateNeutralLosses(
			Ions& ions,
			double mass,
			long position,
			const std::vector<NeutralLossPair>& neutralLosses) const;
			
		virtual double fixMass(double mass) const;
};

class BIonGenerator : public SimpleIonGenerator
{
	public:
		BIonGenerator();

		~BIonGenerator() override = default;
		
	private:
		void generateRadicalIons(
			Ions& ions,
			double mass,
			long position) const override;
			
		double fixMass(double mass) const override;
};

class YIonGenerator : public SimpleIonGenerator
{
	public:
		YIonGenerator();

		~YIonGenerator() override = default;
		
	private:
			
		double fixMass(double mass) const override;
};

class AIonGenerator : public SimpleIonGenerator
{
	public:
		AIonGenerator();

		~AIonGenerator() override = default;
		
	private:
		void generateRadicalIons(
			Ions& ions,
			double mass,
			long position) const override;
			
		double fixMass(double mass) const override;
};

class CIonGenerator : public SimpleIonGenerator
{
	public:
		CIonGenerator();

		~CIonGenerator() override = default;
		
	private:
		void generateRadicalIons(
			Ions& ions,
			double mass,
			long position) const override;
			
		double fixMass(double mass) const override;
};

class ZIonGenerator : public SimpleIonGenerator
{
	public:
		ZIonGenerator();

		~ZIonGenerator() override = default;
		
	private:
		void generateRadicalIons(
			Ions& ions,
			double mass,
			long position) const override;
			
		double fixMass(double mass) const override;
};

class XIonGenerator : public SimpleIonGenerator
{
	public:
		XIonGenerator();

		~XIonGenerator() override = default;

	private:
		void generateRadicalIons(
			Ions& ions,
			double mass,
			long position) const override;

		double fixMass(double mass) const override;
};

class ImmoniumIonGenerator : public SimpleIonGenerator
{
	public:
		ImmoniumIonGenerator();

		~ImmoniumIonGenerator() override = default;
		
	private:
		std::pair<int, int> preProcessMasses(
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

		~PrecursorIonGenerator() override = default;
		
		Ions generate(
			const std::vector<double>& masses,
			long charge,
			const std::vector<NeutralLossPair>& neutralLosses,
			bool radical,
			const std::string& sequence) const override;
};

void chargeIons(const Ions& sourceIons, Ions& target, long chargeState);

inline Ion generateNeutralLossIon(
	const std::string& typeChar,
	const NeutralLossPair neutralLoss,
	double mass,
	long position)
{
	return {
		mass - neutralLoss.second,
		"[" + typeChar + std::to_string(position + 1) + "-" + neutralLoss.first + "][+]",
		position + 1
	};
}

void mergeIonVectors(Ions& target, const Ions& source);

#endif // _PEPFRAG_IONGENERATOR_H
