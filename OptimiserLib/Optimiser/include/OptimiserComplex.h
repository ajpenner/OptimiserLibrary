// OptimiserComplex.h
// Simplex method with contraints
// 2017
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "optimiser.h"
#include "optimisable.h"
#include "constraints.h"
#include "CompoundMOM.h"

#include <limits>
#include <armadillo>
#include <random>

//////////////////////////////////////////////////

namespace eComplexParameters
{
	enum type
	{
		eBeta = 0,
		eGamma,
		eRange
	};
}

namespace eComplexCoordinates
{
	enum type
	{
		eCoordinate = 0,
		eEvaluation
	};
}


class OPTIMISER COptimiserComplex : public BConstraintOptimiser
{
public:
	COptimiserComplex () : m_simplexCount(0), m_dim(0)
	{
		m_ptrGenerator = std::make_unique<std::mt19937>((*m_ptrRandDev)());
		m_ptrDistribution = std::make_unique<std::uniform_real_distribution<>> (0, 1);
	}

	~COptimiserComplex () {}

	void SetParameters (const std::vector<variant>& parameters) override;

	void Optimise (double target) override;

	void AddConstraint (std::shared_ptr<const IConstraint> ptrConstraint) override;

private:
	void SetInitialPoints();
	bool VerifyConstraints(const arma::vec& vec);
	void EnforceRange(arma::vec& vec);
	void EnforceConstraints(arma::vec& vec);

	std::vector<double> m_Multipliers;
	std::vector<std::shared_ptr<const IConstraint>> m_vConstraints;
	std::unique_ptr<std::random_device> m_ptrRandDev;
	std::unique_ptr<std::mt19937> m_ptrGenerator;
	std::unique_ptr<std::uniform_real_distribution<>> m_ptrDistribution;

	std::vector<arma::vec> m_vCoordinates;
	std::vector<double> m_evaluations;
	std::vector<std::pair<double, double>> m_vRange;

	uint8_t m_simplexCount;
	uint8_t m_dim;
	double m_gamma;
	double m_beta;
};

//////////////////////////////////////////////////