// OptimiserRandom.h
//
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

namespace eRandomParameters
{
	enum type
	{
		eScale = 0,
		ePointCount,
		eRange
	};
}

class OPTIMISER COptimiserRandom : public BConstraintOptimiser
{
public:
	COptimiserRandom () : m_iPointCount(0), m_dim(0)
	{
		m_ptrGenerator = std::make_unique<std::mt19937>((*m_ptrRandDev)());
		m_ptrDistribution = std::make_unique<std::uniform_real_distribution<>> (-0.5, 5);
	}

	~COptimiserRandom () {}

	void SetParameters (const std::vector<variant>& parameters) override;

	void Optimise (double target) override;

	void AddConstraint (std::shared_ptr<const IConstraint> ptrConstraint) override;

private:
	void SetRandomPoints();
	bool VerifyConstraints(const arma::vec& vec);
	void EnforceRange(arma::vec& vec);

	double m_rScale;
	size_t m_iPointCount;
	std::vector<double> m_Multipliers;
	std::vector<std::shared_ptr<const IConstraint>> m_vConstraints;
	std::unique_ptr<std::random_device> m_ptrRandDev;
	std::unique_ptr<std::mt19937> m_ptrGenerator;
	std::unique_ptr<std::uniform_real_distribution<>> m_ptrDistribution;

	std::vector<arma::vec> m_vCoordinates;
	std::vector<double> m_evaluations;
	std::vector<std::pair<double, double>> m_vRange;

	uint8_t m_dim;
};

//////////////////////////////////////////////////