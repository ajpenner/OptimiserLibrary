// optimiser.h
//
// 2016
//////////////////////////////////////////////////

#pragma once
// using once shows this error once per header 
#ifdef _WIN32
#pragma warning( disable : 4251 ) // disable annoying dll-interface warning due to using armadillo (same will happen with STL)
#endif

//////////////////////////////////////////////////

#include <armadillo>
#include <boost/variant.hpp>

//////////////////////////////////////////////////

#ifdef _WIN32
#include "../OptimiserDLL/Exporter.h"
#else
#define OPTIMISER
#endif

class IFunction;
// want to remove the arma::vec, this is too much
typedef boost::variant<uint32_t, int, double, arma::vec, std::pair<double,double>> variant;

//////////////////////////////////////////////////

class OPTIMISER IOptimiser
{
public:
	virtual void Optimise(double target = 0) = 0;
	virtual void SetOptimisableFunction(std::shared_ptr<const IFunction>ptrOptimisableFunction) = 0;
	virtual double GetOOF() const = 0;
	virtual void SetTolerance (double rTolerance) = 0;
	virtual arma::vec GetVector() const = 0;
	virtual void SetInitialGuess (const arma::vec& vGuess) = 0;
	virtual void SetMaximumIterations (uint32_t max) = 0;

	virtual double GetInverseOOF () const { return std::numeric_limits<double>::max();  };
	virtual arma::vec GetInverseVector () const { return arma::vec (); };
	virtual void InverseOptimise (double target = 0) {};

	virtual void SetLeftBracket(const arma::vec& v) = 0;
	virtual void SetRightBracket(const arma::vec& v) = 0;

	virtual void SetParameters (const std::vector<variant>& parameters) = 0;
	virtual ~IOptimiser() {};
};

//////////////////////////////////////////////////

class OPTIMISER BOptimiser : public IOptimiser
{
public:
	BOptimiser () : 
		m_ptrOptimisableFunction (nullptr), 
		m_MaxIteration (100), 
		m_OOFValue (std::numeric_limits<double>::max()), 
		m_rTolerance (1e-3), 
		m_MinVector (arma::vec({0})) 
	{}

	BOptimiser (std::shared_ptr<const IFunction> ptrOptimiseable) : 
		m_ptrOptimisableFunction (ptrOptimiseable), 
		m_MaxIteration (100),
		m_OOFValue (std::numeric_limits<double>::max()),
		m_rTolerance (1e-3),
		m_MinVector (arma::vec ({ 0 })) 
	{}

	void SetOptimisableFunction (std::shared_ptr<const IFunction> ptrOptimisableFunction) override
	{
		m_ptrOptimisableFunction = ptrOptimisableFunction;
	}

	double GetOOF () const override
	{
		return m_OOFValue;
	}

	void SetTolerance (double rTolerance) override
	{
		m_rTolerance = rTolerance;
	}

	arma::vec GetVector () const override
	{
		return m_MinVector;
	}

	void SetInitialGuess (const arma::vec& vGuess) override
	{
		m_MinVector = vGuess;
	}

	void SetMaximumIterations (uint32_t max) override 
	{
		m_MaxIteration = max;
	}

	void SetLeftBracket (const arma::vec& v) override
	{
		m_LeftBracket = v;
	}

	void SetRightBracket (const arma::vec& v) override
	{
		m_RightBracket = v;
	}

	virtual ~BOptimiser () {};

protected:
	std::shared_ptr<const IFunction> m_ptrOptimisableFunction;
	uint32_t m_MaxIteration;
	double m_OOFValue;
	double m_rTolerance;
	arma::vec m_MinVector;

	arma::vec m_LeftBracket;
	arma::vec m_RightBracket;
};

//////////////////////////////////////////////////

class OPTIMISER IExtremum
{
public:
	virtual void SetDirectionVector (const arma::vec& vector) = 0;
	virtual void SetLeftBracketRange (double rLeft) = 0;
	virtual void SetRightBracketRange (double rRight) = 0;

	virtual ~IExtremum () {};
};

//////////////////////////////////////////////////

class OPTIMISER BExtremum : public BOptimiser, public IExtremum
{
public:
	BExtremum () {}

	void SetDirectionVector (const arma::vec& vector) override
	{
		BOptimiser::SetInitialGuess(vector);
	}

	void SetLeftBracketRange (double rLeft) override 
	{ 
		m_LeftBracket = (rLeft*m_MinVector).eval(); 
	}

	void SetRightBracketRange (double rRight) override 
	{ 
		m_RightBracket = (rRight*m_MinVector).eval (); 
	}

	virtual ~BExtremum () {};

private:

};

//////////////////////////////////////////////////

class OPTIMISER BBracket : public BOptimiser
{
public:
	virtual arma::vec GetRange() const = 0;

	virtual ~BBracket () {};
};

//////////////////////////////////////////////////

class IConstraint;

class OPTIMISER IConstraintOptimiser
{
public:
	virtual void SetMaxIterations(size_t iCount) = 0;
	virtual void AddConstraint( std::shared_ptr<const IConstraint> ptrConstraint) = 0;
	virtual void SetUnconstrainedOptimiser (std::unique_ptr<IOptimiser> ptrUnconstrainedOptimiser) = 0;

	virtual ~IConstraintOptimiser () {};
};

//////////////////////////////////////////////////

class OPTIMISER BConstraintOptimiser : public IConstraintOptimiser, public BOptimiser
{
public:
	virtual void SetMaxIterations (size_t iCount) 
	{ 
		m_uMaxIterations = iCount; 
	}

	virtual void AddConstraint (std::shared_ptr<const IConstraint> ptrConstraint) = 0;
	
	virtual void SetUnconstrainedOptimiser (std::unique_ptr<IOptimiser> ptrUnconstrainedOptimiser)
	{
		m_ptrUnconstrainedOptimiser.reset( ptrUnconstrainedOptimiser.release() );
	}

	virtual ~BConstraintOptimiser () {};

protected:

	size_t m_uMaxIterations;
	std::unique_ptr<IOptimiser> m_ptrUnconstrainedOptimiser;
};

//////////////////////////////////////////////////
