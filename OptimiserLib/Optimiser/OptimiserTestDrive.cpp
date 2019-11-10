// OptimiserTestDrive.cpp : Defines the entry point for the console application.
//
// 2017
//////////////////////////////////////////////////

#include "stdafx.h"

//////////////////////////////////////////////////

#include "LinearConstraintImpl.h"
#include "NonlinearConstraintImpl.h"

#include "optimisable.h"

// test functions
#include "Rosenbrock.h"
#include "Parabolic.h"
#include "Himmelblau.h"
#include "Deb.h"
#include "FirstConstraint.h"
#include "SecondConstraint.h"
#include "NegativeFirstConstraint.h"

// Factories
#include "DirectionFactory.h"
#include "LineSearchFactory.h"
#include "OptimiserFactory.h"
#include "OptimiserFactoryVM.h"
#include "ExtremumFactory.h"

// Bracketing
#include "BracketBoundingPhase.h"

#include <vector>
#include <iostream>
#include <iomanip>
#include <armadillo>
#include <functional>

//////////////////////////////////////////////////

int main ( int argc, char* argv[] )
{	
	auto ptrParabolic = std::make_shared<CParabolic>(1.9,2.5,-6.65);

	std::unique_ptr<IOptimiser> ptrOptimiserBrent = OptimiserFactory::make_optimiser(eOptimiser::eBrent);
	ptrOptimiserBrent->SetLeftBracket (arma::vec ({-2}));
	ptrOptimiserBrent->SetRightBracket (arma::vec ({4}));
	ptrOptimiserBrent->SetOptimisableFunction(ptrParabolic);
	ptrOptimiserBrent->SetTolerance (1e-8);
	ptrOptimiserBrent->Optimise();
	auto minVec = ptrOptimiserBrent->GetVector ();

	std::unique_ptr<BExtremum> ptrMinimiserBrent = ExtremumFactory::make_extremum(eExtremum::eMinBrent);
	ptrMinimiserBrent->SetDirectionVector(arma::vec ({1}));
	ptrMinimiserBrent->SetLeftBracketRange(-2);
	ptrMinimiserBrent->SetRightBracketRange(4);
	ptrMinimiserBrent->SetOptimisableFunction (ptrParabolic);
	ptrMinimiserBrent->SetTolerance (1e-8);
	ptrMinimiserBrent->Optimise ();
	auto minVec2 = ptrMinimiserBrent->GetVector ();

	auto ptrDeb = std::make_shared<CDeb>(54);
	std::unique_ptr<BExtremum> ptrMinimiserGolden = ExtremumFactory::make_extremum (eExtremum::eGolden);
	ptrMinimiserGolden->SetDirectionVector (arma::vec ({1}));
	ptrMinimiserGolden->SetLeftBracketRange (0);
	ptrMinimiserGolden->SetRightBracketRange (5);
	ptrMinimiserGolden->SetOptimisableFunction (ptrDeb);
	ptrMinimiserGolden->SetTolerance (1e-14);
	ptrMinimiserGolden->Optimise ();
	auto minVec3 = ptrMinimiserGolden->GetVector ();

	CBracketBoundingPhase boundingPhase;
	boundingPhase.SetInitialGuess(arma::vec ({0.6}));
	boundingPhase.SetOptimisableFunction (ptrDeb);
	std::vector<variant> delta = {0.001};
	boundingPhase.SetParameters(delta);
	boundingPhase.Optimise (0);
	auto range = boundingPhase.GetRange();
	
	std::unique_ptr<BExtremum> ptrMinimiserBisection = ExtremumFactory::make_extremum (eExtremum::eBisection);
	ptrMinimiserBisection->SetDirectionVector (arma::vec ({1}));
	ptrMinimiserBisection->SetLeftBracketRange (-2);
	ptrMinimiserBisection->SetRightBracketRange (4);
	ptrMinimiserBisection->SetOptimisableFunction (ptrParabolic);
	ptrMinimiserBisection->SetTolerance (1e-12);
	ptrMinimiserBisection->Optimise ();
	auto minVec4 = ptrMinimiserBisection->GetVector ();
	
	std::unique_ptr<IDirection> ptrDirection = DirectionFactory::make_direction(eDirection::eCauchy);
	std::unique_ptr<ILineSearch> ptrLineSearchMethod = LineSearchFactory::make_lineSearch(eLineSearch::eStrongWolfe);
	ptrLineSearchMethod->SetDirectionMethod(ptrDirection.release());
	ptrLineSearchMethod->SetTolerance (1e-10);

	auto ptrRosenbrock = std::make_shared<CRosenbrock>(2, 4); // A = 2, B = 4
	std::unique_ptr<IOptimiser> ptrOptimiserNewton = OptimiserFactory::make_optimiser (eOptimiser::eNewton);
	ptrOptimiserNewton->SetOptimisableFunction (ptrRosenbrock);
	ptrOptimiserNewton->SetTolerance (1e-9);
	ptrOptimiserNewton->SetInitialGuess(arma::vec ({ 1,2 }));
	ptrOptimiserNewton->Optimise();

	auto ptrHimmelblau = std::make_shared<CHimmelblau>(-11, -7);

	std::unique_ptr<IOptimiser> ptrOptimiserBox = OptimiserFactory::make_optimiser (eOptimiser::eBox);
	ptrOptimiserBox->SetOptimisableFunction (ptrHimmelblau);
	ptrOptimiserBox->SetTolerance (1e-9);
	ptrOptimiserBox->SetLeftBracket (arma::vec ({0,0}));
	ptrOptimiserBox->SetRightBracket (arma::vec ({2,2}));
	ptrOptimiserBox->Optimise ();

	std::unique_ptr<IOptimiser> ptrOptimiserSimplex = OptimiserFactory::make_optimiser (eOptimiser::eSimplex);
	ptrOptimiserSimplex->SetOptimisableFunction (ptrHimmelblau);
	ptrOptimiserSimplex->SetTolerance (1e-9);
	std::vector<variant> simplex = { 0.5, 1.5, arma::vec ({ 0,0 }), arma::vec ({ 2,0 }), arma::vec ({ 1,1 }) };
	ptrOptimiserSimplex->SetParameters(simplex);
	ptrOptimiserSimplex->SetLeftBracket (arma::vec ({ 0,0 }));
	ptrOptimiserSimplex->SetRightBracket (arma::vec ({ 2,2 }));
	ptrOptimiserSimplex->Optimise ();

	std::unique_ptr<IOptimiser> ptrOptimiserHJ = OptimiserFactory::make_optimiser (eOptimiser::eHJ);
	ptrOptimiserHJ->SetOptimisableFunction (ptrHimmelblau);
	ptrOptimiserHJ->SetTolerance (1e-9);
	ptrOptimiserHJ->SetInitialGuess (arma::vec ({0,0}));
	ptrOptimiserHJ->Optimise ();
	
	// bounding phase needs to be handled carefully, sloppy tolerances are good
	std::unique_ptr<IOptimiser> ptrOptimiserPowell = OptimiserFactory::make_optimiser (eOptimiser::ePowell);
	ptrOptimiserPowell->SetOptimisableFunction(ptrHimmelblau);
	ptrOptimiserPowell->SetTolerance (1e-9);
	ptrOptimiserPowell->SetInitialGuess(arma::vec({0,4}));
	ptrOptimiserPowell->Optimise ();

	std::unique_ptr<IOptimiser> ptrOptimiserCauchy = OptimiserFactory::make_optimiser (eOptimiser::eCauchy);
	ptrOptimiserCauchy->SetOptimisableFunction (ptrHimmelblau);
	ptrOptimiserCauchy->SetTolerance (1e-9);
	ptrOptimiserCauchy->SetMaximumIterations(100);
	ptrOptimiserCauchy->SetInitialGuess (arma::vec ({0,0}));
	ptrOptimiserCauchy->Optimise ();

	std::unique_ptr<IOptimiser> ptrOptimiserMarquardt = OptimiserFactory::make_optimiser (eOptimiser::eMarquardt);
	ptrOptimiserMarquardt->SetOptimisableFunction (ptrHimmelblau);
	ptrOptimiserMarquardt->SetTolerance (1e-9);
	std::vector<variant> penalty = { 100.0 };
	ptrOptimiserMarquardt->SetParameters(penalty);
	ptrOptimiserMarquardt->SetInitialGuess (arma::vec ({ 0,0 }));
	ptrOptimiserMarquardt->Optimise (0);

	std::unique_ptr<IOptimiser> ptrOptimiserCG = OptimiserFactory::make_optimiser (eOptimiser::eConjugateGradient);
	ptrOptimiserCG->SetOptimisableFunction (ptrHimmelblau);
	ptrOptimiserCG->SetTolerance (1e-9);
	ptrOptimiserCG->SetInitialGuess (arma::vec ({ 0,0 }));
	ptrOptimiserCG->Optimise (0);

	// Variable metric options
	std::unique_ptr<IOptimiserVariableMetric> ptrOptimiserVMDFP = OptimiserFactoryVM::make_optimiser (eOptimiserVM::eVMDFP);
	ptrOptimiserVMDFP->SetOptimisableFunction (ptrHimmelblau);
	ptrOptimiserVMDFP->SetTolerance (1e-9);
	ptrOptimiserVMDFP->SetInitialGuess (arma::vec ({ 0,0 }));
	ptrOptimiserVMDFP->Optimise ();

	std::unique_ptr<IOptimiserVariableMetric> ptrOptimiserVMBFGS = OptimiserFactoryVM::make_optimiser (eOptimiserVM::eVMBFGS);
	ptrOptimiserVMBFGS->SetOptimisableFunction (ptrHimmelblau);
	ptrOptimiserVMBFGS->SetTolerance (1e-9);
	ptrOptimiserVMBFGS->SetInitialGuess (arma::vec ({ 0,0 }));
	ptrOptimiserVMBFGS->Optimise ();

	std::unique_ptr<IOptimiserVariableMetric> ptrOptimiserVMSR1 = OptimiserFactoryVM::make_optimiser (eOptimiserVM::eVMSR1);
	ptrOptimiserVMSR1->SetOptimisableFunction(ptrHimmelblau);
	ptrOptimiserVMSR1->SetTolerance(1e-9);
	ptrOptimiserVMSR1->SetInitialGuess (arma::vec ({0,0}));
	ptrOptimiserVMSR1->Optimise();

	std::unique_ptr<IOptimiserVariableMetric> ptrOptimiserVMBroyden = OptimiserFactoryVM::make_optimiser (eOptimiserVM::eVMBroyden);
	ptrOptimiserVMBroyden->SetOptimisableFunction (ptrHimmelblau);
	ptrOptimiserVMBroyden->SetTolerance (1e-9);
	std::vector<variant> phi = { 0.5 };
	ptrOptimiserVMBroyden->SetParameters (phi);
	ptrOptimiserVMBroyden->SetInitialGuess (arma::vec ({ 0,0 }));
	ptrOptimiserVMBroyden->Optimise ();

	std::unique_ptr<IOptimiserVariableMetric> ptrOptimiserVMBlend = OptimiserFactoryVM::make_optimiser (eOptimiserVM::eVMBlend);
	ptrOptimiserVMBlend->SetOptimisableFunction (ptrHimmelblau);
	ptrOptimiserVMBlend->SetTolerance (1e-9);
	std::vector<variant> blendRatio = { 0.5 };
	ptrOptimiserVMBlend->SetParameters (blendRatio);
	ptrOptimiserVMBlend->SetInitialGuess (arma::vec ({ 0,0 }));
	ptrOptimiserVMBlend->Optimise ();

	std::unique_ptr<IOptimiserVariableMetric> ptrOptimiserVMPSB = OptimiserFactoryVM::make_optimiser (eOptimiserVM::eVMPSB);
	ptrOptimiserVMPSB->SetOptimisableFunction (ptrHimmelblau);
	ptrOptimiserVMPSB->SetTolerance (1e-9);
	ptrOptimiserVMPSB->SetInitialGuess (arma::vec ({ 0,0 }));
	ptrOptimiserVMPSB->Optimise ();

	// BELOW NEEDS FURTHER INVESTIGATION
/* // Honestly, this is such a mess, I have a hard time seeing how this saves anything
	std::unique_ptr<IOptimiserVariableMetric> ptrOptimiserVMLBFGS = OptimiserFactoryVM::make_optimiser (eOptimiserVM::eVMLBFGS); //INCOMPLETE
	ptrOptimiserVMLBFGS->SetOptimisableFunction (ptrHimmelblau);
	ptrOptimiserVMLBFGS->SetTolerance (1e-9);
	ptrOptimiserVMLBFGS->SetInitialGuess (arma::vec ({ 0,0 }));
	ptrOptimiserVMLBFGS->Optimise ();

	COptimiserHennig optimiser8; // I am not convinced this is working, it converges but the correlation matrix is zero in the end???
	optimiser8.SetOptimisableFunction(ptrHimmelblau);
	optimiser8.SetTolerance(1e-9);
	optimiser8.SetInitialGuess (arma::vec({0,0}));
	optimiser8.Optimise(0);
*/
	// Constrained Optimisation
	// It looks like we might want to penalise coordinate boundary violations stronger than other contraints
	std::unique_ptr<BConstraintOptimiser> ptrOptimiserPenalty = ConstraintOptimiserFactory::make_optimiser (eConstraintOptimiser::ePenalty);
	{
		ptrOptimiserPenalty->SetOptimisableFunction (ptrHimmelblau);
		{
			std::unique_ptr<IOptimiser> ptrOptimiserPowell = OptimiserFactory::make_optimiser (eOptimiser::ePowell);
			ptrOptimiserPenalty->SetUnconstrainedOptimiser (std::move (ptrOptimiserPowell));
			// OptimiserPowell is gone, referencing is no longer possible
		}
		ptrOptimiserPenalty->SetTolerance (1e-6);
		ptrOptimiserPenalty->SetMaxIterations (10);

		auto ptrFC = std::make_shared<CFirstConstraint>();
		ptrFC->SetParameters (arma::vec ({ -5, -26 }));
		auto ptrC1 = std::make_shared<CNonlinearConstraintImpl>(eGE, ptrFC);
		ptrOptimiserPenalty->AddConstraint (ptrC1);

		auto ptrBoundaryX1 = std::make_shared<CLinearConstraintImpl> (eGE, 0, arma::vec ({1,0}));
		ptrOptimiserPenalty->AddConstraint (ptrBoundaryX1);
		auto ptrBoundaryX2 = std::make_shared<CLinearConstraintImpl> (eGE, 0, arma::vec ({ 0,1 }));
		ptrOptimiserPenalty->AddConstraint (ptrBoundaryX2);

		ptrOptimiserPenalty->SetInitialGuess (arma::vec ({ 0,0 }));

		std::vector<variant> parameters = { 0.1, 10.0 };
		ptrOptimiserPenalty->SetParameters (parameters);

		ptrOptimiserPenalty->Optimise ();
	}

	std::unique_ptr<BConstraintOptimiser> ptrOptimiserMOM = ConstraintOptimiserFactory::make_optimiser (eConstraintOptimiser::ePenaltyMOM);
	{
		ptrOptimiserMOM->SetOptimisableFunction (ptrHimmelblau);
		{
			std::unique_ptr<IOptimiser> ptrOptimiserPowell = OptimiserFactory::make_optimiser (eOptimiser::ePowell);
			ptrOptimiserMOM->SetUnconstrainedOptimiser (std::move (ptrOptimiserPowell));
		}
		ptrOptimiserMOM->SetTolerance (1e-6);
		ptrOptimiserMOM->SetMaxIterations (10);

		auto ptrFC = std::make_shared<CFirstConstraint> ();
		ptrFC->SetParameters (arma::vec ({ -5, -26 }));
		auto ptrC1 = std::make_shared<CNonlinearConstraintImpl> (eGE, ptrFC);
		ptrOptimiserMOM->AddConstraint (ptrC1);

		auto ptrBoundaryX1 = std::make_shared<CLinearConstraintImpl> (eGE, 0, arma::vec ({ 1,0 }));
		ptrOptimiserMOM->AddConstraint (ptrBoundaryX1);
		auto ptrBoundaryX2 = std::make_shared<CLinearConstraintImpl> (eGE, 0, arma::vec ({ 0,1 }));
		ptrOptimiserMOM->AddConstraint (ptrBoundaryX2);

		ptrOptimiserMOM->SetInitialGuess (arma::vec ({ 0,0 }));

		// needed a strong penalty for our use of Powell method, it allowed boundary violations too liberally
		std::vector<variant> parameters = { 100.0, 10.0 };
		ptrOptimiserMOM->SetParameters (parameters);

		ptrOptimiserMOM->Optimise ();
	}

	// complex does not work yet
	std::unique_ptr<BConstraintOptimiser> ptrOptimiserComplex = ConstraintOptimiserFactory::make_optimiser (eConstraintOptimiser::eComplex);
	{
		ptrOptimiserComplex->SetOptimisableFunction (ptrHimmelblau);
		{
			std::unique_ptr<IOptimiser> ptrOptimiserPowell = OptimiserFactory::make_optimiser (eOptimiser::ePowell);
			ptrOptimiserComplex->SetUnconstrainedOptimiser (std::move (ptrOptimiserPowell));
		}
		ptrOptimiserComplex->SetTolerance (1e-6);
		ptrOptimiserComplex->SetMaxIterations (10);

		auto ptrFC = std::make_shared<CNegativeFirstConstraint> ();
		ptrFC->SetParameters (arma::vec ({ -5, 26 }));
		auto ptrC1 = std::make_shared<CNonlinearConstraintImpl> (eGE, ptrFC);
		ptrOptimiserComplex->AddConstraint (ptrC1);

		auto ptrSC = std::make_shared<CSecondConstraint> ();
		ptrSC->SetParameters (arma::vec ({ 20, -5 }));
		auto ptrC2 = std::make_shared<CNonlinearConstraintImpl> (eGE, ptrSC);
		ptrOptimiserComplex->AddConstraint (ptrC2);

		auto ptrBoundaryX1 = std::make_shared<CLinearConstraintImpl> (eGE, 0, arma::vec ({ 1,0 }));
		ptrOptimiserComplex->AddConstraint (ptrBoundaryX1);
		auto ptrBoundaryX2 = std::make_shared<CLinearConstraintImpl> (eGE, 0, arma::vec ({ 0,1 }));
		ptrOptimiserComplex->AddConstraint (ptrBoundaryX2);

		ptrOptimiserComplex->SetInitialGuess (arma::vec ({ 0,0 }));

		std::vector<variant> parameters = { 0.5, 1.3, std::make_pair(0.0, 5.0), std::make_pair(0.0, 5.0) };
		ptrOptimiserComplex->SetParameters (parameters);

		ptrOptimiserComplex->Optimise ();
	}

	// this works, and is EXPENSIVE
	std::unique_ptr<BConstraintOptimiser> ptrOptimiserRandom = ConstraintOptimiserFactory::make_optimiser (eConstraintOptimiser::eRandom);
	{
		ptrOptimiserRandom->SetOptimisableFunction (ptrHimmelblau);

		// this is a direct search, it does not use an unconstrained optimiser
		ptrOptimiserRandom->SetTolerance (1e-6);
		ptrOptimiserRandom->SetMaxIterations (10);

		auto ptrFC = std::make_shared<CNegativeFirstConstraint> ();
		ptrFC->SetParameters (arma::vec ({ -5, 26 }));
		auto ptrC1 = std::make_shared<CNonlinearConstraintImpl> (eGE, ptrFC);
		ptrOptimiserRandom->AddConstraint (ptrC1);

		auto ptrSC = std::make_shared<CSecondConstraint> ();
		ptrSC->SetParameters (arma::vec ({ 20, -5 }));
		auto ptrC2 = std::make_shared<CNonlinearConstraintImpl> (eGE, ptrSC);
		ptrOptimiserRandom->AddConstraint (ptrC2);

		auto ptrBoundaryX1 = std::make_shared<CLinearConstraintImpl> (eGE, 0, arma::vec ({ 1,0 }));
		ptrOptimiserRandom->AddConstraint (ptrBoundaryX1);
		auto ptrBoundaryX2 = std::make_shared<CLinearConstraintImpl> (eGE, 0, arma::vec ({ 0,1 }));
		ptrOptimiserRandom->AddConstraint (ptrBoundaryX2);

		ptrOptimiserRandom->SetInitialGuess (arma::vec ({ 3,3 }));

		std::vector<variant> parameters = { 0.25, 100, std::make_pair (0.0, 5.0), std::make_pair (0.0, 5.0) };
		ptrOptimiserRandom->SetParameters (parameters);

		ptrOptimiserRandom->Optimise ();
	}
/*
	std::unique_ptr<BConstraintOptimiser> ptrOptimiserEllipsoid = ConstraintOptimiserFactory::make_optimiser (eConstraintOptimiser::eEllipsoid);
	{
		ptrOptimiserEllipsoid->SetOptimisableFunction(ptrHimmelblau);
		ptrOptimiserEllipsoid->SetTolerance (1e-9);
		ptrOptimiserEllipsoid->SetMaxIterations (100);

		auto ptrC1 = std::make_shared<CLinearConstraintImpl>(eLE, -8, arma::vec ({ -1,0.2 }));
		ptrOptimiserEllipsoid->AddConstraint (ptrC1);

		auto ptrC2 = std::make_shared<CLinearConstraintImpl>(eLE, 4, arma::vec ({ 1,1 }));
		ptrOptimiserEllipsoid->AddConstraint (ptrC2);

		auto ptrC3 = std::make_shared<CLinearConstraintImpl>(eLE, 9, arma::vec ({ 0.3,-1 }));
		ptrOptimiserEllipsoid->AddConstraint (ptrC3);

		ptrOptimiserEllipsoid->SetInitialGuess (arma::vec ({0,0}));

		std::vector<variant> radius = { 13 };
		ptrOptimiserEllipsoid->SetParameters (blendRatio);

		ptrOptimiserEllipsoid->Optimise ();
	}
*/
/*
	IConstraintOptimiser* pOptimiserGRG = ConstraintOptimiserFactory::make_optimiser (eConstraintOptimiser::eGRG).get ();
	{
		CHimmelblau function(-11.0, -7.0);

		pOptimiserGRG->SetOptimisableFunction (&function);
		pOptimiserGRG->SetTolerance (1e-9);
		pOptimiserGRG->SetMaxIterations (100);

//		auto constraint1 = [](const arma::vec& vCoordinates) { return 26.0 - std::pow(vCoordinates[0]-5.0,2) - std::pow(vCoordinates[1],2); };
//		CNonlinearConstraintImpl<decltype(constraint1), decltype(nullptr)> C1 (eGE, constraint1, nullptr, -8);
		auto constraint1 = CHimmelblauConstraint1 ();
		CNonlinearConstraintImpl C1 (eGE, &constraint1, -8);
		pOptimiserGRG->AddConstraint(&C1);

//		auto constraint2 = [](const arma::vec& vCoordinates) { return 20.0 - 4.0*vCoordinates[0] - vCoordinates[1]; };
//		CNonlinearConstraintImpl<decltype(constraint2), decltype(nullptr)> C2 (eGE, constraint2, nullptr, 4);
		auto constraint2 = CHimmelblauConstraint2();
		CNonlinearConstraintImpl C2 (eGE, &constraint2, 4);
		pOptimiserGRG->AddConstraint (&C2);

		auto vGuess8 = arma::vec (2);
		vGuess8 (0) = 1.;
		vGuess8 (1) = 2.;
		pOptimiserGRG->SetInitialGuess (vGuess8);

		//set range for free variables
		auto vBoundX1 = arma::vec(2);
		vBoundX1(0) = 0.0;
		vBoundX1(1) = 5.0;
		std::vector<arma::vec> vBounds;
		vBounds.push_back (vBoundX1);
		vBounds.push_back (vBoundX1); // simplified for example
		auto pConcreteGRG = static_cast<COptimiserGRG*>(pOptimiserGRG);
		pConcreteGRG->SetBounds(vBounds);

		pOptimiserGRG->Optimise ();
	}

	std::cout << std::setprecision(18) << pOptimiserEllipsoid->GetOOF () << std::endl;
	*/
}

//////////////////////////////////////////////////