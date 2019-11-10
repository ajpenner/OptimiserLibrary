// OptimiserHenning.cpp
// Bayesian quasi - Newton method with nonparametric prior for the Hessian
//
// Bayesian quasi - Newton method with a constant model for the Hessian
// (i.e.it assumes the Hessian is the same everywhere).This is a
// generalization of classic quasi - Newton methods, such as DFP, BFGS and
// SR1.
//
// 2016
//////////////////////////////////////////////////

#include "stdafx.h"

//////////////////////////////////////////////////

#include "OptimiserHennig.h"
#include "optimisable.h"

#include "LinearFunction.h"
#include "MinimiserGolden.h"
#include "BracketBoundingPhase.h"

#include <cassert>
#include <iostream>

//////////////////////////////////////////////////

double COptimiserHennig::GetMinimalScaleValue (std::shared_ptr<CLinearFunction> ptrVF, CMinimiserGolden& minimiser, const arma::vec& direction)
{
	arma::vec guess = arma::vec ({ 1 });
	ptrVF->SetVector (direction);
	ptrVF->SetOffset (m_MinVector);

	minimiser.SetOptimisableFunction (ptrVF);

	CBracketBoundingPhase bp;
	std::vector<variant> delta = { 0.5 };
	bp.SetParameters (delta);
	bp.SetOptimisableFunction (ptrVF);
	bp.SetInitialGuess (arma::vec{ 0 });
	bp.Optimise (0);
	auto range = bp.GetRange ();

	minimiser.SetInitialGuess (guess);
	minimiser.SetLeftBracket (arma::vec ({ range.at (0) }));
	minimiser.SetRightBracket (arma::vec ({ range.at (1) }));

	minimiser.Optimise (0);

	return minimiser.GetVector ().at (0);
}

//////////////////////////////////////////////////

arma::mat COptimiserHennig::UpdateHessian()
{
	arma::vec VarS = m_mVar*m_deltaX;

	arma::mat rDenom = m_deltaX.t()*VarS; // may need to add noise
	assert (rDenom.at(0) != 0.0);
	auto rho = 1.0 / rDenom.at(0);

	arma::vec V = m_deltaGrad - m_mHess*m_deltaX;

	arma::mat mP1 = V*VarS.t()*rho;
	arma::mat mP2 = VarS*(m_deltaX.t()*V)*VarS.t()*rho*rho;

	arma::mat thing = VarS*VarS.t()*rho;
	m_mVar -= thing; // after two iterations m_mVar is ALWAYS zero - identically
	return m_mHess + mP1 + mP1.t() - mP2;
}

//////////////////////////////////////////////////

void COptimiserHennig::Optimise (double target)
{
	assert (m_ptrOptimisableFunction);

	const auto iDim = m_ptrOptimisableFunction->GetFunctionDimension ();
	m_mEye = arma::eye (iDim, iDim);
	m_mHess = m_mEye;
	arma::mat mIHess = m_mEye;
	// if I fix m_mVar to identity we recover the PSB method (as expected)
	// if I fix m_mVar to m_mHess we recover the DFP method (as expected)
	m_mVar = 1e4*m_mEye; // why the large constant? The larger the constant the less the original correlation
	m_rNoise = 1e-12;

	auto vf = std::make_shared<CLinearFunction> (m_ptrOptimisableFunction);
	CMinimiserGolden minimiser;
	minimiser.SetTolerance (m_rTolerance);
	auto vGrad = m_ptrOptimisableFunction->CalculateGradient (m_MinVector);
	do
	{
		arma::vec s = -arma::normalise (mIHess*vGrad);
		auto rLambda = GetMinimalScaleValue (vf, minimiser, s);

		m_deltaX = rLambda*s;
		arma::vec newMin = m_MinVector + m_deltaX;
		arma::vec vGradp1 = m_ptrOptimisableFunction->CalculateGradient (newMin);
		m_deltaGrad = vGradp1 - vGrad;

		m_mHess = UpdateHessian();
		mIHess = arma::inv(m_mHess);

		m_MinVector = newMin;
		vGrad = vGradp1;
		m_OOFValue = m_ptrOptimisableFunction->Evaluate (m_MinVector);
	} while ( m_OOFValue - target > m_rTolerance);

	return;
}

//////////////////////////////////////////////////

/*
function [x,fx,i,outs] = constant_Gauss(x0,fx0,dfx0,opts)
% Bayesian quasi-Newton method with a constant model for the Hessian
% (i.e. it assumes the Hessian is the same everywhere). This is a
% generalization of classic quasi-Newton methods, such as DFP, BFGS and
% SR1.
if ~isfield(opts,'SIG'), opts.SIG = 0.5; end % default for line search quality
if ~isfield(opts,'noise'), opts.noise = 0; end
% set up prior:
N = length(x0);
if ~isfield(opts,'B0'), opts.B0 = eye(N); end
if ~isfield(opts,'V0'), opts.V0 = (100 * 100) * eye(N); end
i = opts.length < 0; % initialize resource counter
% run variables
x = x0; fx = fx0; r = -dfx0; u = -r'*r; c = -1/u; ok = 0; % linear search params
B = opts.B0; % prior mean
V = opts.V0; % prior covariance
b = opts.noise * opts.noise; % noise variance
outs.B{1} = B; outs.V{1} = V; outs.X = x;
t = 1; outs.AllEvals.x = []; outs.AllEvals.f = [];
while i < abs(opts.length)
[x,c,fx0,dfx,i, ~, ~, ~, allEvals] = lineSearch(x0,fx0,dfx0,r,u,c,i,opts);
outs.AllEvals.x = [outs.AllEvals.x, allEvals.x];
outs.AllEvals.f = [outs.AllEvals.f, allEvals.f];
if i < 0
i = -i; if ok, ok = 0; else break; end;   % try steepest descent or stop
else
if all(dfx == 0); break; end; % we are at the local minimum.
ok = 1; s = x - x0; y = dfx - dfx0;
outs.y{t} = y; outs.s{t} = s;
% Bayesian regression update:
YBS = y - B * s; VS = V * s; SYBS = s' * YBS; % VS = V * s
% check noise required for positive definiteness
%a = y' * (B \ VS); d = (y' * (B \ YBS)) * (VS' * (B \ VS));
%b1 = -a/2 + sqrt(d-3*a); b2 = -a/2 - sqrt(d-3*a);
%if ~isreal(b1); b1 = 0; b2 = 0; end
%bt = max([0,b1,b2,b]) + 0.5;
SVSB = s' * VS + b;
VSSV = (VS * VS') / SVSB;
T1  = YBS * VS'; T1 = (T1 + T1') / SVSB;
T2  = VSSV * (SYBS / SVSB);
B   = B + T1 - T2;
B   = 0.5 * (B + B');
V = V - VSSV;
end
r = - B \ dfx;  % Newton direction. Analytic inverse to be implemented soon.
if r' * dfx > 0; r = -r;  end
u = r' * dfx; x0 = x; dfx0= dfx; fx = [fx;fx0];  %#ok<AGROW>
t = t + 1;
outs.B{t} = B;
outs.V{t} = V;
outs.X    = [outs.X, x];
end

*/