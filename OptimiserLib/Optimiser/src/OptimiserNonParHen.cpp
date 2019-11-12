// OptimiserNonParHen.cpp
// Bayesian quasi - Newton method with nonparametric prior for the Hessian
//
// Bayesian quasi - Newton method with a constant model for the Hessian
// (i.e.it assumes the Hessian is the same everywhere).This is a
// generalization of classic quasi - Newton methods, such as DFP, BFGS and
// SR1.
//
// 2016
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif

//////////////////////////////////////////////////

#include "OptimiserNonParHen.h"
#include "optimisable.h"
#include <cassert>

//////////////////////////////////////////////////

void COptimiserNonParHen::Optimise (double target)
{
	assert (m_ptrOptimisableFunction);

	const auto iDim = m_ptrOptimisableFunction->GetFunctionDimension ();
	auto mEye = arma::mat (iDim, iDim).eye ();
	auto mHess = mEye;
	auto mVar = mEye;

	auto vGrad = m_ptrOptimisableFunction->CalculateGradient (m_MinVector);

	auto rAlpha = 1.0e-1; // this should satisfy the Wolfe conditions

	auto vVOld = (rAlpha*m_MinVector).eval();

	do
	{
		arma::vec vSigma = m_MinVector - vVOld;//-rAlpha*(mHess*vGrad).eval ();

		arma::vec vXip1 = (m_MinVector + vSigma).eval ();

		auto vGrad = m_ptrOptimisableFunction->CalculateGradient (m_MinVector);
		auto vGradp1 = m_ptrOptimisableFunction->CalculateGradient (vXip1);

		auto rY = (vGradp1 - vGrad).eval ();

		auto vYBS = (rY - mHess*vSigma).eval();
		auto vVS = (mVar*vSigma).eval();
		auto rSYBS = (vSigma.t()*vYBS).eval()[0];

		//SVSB = s' * VS + b; // denom
		auto rSVSB = (vSigma.t()*vVS).eval()[0]; // supposed to add noise
		assert (rSVSB != 0.0);
		auto mVSSV = (vVS*arma::trans (vVS) / rSVSB).eval();

		auto mT1p = (vYBS * vVS.t()).eval();
		auto mT1 = ( (mT1p + mT1p.t()) / rSVSB).eval();
		auto mT2 = (mVSSV * (rSYBS / rSVSB)).eval();

		mHess = (mHess + mT1 - mT2).eval();
		mHess = (0.5*(mHess + mHess.t())).eval(); // mHess should be symmetric so this step is only cautionary
		mVar = (mVar - mVSSV).eval();
		auto mIHess = (arma::inv(mHess)).eval();
		vVOld = m_MinVector;
		m_MinVector = (m_MinVector - mIHess*vGrad).eval ();
	} while ((m_OOFValue = m_ptrOptimisableFunction->Evaluate (m_MinVector)) - target > m_rTolerance);

	return;
}

//////////////////////////////////////////////////
// Second method, move to new file
/*
function[x, fx, i, outs] = nonparametric_Gauss (x0, fx0, dfx0, opts) % #ok<*DEFNU>
% set up prior
i = opts.length < 0; % initialize resource counter
	% run variables
	x = x0; fx = fx0; r = -dfx0; u = -r'*r; c = -1/u; ok = 0; % linear search params
	E = r / norm (r); X0 = x; % unit line search storage
	B = opts.B0; b = max (opts.noise * opts.noise, 1.0e-8);

X0 = x; S = []; Y = []; % storage for past observations
sKs = []; U = []; NewtonFlag = 0;
while i < abs (opts.length)
	[x, c, fx0, dfx, i, A, F, dF, allEvals] = lineSearch (x0, fx0, dfx0, r, u, c, i, opts);
outs.AllEvals.x = [outs.AllEvals.x, allEvals.x];
outs.AllEvals.f = [outs.AllEvals.f, allEvals.f];
%A = A (end); F = F (end); dF = F (:, end);     % test: use only last observation
if i < 0
	i = -i; if ok, ok = 0; else break; end;   % try steepest descent or stop
else
if H >= 0
LS{ min (t,H + 1) }.a = A; LS{ min (t,H + 1) }.F = F; LS{ min (t,H + 1) }.dF = dF;
LS{ min (t,H + 1) }.b = [0; A (1:end - 1)];
%     elseif length (A) > 1
% x0 = x0 + A (end - 1) * E (:, end); X0 (:, end) = x0;
%         LS{ t }.a = A (end) - A (end - 1);
%         LS{ t }.b = 0;
%         LS{ t }.F = F (end); LS{ t }.dF = dF (:, end); dfx0 = dF (:, end - 1);
else
LS{ t }.a = A; LS{ t }.b = [0; A (1:end - 1)]; LS{ t }.F = F; LS{ t }.dF = dF;
%if t > 1; LLS0{ tt0 }.a = A; LLS0{ tt0 }.b = [0; A (1:end - 1)]; end
end;
if all (dfx == 0); ok = 1; break; end; % we are at the local minimum.
ok = 1;
Delta = LS{ t }.a - LS{ t }.b; % A - [0; A (1:end - 1)];
s = E (:, end) * Delta'; y = LS{t}.dF - [dfx0, LS{t}.dF(:,1:end-1)]; %y = dfx - dfx0;
Y = [Y, y]; S = [S, s]; %#ok<AGROW>
%gamma = 1; % mean (diag (y' * y ./ (y'*s)))
% B0 = gamma * eye (N); H0 = eye (N) / gamma;
[sKs, U] = SchurUpdate (sKs .* (U*U'),E,X0,LS,hyp);% update Schur complement of covariance
	%if t > 1; [sks1, U1] = SchurUpdate (ssKs0 .* (UU0*UU0'),EE0,XX0,LLS0,hyp); end
		ks = Projection (x, E, X0, LS, hyp, U);  % update kernel projection
		%EE0 = E; XX0 = X0; LLS0 = LS; ssKs0 = sKs; UU0 = U; tt0 = t;
[sKs, ks, U, E, X0, LS, Y, S, t] = SlimGramMatrix (1e-2, 100, sKs, ks, U, E, X0, LS, Y, S, t); % ensure numerical stability
try
L = chol (sKs + b * eye (size (sKs, 1)));
catch
fprintf 'reached numerical criticality. Switching to BFGS'
p = opts; p.H0 = inv (B); p.length = p.length - sign (p.length) * i;
outs.bfgsswitch = t;
[xbfgs, fxbfgs, ibfgs, outs_bfgs] = BFGS (x, fx0, dfx, p);
fx = [fx; fxbfgs];
x = xbfgs;
return;
end
YBS = bsxfun (@rdivide,Y - B0 * S, U'); % numerical conditioning
	SYBS = bsxfun (@ldivide,U, S' * YBS);   % numerical conditioning
		ksL = ks / L;
T1 = YBS / L * ksL'; T1 = T1 + T1';
T2 = ksL * (L' \ SYBS) / L * ksL';
B = B0 + T1 - T2;
B = 0.5 * (B + B'); % symmetrize, because the method does not explicitly encode this.
	t = t + 1;
%         if size (sKs, 1) > 20; % only if the posterior is beginning to settle
%             % adjust hyperparameters :
%             v = ones (N, 1); % power iteration
%             for i = 1:10;
%                 q = v / norm (v);
%                 v = abs (B) * q;
%                 lambda = q' * v;
%             end
%             hyp.V = 0.9 * hyp.V + 0.1 * diag (lambda * q);
%             %keyboard;
%         end
end
if ok;
Hnew = InverseHessian (ksL, L, YBS, SYBS, H0);
%keyboard;
% r = -B \ dfx;
r = -Hnew * dfx;  % Newton direction.Analytic inverse to be implemented soon.
else r = -dfx; % line search failed.Try gradient - descent.
end
if r' * dfx > 0; r = -r; NewtonFlag = [NewtonFlag; 0]; else NewtonFlag = [NewtonFlag; 1]; end % heuristic way of avoiding indefinite updates
%r = r + 0.001 * randn (N, 1);
% keyboard;
u = r' * dfx; x0 = x; dfx0= dfx; fx = [fx;fx0];  %#ok<AGROW>
outs.NewtonFlag = NewtonFlag;
outs.Y = Y; outs.S = S;
if ok
outs.E = E; outs.X0 = X0; outs.LS = LS; outs.U = U; outs.sKs = sKs;
E = [E, r / norm (r)]; X0 = [X0, x];
end
%EE0 = [EE0, r / norm (r)]; XX0 = [XX0, x]; tt0 = tt0 + 1;
outs.B{ t } = B; outs.X = [outs.X, x]; outs.V{ t } = hyp.V;
% removing old observations from the pipeline :
if H >= 0 && t > H + 1 % incorporate old observations into mean
B0 = outs.B{ t - H }; % move mean to posterior mean from H observations ago
E = E (:, end - H : end); % remove old linesearches from storage
X0 = X0 (:, end - H : end);
if H > 0;
ls = cell (H, 1); for h = 1:H; ls{ h } = LS{ end - H + h }; end; LS = ls;
mt = 0; for h = 1:H; mt = mt + length (LS{ h }.a); end;
sKs = sKs (end - mt + 1:end, end - mt + 1 : end);
S = S (:, end - mt + 1 : end);
Y = Y (:, end - mt + 1 : end);
U = U (end - mt + 1:end);
else LS = {}; sKs = []; S = []; Y = []; U = []; end;
end
end
*/
