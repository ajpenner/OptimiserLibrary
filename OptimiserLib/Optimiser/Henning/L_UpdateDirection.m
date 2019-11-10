function r = L_UpdateDirection(k_cholK,cholK,YBS,SYBS,H0,dfx)
% analytic inverse of the Hessian -- limited memory version.
%
% Martin Kiefel, May 2012

[~,M] = size(k_cholK);

SYBS = (SYBS + SYBS');

V = k_cholK / cholK';
W = YBS - 0.25 .* V * SYBS;

VH = bsxfun(@times,V',H0');
WH = bsxfun(@times,W',H0');

A11 = VH * V; 
A12 = VH * W + eye(M); 
A21 = WH * V + eye(M);
A22 = WH * W;

% more expensive, but numerically more benign:
A = [A11, A12; A21, A22];
r = -H0.*dfx + [VH', WH'] * (A \ ([VH; WH] * dfx));
