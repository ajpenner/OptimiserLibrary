function H = InverseHessian(k_cholK,cholK,YBS,SYBS,H0)
% analytic inverse of the Hessian, constructed in O(N^2 M^2) time. 
%
% Philipp Hennig, February 2012

[~,M] = size(k_cholK);

SYBS = (SYBS + SYBS');

V = k_cholK / cholK';
W = YBS - 0.25 .* V * SYBS;

VH = V' * H0;
WH = W' * H0;

A11 = VH * V; 
A12 = VH * W + eye(M); 
A21 = WH * V + eye(M);
A22 = WH * W;

% more expensive, but numerically more benign:
A = [A11, A12; A21, A22];
H = H0 - [VH', WH'] * (A \ [VH; WH]);

% 8 times cheaper, but badly conditioned
%F1 = A11 - A12 * (A22 \ A21);
%F2 = A22 - A21 * (A11 \ A12);

%T1 = F1 \ VH - A11 \ (A12 /  F2) * WH;
%T2 = F2 \ (WH - (A21 / A11) * VH);

%H  = H0 - [VH', WH'] * [T1; T2];
% keyboard;