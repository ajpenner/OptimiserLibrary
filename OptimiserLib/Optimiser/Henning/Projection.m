function k = Projection(x,E,X0,LS,hyp,U)
% returns the projection ks of the covariance between the integral observations
% along line searches and the observation locations x. x may be a matrix of size
% N x k. Otherwise, it is assumed to be a column vector of length N.
%
% If a numerical conditioner U is provided, the returned k is a normalized form
% k_{nm} = k0 _nm * U_m;
%
% Philipp Hennig, January 2012

%% unpack and pre-compute
sf2  = hyp.sf2; ell = hyp.ell; V = hyp.V;        % unpack kernel hyperparameters
[N,i]= size(E);

isLimitedMemory = min(size(V))==1 && ~all(size(V)==1);

En   = bsxfun(@rdivide,E,ell); X0n = bsxfun(@rdivide,X0,ell); 
xn   = bsxfun(@rdivide,x,ell);
EE   = En' * En; EX0  = En' * X0n; Ex = En' * xn; 
XX = xn' * xn; XX0 = xn' * X0n; X0X0 = X0n' * X0n;

MT   = 0; for j = 1:i; MT = MT + length(LS{j}.a); end
if isLimitedMemory
  VE = bsxfun(@times,V,E);
else
  VE   = V * E;
end

k  = zeros(N,MT,size(x,2));
h0 = 0;
for j = 1:i                                               % for each line search
    % pre-computations for speed
    a2 = EE(j,j);
    upj = LS{j}.a; dnj = LS{j}.b; %[0; LS{j}.a(1:end-1)];
    for h = 1 : length(LS{j}.a)                   % for each ls-location in ls j
        b = EX0(j,j) + dnj(h) * EE(j,j) - Ex(j,:);
        c = diag(XX)' - 2 * (XX0(:,j)' + dnj(h) * Ex(j,:)) ...
            + X0X0(j,j) + dnj(h) * dnj(h) * EE(j,j) + 2 * dnj(h) * EX0(j,j);
        kn = sf2 * exp(-(c-b.*b./a2)/2) * sqrt(pi/(2*a2)) .* ... % 
            ( erf(((upj(h) - dnj(h)) * a2 + b)/sqrt(2*a2)) - erf(b/sqrt(2*a2)) );
        k(:,h0 + h,:) = bsxfun(@times,VE(:,j),kn);
    end
    h0 = h0 + length(LS{j}.a);
end

if nargin > 5 % numerical conditioning
    k = bsxfun(@rdivide,k,U');
end
    
