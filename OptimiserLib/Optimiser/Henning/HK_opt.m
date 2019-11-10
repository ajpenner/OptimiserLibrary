function [X,fX,i,outs] = HK_opt(X,F,opts,varargin)
%% minimize function F, using a Bayesian Quasi-Newton method. 
% The outward interface of this code is adapted from Carl Rasmussen's
% minimize.m :
%
% Usage: [X,fX,i,outs] = HK_opt(X,F,opts,varargin)
% if you do not want to use all the bells and whistles, but just want
% something that works reasonably well, the minimial syntax is
% [X,fX] = HK_opt(X,F,l,varargin), with a scalar l: the maximum number of
% line searches / function evaluations (see below).
% where
%   X    is an initial guess (any type: vector, matrix, cell array, struct)
%   F    is the objective function (function pointer or name)
%   opts are further parameters 
%         if opts is a number, it corresponds to opts.length below
%         opts.length       allowed # linesearches; if -ve: minus # func evals
%         opts.noise        if <=0 scalar: evaluations are without noise. 
%                           if >0 scalar : standard deviation (not variance!) 
%                                          of noise on each element of gradient
%         opts.model        model used for the minimization:
%         'constant_Gauss','DFP','kernel_Gauss', 'nonparametric_Gauss', 'L_NPQN'
%   other, ...     other parameters, passed to the function F
%   X     returned minimizer
%   fX    vector of function values showing minimization progress
%   i     final number of linesearches or function evaluations
%   outs  diagnostic outputs for research purposes
% The function F must take the following syntax [f, df] = F(X, other, ...)
% where f is the function value and df its partial derivatives. The types of X
% and df must be identical (vector, matrix, cell array, struct, etc).
%
% Copyright (C) 2012 by Philipp Hennig & Martin Kiefel, 2012-01-17

%% Copyright Notice:
% Permission is hereby granted, free of charge, to any person OBTAINING A COPY
% OF THIS SOFTWARE AND ASSOCIATED DOCUMENTATION FILES (THE "SOFTWARE"), TO DEAL
% IN THE SOFTWARE WITHOUT RESTRICTION, INCLUDING WITHOUT LIMITATION THE RIGHTS
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

%% top-level code. Much of this stems from minimize.m. See copyright notice below.
if isnumeric(opts), opts = struct('length', opts); end       % convert to struct
if opts.length > 0, 
    opts.S = 'linesearch #'; else opts.S = 'function evaluation #'; end;
x = unwrap(X);                                % convert initial guess to vector
if ~isfield(opts,'method'), 
    %if length(x) > 1000, opts.method = @LBFGS; else opts.method = @BFGS; end; end   
    opts.method = @L_NPQN; % set default method
end
if ~isfield(opts,'verbosity'), opts.verbosity = 2; end   % default 1 line text output
if ~isfield(opts,'MFEPLS'), opts.MFEPLS = 20; end    % Max Func Evals Per Line Search
if ~isfield(opts,'MSR'), opts.MSR = 100; end                % Max Slope Ratio default
f(F, X, varargin{:});                                   % set up the function f
[fx, dfx] = f(x);                      % initial function value and derivatives 
if opts.verbosity, printf('Initial Function Value %4.6e\r', fx); end
if opts.verbosity > 2,
  clf; subplot(211); hold on; xlabel(opts.S); ylabel('function value');
  plot(opts.length < 0, fx, '+'); drawnow;
end
[x, fX, i, outs] = feval(opts.method, x, fx, dfx, opts);  % minimize using direction method 
X = rewrap(X, x);                   % convert answer to original representation
if opts.verbosity, printf('\n'); end

function [x,fx,i,outs] = DFP(x0,fx0,dfx0,opts)
% classic DFP method
if ~isfield(opts,'noise'), opts.noise = 0; end
if ~isfield(opts,'SIG'), opts.SIG = 0.5; end % default for line search qualityif ~isfield(opts,'V0'), opts.V0 = 1 * eye(N); end
if ~isfield(opts,'B0'); opts.B0 = 1 * eye(length(x0)); end
i = opts.length < 0; % initialize resource counter
% run variables
x = x0; fx = fx0; r = -dfx0; u = -r'*r; c = -1/u; ok = 0; % linear search params
B = opts.B0; % prior mean
outs.B{1} = B; outs.X = x; outs.AllEvals.x = []; outs.AllEvals.f = [];
t = 1;
while i < abs(opts.length)
    [x,c,fx0,dfx,i, ~, ~, ~, allEvals] = lineSearch(x0,fx0,dfx0,r,u,c,i,opts);
    outs.AllEvals.x = [outs.AllEvals.x, allEvals.x];
    outs.AllEvals.f = [outs.AllEvals.f, allEvals.f];
    if i < 0
        i = -i; if ok, ok = 0; else break; end;   % try steepest descent or stop
    else
        if all(dfx == 0); ok = 1; break; end; % we are at the local minimum.
        ok = 1; s = x - x0; y = dfx - dfx0;
        outs.y{t} = y; outs.s{t} = s;
        % Bayesian regression update:
        YBS = y - B * s; VS = y; SYBS = s' * YBS;
        SVSB = s' * VS + opts.noise * opts.noise; 
        VSSV = (VS * VS') / SVSB;
        T1  = YBS * VS'; T1 = (T1 + T1') / SVSB;
        T2  = VSSV * (SYBS / SVSB);
        B   = B + T1 - T2;
    end
    r = - B \ dfx;  % Newton direction. Analytic inverse to be implemented soon.
    u = r' * dfx; x0 = x; dfx0= dfx; fx = [fx;fx0];  %#ok<AGROW>
    t = t + 1;
    outs.B{t} = B;
    outs.X = [outs.X, x];
end
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
        YBS = y - B * s; 
		VS = V * s; 
		SYBS = s' * YBS; % VS = V * s
        % check noise required for positive definiteness
        %a = y' * (B \ VS); d = (y' * (B \ YBS)) * (VS' * (B \ VS));
        %b1 = -a/2 + sqrt(d-3*a); b2 = -a/2 - sqrt(d-3*a);
        %if ~isreal(b1); b1 = 0; b2 = 0; end
        %bt = max([0,b1,b2,b]) + 0.5;
        SVSB = s' * VS + b; 
        VSSV = (VS * VS') / SVSB;
        T1  = YBS * VS'; 
		T1 = (T1 + T1') / SVSB;
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

function [x,fx,i,outs] = nonparametric_Gauss(x0,fx0,dfx0,opts) %#ok<*DEFNU>
% Bayesian quasi-Newton method with nonparametric prior for the Hessian
if ~isfield(opts,'SIG'), opts.SIG = 0.5; end
if ~isfield(opts,'noise'), opts.noise = 0; end
% set up prior
N = length(x0);
if ~isfield(opts,'B0') && ~isfield(opts,'H0'); opts.B0 = 1 * eye(N); opts.H0 = 1 * eye(N);
elseif (~isfield(opts,'B0') && isfield(opts,'H0')); opts.B0 = inv(opts.H0);
elseif (isfield(opts,'B0') && ~isfield(opts,'H0')); opts.H0 = inv(opts.B0);
end
%if ~isfield(opts,'H0'); opts.B0 = 1 * eye(N); end
if ~isfield(opts,'V0'); opts.V0 = 1 * eye(N); end
if ~isfield(opts,'hyp'); opts.hyp = log(ones(N+1,1)); end % hyperparameters to be learned in the future
if ~isfield(opts,'H'); opts.H = -3; end % number of linesearches to keep around
hyp.ell = exp(opts.hyp(1:N)); hyp.sf2 = exp(2*opts.hyp(end)); hyp.V = opts.V0;
H = opts.H;
i = opts.length < 0; % initialize resource counter
% run variables
x = x0; fx = fx0; r = -dfx0; u = -r'*r; c = -1/u; ok = 0; % linear search params
E     = r / norm(r); X0 = x; % unit line search storage
B = opts.B0; b = max(opts.noise * opts.noise,1.0e-8);
B0= opts.B0; H0 = opts.H0;
outs.B{1} = B; outs.X = x; outs.V{1} = hyp.V; t = 1;
outs.AllEvals.x = []; outs.AllEvals.f = [];
X0 = x; S = []; Y = []; % storage for past observations 
sKs = []; U = []; NewtonFlag = 0;
while i < abs(opts.length)
    [x,c,fx0,dfx,i,A,F,dF,allEvals] = lineSearch(x0,fx0,dfx0,r,u,c,i,opts);
    outs.AllEvals.x = [outs.AllEvals.x, allEvals.x];
    outs.AllEvals.f = [outs.AllEvals.f, allEvals.f];
    %A = A(end); F = F(end); dF = F(:,end);     % test: use only last observation
    if i < 0 
        i = -i; if ok, ok = 0; else break; end;   % try steepest descent or stop
    else
        if H >= 0
            LS{min(t,H+1)}.a = A; LS{min(t,H+1)}.F = F; LS{min(t,H+1)}.dF = dF;
            LS{min(t,H+1)}.b = [0; A(1:end-1)];
            %     elseif length(A) > 1
            %         x0 = x0 + A(end-1) * E(:,end); X0(:,end) = x0;
            %         LS{t}.a = A(end) - A(end-1);
            %         LS{t}.b = 0;
            %         LS{t}.F = F(end); LS{t}.dF = dF(:,end); dfx0 = dF(:,end-1);
        else
            LS{t}.a = A; LS{t}.b = [0;A(1:end-1)]; LS{t}.F = F; LS{t}.dF = dF;
            %if t > 1; LLS0{tt0}.a = A; LLS0{tt0}.b = [0;A(1:end-1)]; end
        end;
        if all(dfx == 0); ok = 1; break; end; % we are at the local minimum.
        ok  = 1;   
        Delta = LS{t}.a - LS{t}.b; % A - [0;A(1:end-1)];
        s = E(:,end) * Delta'; y = LS{t}.dF - [dfx0, LS{t}.dF(:,1:end-1)]; %y = dfx - dfx0;
        Y   = [Y,y]; S = [S,s]; %#ok<AGROW>
        %gamma = 1; % mean(diag(y' * y ./ (y'*s)))
        %B0 = gamma * eye(N); H0 = eye(N) / gamma;
        [sKs,U] = SchurUpdate(sKs .* (U*U'),E,X0,LS,hyp);% update Schur complement of covariance
        %if t > 1; [sks1,U1] = SchurUpdate(ssKs0 .* (UU0*UU0'),EE0,XX0,LLS0,hyp); end
        ks  = Projection(x,E,X0,LS,hyp,U);  % update kernel projection
        %EE0 = E; XX0 = X0; LLS0 = LS; ssKs0 = sKs; UU0 = U; tt0 = t;
        [sKs,ks,U,E,X0,LS,Y,S,t] = SlimGramMatrix(1e-2,100,sKs,ks,U,E,X0,LS,Y,S,t); % ensure numerical stability
        try
            L   = chol(sKs + b * eye(size(sKs,1)));
        catch
            fprintf 'reached numerical criticality. Switching to BFGS'
            p = opts; p.H0 = inv(B); p.length = p.length - sign(p.length) * i;
            outs.bfgsswitch = t;
            [xbfgs, fxbfgs, ibfgs, outs_bfgs] = BFGS(x, fx0, dfx, p);
            fx = [fx; fxbfgs];
            x = xbfgs;
            return;
        end
        YBS = bsxfun(@rdivide,Y - B0 * S,U'); % numerical conditioning
        SYBS = bsxfun(@ldivide,U,S' * YBS);   % numerical conditioning
        ksL = ks / L;
        T1  = YBS / L * ksL'; T1 = T1 + T1';
        T2  = ksL * (L' \ SYBS) / L * ksL';
        B   = B0 + T1 - T2;
        B   = 0.5 * (B + B'); % symmetrize, because the method does not explicitly encode this.
        t = t + 1;
%         if size(sKs,1) > 20; % only if the posterior is beginning to settle
%             % adjust hyperparameters:
%             v = ones(N,1); % power iteration
%             for i = 1:10;
%                 q = v / norm(v);
%                 v = abs(B) * q;
%                 lambda = q' * v;
%             end
%             hyp.V = 0.9 * hyp.V + 0.1 * diag(lambda * q);
%             %keyboard;
%         end
    end
    if ok; 
        Hnew = InverseHessian(ksL,L,YBS,SYBS,H0);
        %keyboard;
        % r = - B \ dfx;
        r = - Hnew * dfx;  % Newton direction. Analytic inverse to be implemented soon.
    else r = - dfx; % line search failed. Try gradient-descent.
    end  
    if r' * dfx > 0; r = -r; NewtonFlag = [NewtonFlag; 0]; else NewtonFlag = [NewtonFlag; 1]; end % heuristic way of avoiding indefinite updates
    %r = r + 0.001 * randn(N,1);
    % keyboard;
    u = r' * dfx; x0 = x; dfx0= dfx; fx = [fx;fx0];  %#ok<AGROW>
    outs.NewtonFlag = NewtonFlag;
    outs.Y = Y; outs.S  = S;
    if ok
        outs.E = E; outs.X0 = X0; outs.LS = LS; outs.U = U; outs.sKs = sKs;
        E = [E, r / norm(r)]; X0 = [X0, x];
    end
    %EE0 = [EE0, r / norm(r)]; XX0 = [XX0, x]; tt0 = tt0 + 1;
    outs.B{t} = B; outs.X = [outs.X, x]; outs.V{t} = hyp.V;
    % removing old observations from the pipeline:
    if H >= 0 && t > H + 1 % incorporate old observations into mean
        B0 = outs.B{t-H}; % move mean to posterior mean from H observations ago
        E  = E (:,end-H:end); % remove old linesearches from storage
        X0 = X0(:,end-H:end);
        if H > 0; 
            ls = cell(H,1); for h = 1:H; ls{h} = LS{end-H+h}; end; LS = ls;
            mt = 0; for h = 1:H; mt = mt + length(LS{h}.a); end;
            sKs = sKs(end-mt+1:end,end-mt+1:end); 
            S   = S(:,end-mt+1:end);
            Y   = Y(:,end-mt+1:end);
            U   = U(end-mt+1:end);
        else LS = {}; sKs = []; S = []; Y = []; U = []; end;
    end
end

function [x,fx,i,outs] = L_NPQN(x0, fx0, dfx0, opts)
% clean, fast implementation of limited-memory nonparametric quasi-Newton method
% this is essentially the same as nonparametric_Gauss above, but assuming a
% diagonal prior mean. The code has also been cleaned up a bit.
%
% Philipp Hennig, February 2012
% Bayesian quasi-Newton method with nonparametric prior for the Hessian

%% setup
N = length(x0);

if ~isfield(opts,'SIG'), opts.SIG = 0.5; end
if ~isfield(opts,'noise'), opts.noise = 0; end
if ~isfield(opts,'mem'), opts.mem = min(25, N); end

% prior mean:
if ~isfield(opts,'B0') && ~isfield(opts,'H0'); opts.B0 = ones(N,1); opts.H0 = ones(N,1);
elseif (~isfield(opts,'B0') && isfield(opts,'H0')); opts.B0 = 1./opts.H0;
elseif (isfield(opts,'B0') && ~isfield(opts,'H0')); opts.H0 = 1./opts.B0;
end
assert (isequal(size(opts.B0),[N,1])); assert (isequal(size(opts.H0),[N,1]));

% prior covariance:
if ~isfield(opts,'V0'); opts.V0 = ones(N,1); end  % output covariances
if ~isfield(opts,'hyp'); opts.hyp = log(ones(N+1,1)); end % input covariances
assert (isequal(size(opts.V0),[N,1]));

% run variables
i = opts.length < 0; % initialize resource counter
x = x0; fx = fx0; r = -dfx0; u = -r'*r; c = -1/u; ok = 0; % linear search params
E     = r ./ norm(r); X0 = x; % unit line search storage
B0    = opts.B0; H0 = opts.H0; % prior mean (and its inverse), as diagonal mats.
b = max(opts.noise * opts.noise,1.0e-8); % noise
hyp.ell = exp(opts.hyp(1:N)); hyp.sf2 = exp(2*opts.hyp(end)); hyp.V = opts.V0;
outs.X = x; t = 1;
X0 = x; S = []; Y = []; % storage for past observations 
K = []; U = []; NewtonFlag = 0;
while i < abs(opts.length)
    [x,c,fx0,dfx,i,A,F,dF] = lineSearch(x0,fx0,dfx0,r,u,c,i,opts);
    if i < 0 % line search failed
        i = -i; if ok, ok = false; else break; end; % try steepest descent, or stop.
    else
        LS{t}.a = A; LS{t}.b = [0;A(1:end-1)]; LS{t}.F = F; LS{t}.dF = dF; %store
        if all(abs(dfx) < 1e-10); ok = 1; break; end; % successully found minimimum.
        ok = true;  % this line search went fine.
        S = [ S, E(:,end) * (LS{t}.a - LS{t}.b)' ]; % update locations
        Y = [ Y, LS{t}.dF - [dfx0, LS{t}.dF(:,1:end-1)] ]; % update observations
        [K,U] = SchurUpdate(K .* (U * U'), E, X0, LS, hyp); % Gram matrix
        k     = Projection(x, E, X0, LS, hyp, U); % kernel projection
        % limit computation costs, ensure numerical stability:
        [K,k,U,E,X0,LS,Y,S,t] = SlimGramMatrix(1e-2,opts.mem,K,k,U,E,X0,LS,Y,S,t);
        try
            L = chol(K + b * eye(size(K)));
        catch  %#ok<CTCH>
            fprintf 'reached numerical criticality. Switching to L-BFGS'
            p = opts; p.length = p.length - sign(p.length) * i;
            [xbfgs,fxbfgs,~,outs_bfgs] = LBFGS(x,fx0,dfx,p);
            keyboard;
            x = xbfgs; outs.X = [outs.X, outs_bfgs.X]; fx = [fx; fxbfgs];
            return;
        end
        YBS  = bsxfun(@rdivide,Y - bsxfun(@times,B0,S), U'); % numerical conditioning
        SYBS = bsxfun(@ldivide,U, S' * YBS);    % numerical conditioning
        kL   = k / L;
        t    = t + 1;
    end
    if ok; r  = L_UpdateDirection(kL,L,YBS,SYBS,opts.H0,dfx); else r = -dfx; end;
    if r' * dfx > 0; r = -r; NewtonFlag = [NewtonFlag, 0]; ...
    else NewtonFlag = [NewtonFlag, 1]; end; % heuristic avoidance of nonposdef update.
    u = r' * dfx; x0 = x; dfx0 = dfx; fx = [fx;fx0];
    if ok
        E = [E, r/norm(r)]; X0 = [X0, x];
    end 
end
outs.NewtonFlag = NewtonFlag; outs.Y = Y; outs.S = S; outs.X = X0;
outs.X0 = X0; outs.E = E; outs.LS = LS; outs.U = U; outs.K = K;

%% Carl Rasmussen's linesearch and utilities. For all code below this line, 
% the following  copyright notice applies:
% Copyright (C) 1996 - 2011 by Carl Edward Rasmussen, 2011-10-06.
% Permission is hereby granted, free of charge, to any person OBTAINING A COPY
% OF THIS SOFTWARE AND ASSOCIATED DOCUMENTATION FILES (THE "SOFTWARE"), TO DEAL
% IN THE SOFTWARE WITHOUT RESTRICTION, INCLUDING WITHOUT LIMITATION THE RIGHTS
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
function [x, fx, i, outs] = BFGS(x0, fx0, dfx0, p)
if ~isfield(p, 'SIG'), p.SIG = 0.5; end       % default for line search quality
i = p.length < 0;                                 % initialize resource counter
if isfield(p, 'H0'); H = p.H0; else H = eye(length(x0)); end
x = x0; fx = fx0; r = -dfx0; s = -r'*r; b = -1/s; ok = 0; 
outs.X = x; outs.AllEvals.x = []; outs.AllEvals.f = [];
while i < abs(p.length)
  [x, b, fx0, dfx, i, ~, ~, ~, allEvals] = lineSearch(x0, fx0, dfx0, r, s, b, i, p); 
  outs.AllEvals.x = [outs.AllEvals.x, allEvals.x];
  outs.AllEvals.f = [outs.AllEvals.f, allEvals.f];
  if i < 0
    i = -i; if ok, ok = 0; else break; end;              % try steepest or stop
  else
    if all(dfx == 0); ok = 1; break; end;
    ok = 1; t = x - x0; y = dfx - dfx0; ty = t'*y; Hy = H*y;
    H = H + (ty+y'*Hy)/ty^2*t*t' - 1/ty*Hy*t' - 1/ty*t*Hy';       % BFGS update
  end
   r = -H*dfx; 
   if ~ok; r = -dfx; end; % unstable, gradient descent
   % keyboard;
   s = r'*dfx; x0 = x; dfx0 = dfx; fx = [fx; fx0];
   outs.X = [outs.X,x];
end

function [x, fx, i, outs] = LBFGS(x0, fx0, dfx0, p);
if ~isfield(p, 'SIG'), p.SIG = 0.5; end       % default for line search quality
n = length(x0); k = 0; ok = 0; x = x0; fx = fx0; bs = -1/p.MSR;
if isfield(p, 'mem'), m = p.mem; else m = min(100, n); end    % set memory size
a = zeros(1, m); t = zeros(n, m); y = zeros(n, m);            % allocate memory
i = p.length < 0;                                 % initialize resource counter
outs.X = x;
outs.AllEvals.x = []; outs.AllEvals.f = [];
while i < abs(p.length)
  q = dfx0;
  for j = rem(k-1:-1:max(0,k-m),m)+1
    a(j) = t(:,j)'*q/rho(j); q = q-a(j)*y(:,j);
  end
  if k == 0, r = -q/(q'*q); else r = -t(:,j)'*y(:,j)/(y(:,j)'*y(:,j))*q; end
  for j = rem(max(0,k-m):k-1,m)+1
    r = r-t(:,j)*(a(j)+y(:,j)'*r/rho(j));
  end
  s = r'*dfx0; if s >= 0, r = -dfx0; s = r'*dfx0; k = 0; ok = 0; end
  b = bs/min(bs,s/p.MSR);              % suitable initial step size (usually 1)
  if isnan(r) | isinf(r)                                % if nonsense direction
    i = -i;                                              % try steepest or stop
  else
    [x, b, fx0, dfx, i, ~, ~, ~, allEvals] = lineSearch(x0, fx0, dfx0, r, s, b, i, p); 
    outs.AllEvals.x = [outs.AllEvals.x, allEvals.x]; 
    outs.AllEvals.f = [outs.AllEvals.f, allEvals.f];
  end
  if i < 0                                              % if line search failed
    i = -i; if ok, ok = 0; k = 0; else break; end        % try steepest or stop
  else
    j = rem(k,m)+1; t(:,j) = x-x0; y(:,j) = dfx-dfx0; rho(j) = t(:,j)'*y(:,j);
    ok = 1; k = k+1; bs = b*s;
  end
  x0 = x; dfx0 = dfx; fx = [fx; fx0];                  % replace and add values
  outs.X = [outs.X, x];
end

function [x, a, fx, df, i, alpha, fX, dF, allEvals] = lineSearch(x0, f0, df0, d, s, a, i, p)
allEvals.f = []; 
allEvals.x = [];
alpha = []; 
fX = []; 
dF = [];           % store function evaluations for inference
INT = 0.1; 
EXT = 5.0;                    % interpolate and extrapolation limits
	if p.length < 0
		LIMIT = min(p.MFEPLS, -i-p.length); 
	else 
		LIMIT = p.MFEPLS; 
	end
p0.x = 0.0; 
p0.f = f0; 
p0.df = df0; 
p0.s = s; 
p1 = p0;         % init p0 and p1
j = 0; 
p3.x = a; 
wp(p0, p.SIG, p.SIG/2);   % set step & Wolfe-Powell conditions
	if p.verbosity > 2
  		A = [-a a]/5; nd = norm(d);
  		subplot(212); 
  		hold off; 
  		plot(0, f0, 'k+'); 
  		hold on; 
  		plot(nd*A, f0+s*A, 'k-');
  		xlabel('distance in line search direction'); 
  		ylabel('function value');
	end

	while 1                               % keep extrapolating as long as necessary
  		ok = 0; 
  		while ~ok && j < LIMIT
    		j = j+1;
    		try           % try, catch and bisect to safeguard extrapolation evaluation
      		[p3.f p3.df] = f(x0+p3.x*d); 
	  		p3.s = p3.df'*d; 
	  		ok = 1; 
	  		allEvals.f = [allEvals.f, p3.f];
	  		allEvals.x = [allEvals.x, x0 + p3.x * d];
      		if p3.f < f0
          		alpha = [alpha; norm(p3.x*d)];
		  		fX = [fX, p3.f]; 
		  		dF = [dF, p3.df]; % update storage
      		end
      		if isnan(p3.f+p3.s) || isinf(p3.f+p3.s)
		        error('Objective function returned Inf or NaN','');
      		end;
    		catch
	  			if p.verbosity > 1, printf('\n'); 
	  				warning(lasterr); 
	  			end % warn or silence
	  			p3.x = (p1.x+p3.x)/2; 
	  			ok = 0; 
	  			p3.f = NaN;             % bisect, and retry
    		end
  		end
  		if p.verbosity > 2
		    plot(nd*p3.x, p3.f, 'b+'); 
			plot(nd*(p3.x+A), p3.f+p3.s*A, 'b-'); 
			drawnow
  		end
  		if wp(p3) || j >= LIMIT
      		break; 
  		end                                    % done?
  		p0 = p1; 
 		p1 = p3;                                 % move points back one unit
  		a = p1.x-p0.x; 
  		b = minCubic(a, p1.f-p0.f, p0.s, p1.s);  % cubic extrapolation
  		if ~isreal(b) || isnan(b) || isinf(b) || b < a || b > a*EXT     % if b is bad
    		b = a*EXT;                                % then extrapolate maximum amount
  		end
  		p3.x = p0.x+max(b, p1.x+INT*a);             % move away from current, off-set 
	end

	while 1                               % keep interpolating as long as necessary
	  if p1.f > p3.f
	  	p2 = p3; 
	  else 
	  	p2 = p1; 
	  end          % make p2 the best so far
  	  if wp(p2) > 1 || j >= LIMIT
	  	break; 
	  end                                % done?
	  a = p3.x-p1.x; 
	  b = minCubic(a, p3.f-p1.f, p1.s, p3.s);  % cubic interpolation
  	  if ~isreal(b) || isnan(b) || isinf(b) || b < 0 || b > a         % if b is bad
    		b = a/2;                                                      % then bisect
  	  end
  	  p2.x = p1.x+min(max(b, INT*a), (1-INT)*a);  % move away from current, off-set 
 	  [p2.f p2.df] = f(x0+p2.x*d); 
	  j = j+1; 
	  p2.s = p2.df'*d;
  	  allEvals.f = [allEvals.f, p2.f]; allEvals.x = [allEvals.x, x0 + p2.x * d];
  	  if p2.f < f0
      	alpha = [alpha; norm(p2.x*d)]; 
		fX = [fX, p2.f]; 
		dF = [dF, p2.df];  % update storage
  	  end
  	  if p.verbosity > 2
    	plot(nd*p2.x, p2.f, 'r+'); plot(nd*(p2.x+A), p2.f+p2.s*A, 'r'); drawnow
  	  end
  	  if wp(p2) > -1 && p2.s > 0 || wp(p2) < -1
	  	p3 = p2; 
	  else 
	  	p1 = p2; 
	  end % bracket
	end

	x = x0+p2.x*d; 
	fx = p2.f; 
	df = p2.df; 
	a = p2.x;        % return the value found
	if p.length < 0
		i = i+j; 
	else 
		i = i+1; 
	end % count func evals or line searches
	if p.verbosity
		printf('\r%s %6i;  value %4.6e\r', p.S, i, fx); 
	end 
	if wp(p2) < 2
		i = -i; 
	end                                   % indicate faliure 
	
	if p.verbosity > 2
		if i>0
			plot(norm(d)*p2.x, fx, 'go'); 
		end
  		subplot(211); 
		plot(abs(i), fx, '+'); 
		drawnow;
	end
	
function x = minCubic(x, df, s0, s1)  % return minimizer of approximating cubic
A = -6*df+3*(s0+s1)*x; 
B = 3*df-(2*s0+s1)*x; 
x = -s0*x*x/(B+sqrt(B*B-A*s0*x));

function y = wp(p, SIG, RHO)
persistent a b c sig rho;
if nargin == 3    % if three arguments, then set up the Wolfe-Powell conditions
  a = RHO*p.s; 
  b = p.f; 
  c = -SIG*p.s; 
  sig = SIG; 
  rho = RHO; 
  y= 0;
else
  if p.f > a*p.x+b                                  % function value too large?
    if a > 0
		y = -1; 
	else 
		y = -2; 
	end                  
  else
    if p.s < -c, y = 0; elseif p.s > c, y = 1; else y = 2; end
    if sig*abs(p.s) > c, a = rho*p.s; b = p.f-a*p.x; c = sig*abs(p.s); end
  end
end

function [fx, dfx] = f(varargin)
persistent F p;
if nargout == 0
  p = varargin; 
  if ischar(p{1})
  	F = str2func(p{1}); 
  else 
  	F = p{1}; 
  end
else
  [fx, dfx] = F(rewrap(p{2}, varargin{1}), p{3:end}); 
  dfx = unwrap(dfx);
end

function v = unwrap(s)   % extract num elements of s (any type) into v (vector) 
v = [];   
if isnumeric(s)
  v = s(:);                        % numeric values are recast to column vector
elseif isstruct(s)
  v = unwrap(struct2cell(orderfields(s))); % alphabetize, conv to cell, recurse
elseif iscell(s)                                      % cell array elements are
  for i = 1:numel(s)
  	v = [v; unwrap(s{i})]; 
  end % handled sequentially
end                                                   % other types are ignored

function [s v] = rewrap(s, v)    % map elements of v (vector) onto s (any type)
if isnumeric(s)
  if numel(v) < numel(s)
    error('The vector for conversion contains too few elements')
  end
  s = reshape(v(1:numel(s)), size(s));            % numeric values are reshaped
  v = v(numel(s)+1:end);                        % remaining arguments passed on
elseif isstruct(s) 
  [s p] = orderfields(s); p(p) = 1:numel(p);      % alphabetize, store ordering
  [t v] = rewrap(struct2cell(s), v);                 % convert to cell, recurse
  s = orderfields(cell2struct(t,fieldnames(s),1),p);  % conv to struct, reorder
elseif iscell(s)
  for i = 1:numel(s)             % cell array elements are handled sequentially
    [s{i} v] = rewrap(s{i}, v);
  end
end                                             % other types are not processed

function printf(varargin)
fprintf(varargin{:}); 
if exist('fflush','builtin')
	fflush(stdout); 
end
