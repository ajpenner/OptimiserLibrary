function [sKs,ks,U,E,X0,LS,Y,S,t] = SlimGramMatrix(tolerance,MEMLIM,sKs,ks,U,E,X0,LS,Y,S,t)
% Philipp Hennig, February 2012

% check pre condition
ls = 0;
for l=1:length(LS)
  ls = ls + length(LS{l}.a);
end
assert(size(sKs,1) == ls);

M = size(sKs,1);
if M > 2 % do nothing if there are only a handful of data around
    % similarity of colums (rows) in the Gram matrix:
    Diff = max(abs(bsxfun(@minus,reshape(sKs,[M,1,M]),reshape(sKs,[1,M,M]))),[],3);
    Diffn = max(abs(bsxfun(@plus,reshape(sKs,[M,1,M]),reshape(sKs,[1,M,M]))),[],3);
    l = 1;                                                % index of linesearch
    e = 1;                                            % index within linesearch
    m = 1;
    while m <= size(sKs,1) % Matlab doesn't like fiddling with for-loop variables
        if e > length(LS{l}.a); l = l+1; e = 1; end;              % bookkeeping
        for i = m+1:size(sKs,1)
            if Diff(m,i) < tolerance || Diffn(m,i) < tolerance ... % if this evaluation resembles another one
                    || size(sKs,1) > MEMLIM
                fprintf 'X'
                sKs(:,m) = []; sKs(m,:) = []; % remove this entry from sKs ...
                ks(:,m) = [];                 % ... from ks ...
                LS{l}.a(e) = [];              % ... from LS ...
                LS{l}.b(e) = [];
                LS{l}.F(e) = [];
                LS{l}.dF(:,e)= [];
                Y(:,m) = [];                  % ... from Y ...
                S(:,m) = [];                  % ... from S ...
                U(m) = [];                    % ... and from U.
                Diff(m,:) = []; Diff(:,m) = [];        % also from our own Diff
                Diffn(m,:) = []; Diffn(:,m) = [];        % also from our own Diff
                e = e-1;                                          % bookkeeping
                if isempty(LS{l}.a) % if an entire linesearch is gone ...
                    % ... then remove that linesearch from LS ...
                    LS(l) = [];
                    E(:,l)  = [];  % ... from E ...
                    X0(:,l) = [];   % ... and from X0.
                    t = t - 1;
                    e = 0;                                        % bookkeeping
                end
                m = m-1; break                            % this column is gone
            end
        end
        m = m + 1; e = e + 1;                                     % bookkeeping
    end
end