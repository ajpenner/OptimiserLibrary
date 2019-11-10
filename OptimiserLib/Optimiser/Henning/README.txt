minimize function F, using a Bayesian Quasi-Newton method.
The outward interface of this code is adapted from Carl Rasmussen's
minimize.m :

Usage: [X,fX,i,outs] = HK_opt(X,F,opts,varargin)
if you do not want to use all the bells and whistles, but just want
something that works reasonably well, the minimial syntax is
[X,fX] = HK_opt(X,F,l,varargin), with a scalar l: the maximum number of
line searches / function evaluations (see below).
where
  X    is an initial guess (any type: vector, matrix, cell array, struct)
  F    is the objective function (function pointer or name)
  opts are further parameters
        if opts is a number, it corresponds to opts.length below
        opts.length       allowed # linesearches; if -ve: minus # func evals
        opts.noise        if <=0 scalar: evaluations are without noise.
                          if >0 scalar : standard deviation (not variance!)
                                         of noise on each element of gradient
        opts.model        model used for the minimization:
        'constant_Gauss','DFP','kernel_Gauss', 'nonparametric_Gauss', 'L_NPQN'
  other, ...     other parameters, passed to the function F
  X     returned minimizer
  fX    vector of function values showing minimization progress
  i     final number of linesearches or function evaluations
  outs  diagnostic outputs for research purposes
The function F must take the following syntax [f, df] = F(X, other, ...)
where f is the function value and df its partial derivatives. The types of X
and df must be identical (vector, matrix, cell array, struct, etc).

Version 1.0
Copyright (C) 2012 by Philipp Hennig & Martin Kiefel, 2012-01-17
