function firingrate = synthRFH_getSigmoid( ...
  potential, maxrate, threshlevel, threshdeviation )

% function firingrate = synthRFH_getSigmoid( ...
%   potential, maxrate, threshlevel, threshdeviation )
%
% This converts a cell-body potential (V in Robinson's model) to a firing
% rate (Q in Robinson's model).
%
% This is described in eq. 1 of Robinson 2002:
% https://journals.aps.org/pre/abstract/10.1103/PhysRevE.65.041924
% And in terms of sigma, not sigma-prime, in eq. m4 of Freyer 2011:
% https://www.jneurosci.org/content/31/17/6353.short
%
% "potential" is the cell body potential (V) to convert. This may be a
%   vector or matrix, to convert several potentials.
% "maxrate" is the maximum firing rate (Q_max).
% "threshlevel" is the average neuron threshold (theta).
% "threshdeviation" is the standard deviation of the average neuron threshold
%   (sigma - not sigma prime!).
%
% "firingrate" is the resulting firing rate (Q, or S(V)).


sigmaprime = threshdeviation * sqrt(3) / pi;

firingrate = ...
  maxrate ./ ( 1 + exp( - (potential - threshlevel) / sigmaprime ) );


% Done.
end


%
% This is the end of the file.
