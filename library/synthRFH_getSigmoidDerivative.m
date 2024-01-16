function dQdV = synthRFH_getSigmoidDerivative( ...
  potential, maxrate, threshlevel, threshdeviation )

% function dQdV = synthRFH_getSigmoidDerivative( ...
%   potential, maxrate, threshlevel, threshdeviation )
%
% This gets the derivative of the Robinson 2002 activation function with
% respect to cell-body potential, for a given potential.
%
% "potential" is the cell body potential (V) to evaluate the derivative at.
%   This may be a vector or matrix, to evaluate several potentials.
% "maxrate" is the maximum firing rate (Q_max).
% "threshlevel" is the average neuron threshold (theta).
% "threshdeviation" is the standard deviation of the average neuron threshold
%   (sigma - not sigma prime!).
%
% "dQdV" is the derivative of the firing rate with respect to potential,
%   at the specified potential.


% The full derivation is in the model guide.
% Short version:
% Q(V) = Q_max L( (V - Vth) / sigmaprime )
% L'(V) is known to equal L(V) * ( 1 - L(V) )
% Applying the chain rule gets Q'(V) in terms of Q(V).


sigmaprime = threshdeviation * sqrt(3) / pi;

firingrate = synthRFH_getSigmoid( ...
  potential, maxrate, threshlevel, threshdeviation );

dQdV = firingrate .* ( 1 - firingrate ./ maxrate );
dQdV = dQdV / sigmaprime;


% Done.
end


%
% This is the end of the file.
