function potential = synthRFH_getSigmoidInverse( ...
  firingrate, maxrate, threshlevel, threshdeviation )

% function potential = synthRFH_getSigmoidInverse( ...
%   firingrate, maxrate, threshlevel, threshdeviation )
%
% This converts a firing rate (Q in Robinson's model) back to a cell-body
% potential (V in Robinson's model), performing the inverse of
% synthRFH_getSigmoid().
%
% "firingrate" is the firing rate (Q, or S(V)). This may be a vector or
%   matrix, to convert several firing rates into potentials.
% "maxrate" is the maximum firing rate (Q_max).
% "threshlevel" is the average neuron threshold (theta).
% "threshdeviation" is the standard deviation of the average neuron threshold
%   (sigma - not sigma prime!).
%
% "potential" is the cell body potential (V).


sigmaprime = threshdeviation * sqrt(3) / pi;

scratch = maxrate ./ firingrate;
potential = threshlevel - sigmaprime * log(scratch - 1);


% Done.
end


%
% This is the end of the file.
