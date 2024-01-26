function [ firingrates potentials ] = ...
  synthRFH_estimateOperatingPointExponential( ...
    modelparams, intcouplings, startpotentials )

% function [ firingrates potentials ] = ...
%   synthRFH_estimateOperatingPointExponential( ...
%     modelparams, intcouplings, startpotentials )
%
% This attempts to estimate the DC operating point of a Robinson neural model.
%
% Per the model guide, operating points with firing rates much less than the
% maximum are solutions to the equation:
%
% potentials = intcouplings * Q_0 * exp( potentials / sigmaprime )
%
% This function does a brute-force gradient descent search for operating
% points using "fsolve". This only finds one point; several points may
% exist.
%
% NOTE - Operating point firing rates _must_ be examined to confirm that
% they are much less than modelparams.qmax. If they are not several times
% smaller, the operating point is not correct.
%
% "modelparams" is a model parameter structure with the fields described in
%   MODELPARAMS.txt.
% "intcouplings" is a 4x4 matrix indexed by (destination,source) that
%   provides the coupling weights (in mV*s) between excitatory, inhibitory,
%   specific nucleus, and reticular nucleus neural populations.
% "startpotentials" is a vector containing cell potentials for the excitatory,
%   inhibitory, specific nucleus, and reticular nucleus populations used as
%   a starting point for further optimization. Set this to [] to call
%   synthRFH_estimateOperatingPointLinear() to generate starting potentials.
%
% "firingrates" is a vector containing firing rates for the excitatory,
%   inhibitory, specific nucleus, and reticular nucleus populations.
% "potentials" is a vector containing cell potentials for the excitatory,
%   inhibitory, specific nucleus, and reticular nucleus populations.


% Get a starting point if we weren't given one.

if isempty(startpotentials)
  [ scratch startpotentials ] = ...
    synthRFH_estimateOperatingPointLinear( modelparams, intcouplings );
end


% Set up the function handle for fsolve().

sigmaprime = modelparams.threshsigma * sqrt(3) / pi;
qnought = modelparams.qmax * exp( - modelparams.threshlevel / sigmaprime );

anonfunc = @(testvec) helper_calcSolutionError( ...
  qnought, sigmaprime, intcouplings, testvec );


% Call fsolve() and see what it gives us.

options = optimoptions('fsolve', 'Display', 'off');

potentials = fsolve(anonfunc, startpotentials, options);

firingrates = synthRFH_getSigmoid( potentials, ...
  modelparams.qmax, modelparams.threshlevel, modelparams.threshsigma );


% Done.
end


%
% Helper Functions

function solutionerrors = helper_calcSolutionError( ...
  qnought, sigmaprime, intcouplings, testpotentials )

  % We're approximating the potentials using:
  % potentials = intcouplings * Q_0 * exp( potentials / sigmaprime )
  % The error is the amount by which this is off, for each component.

  scratch = exp( testpotentials / sigmaprime );
  scratch = intcouplings * qnought * scratch;

  solutionerrors = scratch - testpotentials;

end



%
% This is the end of the file.
