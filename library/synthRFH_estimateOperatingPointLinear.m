function [ firingrates potentials ] = ...
  synthRFH_estimateOperatingPointLinear( modelparams, intcouplings )

% function [ firingrates potentials ] = ...
%   synthRFH_estimateOperatingPointLinear( modelparams, intcouplings )
%
% This attempts to estimate the DC operating point of a Robinson neural model.
%
% Per the model guide, operating points with firing rates much less than the
% maximum are solutions to the equation:
%
% potentials = intcouplings * Q_0 * exp( potentials / sigmaprime )
%
% This function uses a linear approximation to exp(x) to estimate operating
% points for potentials that are small compared to sigmaprime. NOTE - This
% is not a robus assumption! The operating point _must_ be examined to
% confirm that this condition holds. If it doesn't hold, the estimated
% operating point is not correct.
%
% "modelparams" is a model parameter structure with the fields described in
%   MODELPARAMS.txt.
% "intcouplings" is a 4x4 matrix indexed by (destination,source) that
%   provides the coupling weights (in mV*s) between excitatory, inhibitory,
%   specific nucleus, and reticular nucleus neural populations.
%
% "firingrates" is a vector containing firing rates for the excitatory,
%   inhibitory, specific nucleus, and reticular nucleus populations.
% "potentials" is a vector containing cell potentials for the excitatory,
%   inhibitory, specific nucleus, and reticular nucleus populations.

% Per the model guide, we're approximating V as:
% V = N Q_0 [ vec(1) + (1/sigmaprime) V ]

% Solving for V gives:
% V = [ (1/Q_0) N^-1 - (1/sigmaprime) I ]^-1 * vec(1)
% But that doesn't work, since the matrix inverse can fail.

% It can also be solved to have the form Ax = b:
% [ (1/Q_0) I - (1/sigmaprime) N ] V = N * \vec(1)


sigmaprime = modelparams.threshsigma * sqrt(3) / pi;

qnought = modelparams.qmax * exp( - modelparams.threshlevel / sigmaprime );

popcount = size(intcouplings,1);

identity = eye(popcount);
onesvec = ones(popcount,1);

matrix_a = identity / qnought - intcouplings / sigmaprime;
vector_b = intcouplings * onesvec;


% Use Matlab's helper to solve this linear system.
% This handles special cases better than matrix inverse.

potentials = matrix_a \ vector_b;

firingrates = synthRFH_getSigmoid( potentials, ...
  modelparams.qmax, modelparams.threshlevel, modelparams.threshsigma );


% FIXME - Sanity check, to catch unsolvable systems.

scratch = potentials / sigmaprime + onesvec;
scratch = intcouplings * qnought * scratch;
solutionerror = norm(scratch - potentials) / norm(potentials);

if (~isfinite(solutionerror)) || (solutionerror > 0.03)
  disp('###  Couldn''t solve linear operating point equations!');
  disp(sprintf( 'Reconstruction error:  %.1f %%', 100 * solutionerror ));
end


% Done.
end


%
% This is the end of the file.
