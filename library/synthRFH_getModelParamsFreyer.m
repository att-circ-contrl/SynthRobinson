function [ modelparams intcouplings ] = synthRFH_getModelParamsFreyer()

% function [ modelparams intcouplings ] = synthRFH_getModelParamsFreyer()
%
% This returns model and coupling parameters for use with
% synthRFH_stepCortexThalamus() and related functions.
%
% Values are the ones used in Freyer 2011 (Table 1):
% https://www.jneurosci.org/content/31/17/6353.short
%
% These have slightly adjusted alpha and beta time constants and
% substantially different coupling weights vs Robinson 2002.
%
% No arguments.
%
% "modelparams" is a model parameter structure with the fields described
%   in MODELPARAMS.txt.
% "intcouplings" is a 4x4 matrix indexed by (destination,source) that
%   provides the coupling weights (in mV*s) between excitatory, inhibitory,
%   specific nucleus, and reticular nucleus neurons.


modelparams = struct();

% Paramters for converting potentials to firing rates (sigmoid parameters).
modelparams.qmax = 250;         % 1/sec
modelparams.threshlevel = 15;   % mV
modelparams.threshsigma = 6;    % mV

% Parameters for neural population dynamics (second order DE weights).
modelparams.alpha = 60;    % 1/sec
modelparams.beta = 240;    % 1/sec
modelparams.gamma = 100;   % 1/sec

% Parameters for cortico-thalamic circuit dynamics.
modelparams.halfdelay_ms = 40;   % ms

% Noise parameters.
% FIXME - Mean isn't documented!
modelparams.noisemean = 0;            % mV
modelparams.noisesigma = 0.56;        % mV
modelparams.noisemultfactor = 0.64;   % dimensionless


% Internal coupling parameters.

intcouplings = ...
[ 1.06 -1.8 2.2   0 ; ...
  1.06 -1.8 2.2   0 ; ...
  2.28  0   0    -0.845 ; ...
  0.91  0   0.41  0 ];


% Additional coupling parameters.
% Noise couples to the specific nucleus, population mixtures to excitatory
% cortical neurons.

modelparams.noisecoupling = 1.2;
modelparams.mixturecoupling = 0;  % Only Hindriks used population mixing.


% Done.
end


%
% This is the end of the file.
