function [ modelparams intcouplings ] = synthRFH_getModelParamsHindriks()

% function [ modelparams intcouplings ] = synthRFH_getModelParamsHindriks()
%
% This returns model and coupling parameters for use with
% synthRFH_stepCortexThalamus() and related functions.
%
% Values are the ones used in Hindriks 2023:
% https://www.nature.com/articles/s42003-023-04648-x
% https://github.com/Prejaas/amplitudecoupling
%
% These are identical to the Robinson 2002 values, except with a smaller
% noise coupling coefficient.
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
modelparams.alpha = 50;    % 1/sec
modelparams.beta = 200;    % 1/sec
modelparams.gamma = 100;   % 1/sec

% Parameters for cortico-thalamic circuit dynamics.
modelparams.halfdelay_ms = 40;  % ms

% Noise parameters.
modelparams.noisemean = 0;   % Paper says 0.5, code says 0. Only 0 works.
modelparams.noisesigma = 0.1;
modelparams.noisemultfactor = 0.3;   % Effects above 0.2; swept from 0 to 0.4.


% Internal coupling parameters.

intcouplings = ...
[ 1.2 -1.8 1.2  0 ; ...
  1.2 -1.8 1.2  0 ; ...
  1.2  0   0   -0.8 ; ...
  0.4  0   0.2  0 ];


% Additional coupling parameters.
% Noise couples to the specific nucleus, population mixtures to excitatory
% cortical neurons.

modelparams.noisecoupling = 0.5;

% FIXME - This is documented as swept from 0.05 to 0.08 in the paper, but
% it's ambiguous as to whether that's the weight (used here) or the
% connectivity (in the mixing matrix).
% The code has a stubbed-out weight value of 0.84.

modelparams.mixturecoupling = 0.07;  % Effects above 0.06.


% Done.
end


%
% This is the end of the file.
