function [ modelparams intcouplings ] = synthRFH_getModelParamsRobinson()

% function [ modelparams intcouplings ] = synthRFH_getModelParamsRobinson()
%
% This returns model and coupling parameters for use with
% synthRFH_stepCortexThalamus() and related functions.
%
% Values are the ones used in Robinson 2002 (Table 1):
% https://journals.aps.org/pre/abstract/10.1103/PhysRevE.65.041924
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
% NOTE - Robinson didn't specify v_sn; he specified v_sn * phi_n (mean).
% NOTE - Robinson didn't vary noise; he used a constant background and
% then performed a perturbation analysis.
modelparams.noisemean = 1.0;
modelparams.noisesigma = 0.0;
modelparams.noisemultfactor = 0.0;  % Introduced by Freyer, so zero here.


% Internal coupling parameters.
% NOTE - Robinson set the inhibitory potential to be equal to the
% excitatory potential, rather than modelling inhibitory neurons separately.
% Freyer 2011 and Hindriks 2023 duplicated the couplings to get the same
% behavior.

intcouplings = ...
[ 1.2 -1.8 1.2  0 ; ...
  1.2 -1.8 1.2  0 ; ...
  1.2  0   0   -0.8 ; ...
  0.4  0   0.2  0 ];


% Additional coupling parameters.
% Noise couples to the specific nucleus, population mixtures to excitatory
% cortical neurons.

% NOTE - Robinson didn't specify v_sn; he specified v_sn * phi_n.
modelparams.noisecoupling = 1.0;
modelparams.mixturecoupling = 0;  % Only Hindriks used population mixing.


% Done.
end


%
% This is the end of the file.
