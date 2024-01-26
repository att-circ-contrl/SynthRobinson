function [ rategradients potentialgradients ] = ...
  synthRFH_getOperatingPointGradient( ...
    modelparams, intcouplings, testpotentials, couplingstep, zerohandling )

% function [ rategradients potentialgradients ] = ...
%   synthRFH_getOperatingPointGradient( ...
%     modelparams, intcouplings, testpotentials, couplingstep, zerohandling )
%
% This function attempts to estimate the gradient of the DC operating point
% of a Robinson neural model with respect to the internal coupling matrix.
%
% This assumes firing rates that are much less than the maximum rate. Per the
% model guide, operating points under these conditions are solutions to:
%
% potentials = intcouplings * Q_0 * exp( potentials / sigmaprime )
%
% The gradient of this function is evaluated numerically, by perturbing the
% coupling matrix and finding operating points for each perturbed version.
%
% "modelparams" is a model parameter structure with the fields described in
%   MODELPARAMS.txt.
% "intcouplings" is a 4x4 matrix indexed by (destination, source) that
%   provides the coupling weights (in mV*s) between excitatory, inhibitory,
%   specific nucleus, and reticular nucleus neural populations.
% "testpotentials" is a vector containing cell potentials for the excitatory,
%   inhibitory, specific nucleus, and reticular nucleus populations at the
%   test point to analyze. This is used as the starting point for the
%   operating point search.
% "couplingstep" is a scalar indicating the amount by which each coupling
%   value should be perturbed when evaluating the gradient numerically. This
%   should be much smaller than the magnitude of nonzero coupling values.
% "zerohandling" is 'all' to compute the gradient with respect to all
%   coupling values, and 'nonzero' to compute the gradient with respect to
%   all coupling values with magnitudes larger than "couplingstep" (storing
%   NaN as the gradient for coupling values that are smaller than this).
%
% "rategradients" is a cell array with 4 cells, containing gradient matrices
%   with respect to "intcouplings" for the excitatory, inhibitory, specific
%   nucleus, and reticular nucleus firing rates.
% "potentialgradients" is a cell array with 4 cells, containing gradient
%   matrices with respect to "intcouplings" for the excitatory, inhibitory,
%   specific nucleus, and reticular nucleus cell potentials.


% Get magic values.

regioncount = size(intcouplings,1);

sigmaprime = modelparams.threshsigma * sqrt(3) / pi;
qnought = modelparams.qmax * exp( - modelparams.threshlevel / sigmaprime );



% Make sane default output.

rategradients = {};
potentialgradients = {};

for outidx = 1:regioncount
  rategradients{outidx} = NaN(regioncount,regioncount);
  potentialgradients{outidx} = NaN(regioncount,regioncount);
end



% Iterate across the coupling dimensions.
% Output iteration happens in an inner loop.

for dstidx = 1:regioncount
  for srcidx = 1:regioncount

    thisweight = intcouplings(dstidx,srcidx);

    wantcalc = true;
    if strcmp('nonzero', zerohandling)
      wantcalc = abs(thisweight) >= couplingstep;
    end

    if wantcalc

      % Recompute the operating point for each perturbed location.

      couplingsneg = intcouplings;
      couplingsneg(dstidx,srcidx) = thisweight - couplingstep;
      couplingspos = intcouplings;
      couplignspos(dstidx,srcidx) = thisweight + couplingstep;

      [ scratch potentialsneg ] = ...
        synthRFH_estimateOperatingPointExponential( ...
          modelparams, couplingsneg, testpotentials );
      [ scratch potentialspos ] = ...
        synthRFH_estimateOperatingPointExponential( ...
          modelparams, couplingspos, testpotentials );


      % Get derivatives by comparing these operating points.

      thispotentialderiv = ...
        (potentialspos - potentialsneg) / (2 * couplingstep);

      % Chain rule.
      dQdV = synthRFH_getSigmoidDerivative( ...
        testpotentials, modelparams.qmax, ...
        modelparams.threshlevel, modelparams.threshsigma );
      thisratederiv = thispotentialderiv .* dQdV;

      for outidx = 1:regioncount
        rategradients{outidx}(dstidx,srcidx) = ...
          thisratederiv(outidx);
        potentialgradients{outidx}(dstidx,srcidx) = ...
          thispotentialderiv(outidx);
      end

    end

  end
end



% Done.
end


%
% This is the end of the file.
