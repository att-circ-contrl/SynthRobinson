function edgegaingradients = synthRFH_getEdgeGainGradients( ...
  modelparams, intcouplings, firingrates, rategradients )

% function edgegaingradients = synthRFH_getEdgeGainGradients( ...
%   modelparams, intcouplings, firingrates, rategradients )
%
% This function calculates the gradient of the small-signal gains of each
% network edge in a Robinson neural model with respect to the internal
% coupling matrix, at a specified operating point.
%
% This assumes firing rates that are much less than the maximum rate, and
% assumes that the gradients of the firing rates are already known.
%
% "modelparams" is a model parameter structure with the fields described in
%   MODELPARAMSROBINSON.txt.
% "intcouplings" is a 4x4 matrix indexed by (destination, source) that
%   provides the coupling weights (in mV*s) between excitatory, inhibitory,
%   specific nucleus, and reticular nucleus neural populations.
% "firingrates" is a vector specifying the operating point firing rates of
%   excitatory, inhibitory, specific nucleus, and reticular nucleus neural
%   populations.
% "rategradients" is a cell array with 4 cells, containing gradient matrices
%   with respect to "intcouplings" for the excitatory, inhibitory, specific
%   nucleus, and reticular nucleus firing rates.
%
% "edgegaingradients" is a 4x4 cell array indexed by (destination, source)
%   that contains the gradient with respect to "intcouplings" of the
%   small-signal firing rate gains between each source and destination for
%   the excitatory, inhibitory, specific nucleus, and reticular nucleus
%   neural populations.


% Get relevant metadata.

regioncount = size(intcouplings,1);

sigmaprime = modelparams.threshsigma * sqrt(3) / pi;
qmax = modelparams.qmax;


% From the manual:
% del_N[ G_ab ] = Q'(V_a) del_N[ nu_ab ] + nu_ab del_N[ Q'(V_a) ]
% del_N[ nu_ab ] is a zero matrix with a single cell equal to 1.
% del_N[ Q'(V_a) ] = (1/sigmaprime) (1 - 2*phi_a/Q_max) del_N[ phi_a ]
% ...and we're given the gradient of phi_a as one of the arguments.


potentials = synthRFH_getSigmoidInverse( firingrates, ...
  modelparams.qmax, modelparams.threshlevel, modelparams.threshsigma );

qprime = synthRFH_getSigmoidDerivative( potentials, ...
  modelparams.qmax, modelparams.threshlevel, modelparams.threshsigma );


edgegaingradients = {};

for dstidx = 1:regioncount

  % Most of the factors in each term just depend on source, not destination.

  thisqprime = qprime(dstidx);

  thisratephi = firingrates(dstidx);
  thisrategrad = rategradients{dstidx};

  % Compute nu_ab del_N[ Q'(V_a) ].
  gradfactor = 2 * thisratephi / qmax;
  gradfactor = (1 - gradfactor) / sigmaprime;
  gradfactor = thisrategrad * gradfactor;

  for srcidx = 1:regioncount

    thisweight = intcouplings(dstidx,srcidx);

    thisgaingrad = thisweight * gradfactor;
    thisgaingrad(dstidx,srcidx) = thisgaingrad(dstidx,srcidx) + thisqprime;

    edgegaingradients{dstidx,srcidx} = thisgaingrad;

  end
end




% Done.
end


%
% This is the end of the file.
