function edgegains = synthRFH_getEdgeGains( ...
  modelparams, intcouplings, firingrates )

% function edgegains = synthRFH_getEdgeGains( ...
%   modelparams, intcouplings, firingrates )
%
% This function estimates the small-signal gain of each network edge in a
% Robinson neural model, at a specified operating point.
%
% "modelparams" is a model parameter structure with the fields described in
%   MODELPARAMS.txt.
% "intcouplings" is a 4x4 matrix indexed by (destination, source) that
%   provides the coupling weights (in mV*s) between excitatory, inhibitory,
%   specific nucleus, and reticular nucleus neural populations.
% "firingrates" is a vector specifying the operating point firing rates of
%   excitatory, inhibitory, specific nucleus, and reticular nucleus neural
%   populations.
%
% "edgegains" is a 4x4 matrix indexed by (destination, source) that contains
%   the small-signal firing rate gains between each source and destination
%   for the excitatory, inhibitory, specific nucleus, and reticular nucleus
%   neural populations.


regioncount = size(intcouplings,1);
sigmaprime = modelparams.threshsigma * sqrt(3) / pi;
qmax = modelparams.qmax;


% From the manual, and from Robinson 2002:
% G_ab = ( nu_ab / sigmaprime ) phi_a ( 1 - phi_b/Q_max )

edgegains = zeros(regioncount,regioncount);

for srcidx = 1:regioncount
  for dstidx = 1:regioncount
    dstrate = firingrates(dstidx);
    edgegains(dstidx,srcidx) = intcouplings(dstidx,srcidx) ...
      * dstrate * (1 - (dstrate/qmax)) / sigmaprime;
  end
end


% Done.
end


%
% This is the end of the file.
