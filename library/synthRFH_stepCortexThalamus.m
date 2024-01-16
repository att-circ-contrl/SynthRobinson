function statefuture = synthRFH_stepCortexThalamus( modelparams, ...
  timestep, statepresent, statepast, intcouplings, extrates, extcouplings )

% function statefuture = synthRFH_stepCortexThalamus( modelparams, ...
%   timestep, statepresent, statepast, intcouplings, extrates, extcouplings )
%
% This generates the future state of a cortico-thalamic loop given the
% present state, using the model from Robinson 2002 with augmented input
% per Freyer 2011 and Hindriks 2023:
%
% https://journals.aps.org/pre/abstract/10.1103/PhysRevE.65.041924
% https://www.jneurosci.org/content/31/17/6353.short
% https://www.nature.com/articles/s42003-023-04648-x
%
% This simulates N populations of excitatory and inhibatory neurons in the
% cortex, and N populations of neurons in the specific and reticular nuclei
% in the thalamus (per Freyer 2011). Population bins are independent of
% each other (communication is handled by the caller).
%
% NOTE - This uses the forward Euler method for evolving system state.
% This is numerically stable if and only if the time step is much smaller
% than the time scales of any system dynamics. Make this much smaller than
% you think you need to.
%
% "modelparams" is a structure with the following fields (as described in
%   MODELPARAMSROBINSON.txt):
%   "qmax" is the maximum firing rate (1/sec); typically 250.
%   "threshlevel" is the average neuronal threshold (mV); typ. 15.
%   "threshsigma" is the standard deviation of the neuronal threshold (mV);
%     typically 6.0.
%   "alpha" is the inverse decay time (1/sec); typically 50.
%   "beta" is the inverse rise time (1/sec); typically 200.
%   "gamma" is the inverse within-cortex propagation time (1/sec); typ. 100.
% "timestep" is the amount of time to advance the simulation by (in seconds).
%   NOTE - This must be much smaller than system dynamics timescales!
% "statepresent" is a structure with the following fields:
%   "potentials" is a 4xN matrix containing Ve_k, Vi_k, Vs_k, and Vr_k from
%     the Robinson model.
%   "velocities" is a 4xN matrix containing the approximate first time
%     derivative of "potentials".
%   "cortexrates" is a 1xN matrix containing the gamma-damped firing rate
%     of cortex excitatory neurons (phi_e).
%   "cortexvelocities" is a 1xN matrix containing the approxmate first time
%     derivative of "cortexrates".
% "statepast" is a structure with the same fields as "statepresent", taken
%   from a past stimulation step. This should be delayed by the one-way
%   cortex/thalamus communication time (one half of the round-trip delay).
% "intcouplings" is a 4x4 matrix indexed by (destination,source) that
%   provides the coupling weights (in mV*s) between excitatory cortex
%   neurons (1), inhibitory cortex neurons (2), specific nucleus neurons (3),
%   and reticular nucleus neurons (4). Typical couplings range from -2 to +2.
% "extrates" is a MxN matrix containing the firing rates of M external inputs
%   to the system. These may represent noise (as with Freyer 2011) or
%   cross-coupled activity within the cortex (as with Hindriks 2023).
%   This may be empty (M = 0).
% "extcouplings" is a 4xM matrix indexed by (destination,source) that
%   provides the coupling weights (in mV*s) between the internal model
%   neurons (destinations) and the system's external inputs (sources).
%   This may be empty (M = 0).
%
% "statefuture" is a copy of "statepresent" advanced by one time step.


% Copy relevant data, for convenience.

potentials = statepresent.potentials;
velocities = statepresent.velocities;

pastpotentials = statepast.potentials;

cortexrates = statepresent.cortexrates;
cortexvelocities = statepresent.cortexvelocities;

pastcortexrates = statepast.cortexrates;



%
% Get firing rates. We need both V_a and phi_a.

% NOTE - We're special-casing signals sourced from cortex excitatory neurons.
% These use the gamma-delayed firing rates, rather that the rate computed
% from Ve.
%
% The references only do this for excitatory cortex neurons, not inhibitory
% cortex neurons or thalamus neurons.

firingrates = synthRFH_getSigmoid( potentials, ...
  modelparams.qmax, modelparams.threshlevel, modelparams.threshsigma );

pastfiringrates = synthRFH_getSigmoid( pastpotentials, ...
  modelparams.qmax, modelparams.threshlevel, modelparams.threshsigma );

% Copy the gamma-delayed rates as phi_e.
% Save the non-delayed rate, as we need it to compute the delayed rate.
immedcortexrates = firingrates(1,:);
firingrates(1,:) = cortexrates;
pastfiringrates(1,:) = pastcortexrates;

% Make note of our population size.
popcount = length(cortexrates);

% Make note of our external signal count. This may be zero!
extcount = size(extrates,1);



%
% Get second derivatives using the model.

% Per Robinson 2002, this uses the following equation (where ' is the time
% derivative):
%
% V''_a = (alpha * beta) [ -V_a + sum_b[ weight_ab * phi_b ] ]
%    - (alpha + beta) V'_a
%
% phi_b is either the present or past firing rate of area b, depending
% on whether a and b are both in the cortex/thalamus or if one is in each.
%
% We're adding external inputs as additional regions, per Hindriks 2023.
%
% Per Robinson 2002, this uses the following delay equation for the
% excitatory firing rates in the cortex (phi_eout):
%
% phi''_eout = gamma^2 [ - phi_eout + phi_e ] - 2 * gamma * phi'_eout


cortexindices = [ 1 2 ];


% First, build the term that's inside the square brackets above. This
% has the present potential and the firing rate contributions.

% Leading -V_a term.
accelpotentials = - potentials;

% Internal signals.

for dstidx = 1:4
  % Remaining terms for internal phi_b.
  for srcidx = 1:4
    if ismember(dstidx, cortexindices) == ismember(srcidx, cortexindices)
      % Same region. Use the present state.
      accelpotentials(dstidx,:) = accelpotentials(dstidx,:) ...
        + intcouplings(dstidx,srcidx) * firingrates(srcidx,:);
    else
      % Different region. Use the delayed state.
      accelpotentials(dstidx,:) = accelpotentials(dstidx,:) ...
        + intcouplings(dstidx,srcidx) * pastfiringrates(srcidx,:);
    end
  end
end

% External signals (external phi_b terms).

for dstidx = 1:4
  for srcidx = 1:extcount
    accelpotentials(dstidx,:) = accelpotentials(dstidx,:) ...
      + extcouplings(dstidx,srcidx) * extrates(srcidx,:);
  end
end


% Next, multiply by (alpha * beta) and add the velocity term at the end.

accelpotentials = accelpotentials * modelparams.alpha * modelparams.beta ...
  - velocities * ( modelparams.alpha + modelparams.beta );


% Finally, get the damped/delayed cortex excitatory firing rates.

accelcortexrates = ...
  modelparams.gamma * modelparams.gamma * (immedcortexrates - cortexrates) ...
  - 2 * modelparams.gamma * cortexvelocities;



%
% Generate the future state from the present state and acceleration.

% NOTE - This is only numerically stable if the time step is very much
% smaller than the time scale of system dynamics!


% Blind copy to keep any user-added fields.
statefuture = statepresent;

% Forward Euler method. Easy to apply but has the stability issue noted.

statefuture.potentials = potentials + timestep * velocities;
statefuture.velocities = velocities + timestep * accelpotentials;

statefuture.cortexrates = cortexrates + timestep * cortexvelocities;
statefuture.cortexvelocities = ...
  cortexvelocities + timestep * accelcortexrates;



% Done.
end


%
% This is the end of the file.
