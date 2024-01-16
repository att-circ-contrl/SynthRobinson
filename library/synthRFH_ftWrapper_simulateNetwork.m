function ftdata = synthRFH_ftWrapper_simulateNetwork( ...
  trialcount, triggertime, poplabels, wantprogress, ...
  duration, startup, timestep, modelparams, intcouplings, ...
  popcount, cortexmixing, cortexdelays_ms )

% function ftdata = synthRFH_ftWrapper_simulateNetwork( ...
%   trialcount, triggertime, poplabels, wantprogress, ...
%   duration, startup, timestep, modelparams, intcouplings, ...
%   popcount, cortexmixing, cortexdelays_ms )
%
% This is a wrapper for synthRFH_simulateNetwork().
% See that function's documentation for details.
%
% This simulates cortex and thalamus neural activity, using the model from
% Robinson 2002 with augmented input per Freyer 2011 and Hindriks 2023:
%
% https://journals.aps.org/pre/abstract/10.1103/PhysRevE.65.041924
% https://www.jneurosci.org/content/31/17/6353.short
% https://www.nature.com/articles/s42003-023-04648-x
%
% "trialcount" is the number of Field Trip trials to simulate.
% "triggertime" is the trigger offset from the start of simulation (seconds).
% "poplabels" is a cell array containing Field Trip channel names for the
%   neural populations, or {} to automatically generate channel names.
% "wantprogress" is true to write progress messages to the console, false
%   otherwise.
%
% Remaining arguments are per synthRFH_simulateNetwork().
%
% "ftdata" is a ft_datatype_raw structure containing trial data, including
%   a header (hdr) and a config structure with trial definitions (cfg.trl).


% No starting banner; let the caller handle that.


trialwaves = {};

for tidx = 1:trialcount
  if wantprogress
    disp(sprintf( '.. Generating trial %d...', tidx ));
  end

  allrates = synthRFH_simulateNetwork( ...
    duration, startup, timestep, modelparams, intcouplings, ...
    popcount, cortexmixing, cortexdelays_ms );

  trialwaves{tidx} = reshape( allrates(1,:,:), popcount, [] );
end

if wantprogress
  disp('.. Aggregating.');
end

samprate = round( 1 / timestep );
ftdata = synthRFH_makeFTDataFromMatrices( trialwaves, samprate, ...
  round( triggertime * samprate ), [], 'CH' );

% Overwrite the default channel labels if we were given labels.
if ~isempty(poplabels)
  ftdata.hdr.label = poplabels;
  ftdata.label = chanlabels;
end


% No ending banner; let the caller handle that.


% Done.
end


%
% This is the end of the file.
