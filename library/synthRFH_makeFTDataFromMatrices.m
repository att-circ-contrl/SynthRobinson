function ftdata = synthRFH_makeFTDataFromMatrices( ...
  wavedata, samprate, trigoffsetsamps, trigtimes, labelprefix )

% function ftdata = synthRFH_makeFTDataFromMatrices( ...
%   wavedata, samprate, trigoffsetsamps, trigtimes, labelprefix )
%
% This builds a ft_datatype_raw structure containing supplied waveform data.
%
% "wavedata" is either a Nchans x Nsamples x Ntrials matrix or a cell array
%   with Ntrials cells, each containing a Nchans x Nsamples matrix.
% "samprate" is the sampling rate.
% "trigoffsetsamps" specifies when the trigger time is within each trial.
%   This is 0 if the first sample in the trial is at the trigger, positive
%   if a later sample in the trial is the trigger, and negative if the
%   trigger occurred before the first sample in the trial.
% "trigtimes" is a vector with trigger times for each trial. This is used for
%   constructign trial definitions. If necessary, an offset is added to
%   guarantee that all samples have positive indices in the global data.
%   Specify [] to automatically build trigger times for trial definitions.
% "labelprefix" is a character vector containing a prefix to use when
%   constructing channel labels.
%
% "ftdata" is a ft_datatype_raw structure containing trial data. This
%   includes a header and a config structure with cfg.trl.


% Cheat, and convert a cell array to a matrix, to guarantee consistent
% sizes.

if iscell(wavedata)
  chancount = inf;
  sampcount = inf;
  trialcount = length(wavedata);

  % Find the minimum record size.
  for tidx = 1:trialcount
    thismatrix = wavedata{tidx};

    thischancount = size(thismatrix,1);
    thissampcount = size(thismatrix,2);

    chancount = min(chancount, thischancount);
    sampcount = min(sampcount, thissampcount);
  end

  scratch = zeros(chancount, sampcount, trialcount);

  % Crop all records to the minimum size and copy them.
  for tidx = 1:trialcount
    thismatrix = wavedata{tidx};
    scratch(1:chancount,1:sampcount,tidx) = ...
      thismatrix(1:chancount,1:sampcount);
  end

  wavedata = scratch;
  scratch = [];
end


% Convert the matrix to a cell array of trials.

chancount = size(wavedata,1);
sampcount = size(wavedata,2);
trialcount = size(wavedata,3);

trialdata = {};

for tidx = 1:trialcount
  trialdata{tidx} = reshape( wavedata(:,:,tidx), chancount, sampcount );
end


% Build the canonical time series.

timeseries = 0:(sampcount-1);
timeseries = (timeseries - trigoffsetsamps) / samprate;

timecells = {};
for tidx = 1:trialcount
  timecells{tidx} = timeseries;
end


% Get labels.

chanlabels = {};

for cidx = 1:chancount
  chanlabels{cidx} = sprintf( '%s_%03d', labelprefix, cidx );
end

if ~iscolumn(chanlabels)
  chanlabels = transpose(chanlabels);
end



%
% Build trial definitions.


% Get trial start locations.

if isempty(trigtimes)
  trigtimes_samps = 0:(trialcount-1);
  trigtimes_samps = trigtimes_samps * sampcount;
else
  trigtimes_samps = round(trigtimes * samprate);
end

% Earliest trigger sample index.
scratch = min(trigtimes_samps);

% Earliest start-of-trial sample index.
scratch = scratch - trigoffsetsamps;

% FIXME - Assume Field Trip counts samples starting at 1, rather than 0.
% Make sure the earliest sample is no earlier than that.
if scratch < 1
  trigtimes_samps = trigtimes_samps + 1 - scratch;
end


% Build the rest of the trial definition matrix.
% Tuples are [ start sample, end sample, trig location relative to start ].

trialdefs = zeros(trialcount, 3);

trialdefs(:,1) = trigtimes_samps;
trialdefs(:,2) = trigtimes_samps + sampcount - 1;
trialdefs(:,3) = trigoffsetsamps;



%
% Assemble FT structures.


ftheader = struct();

ftheader.Fs = samprate;
ftheader.nChans = chancount;
ftheader.nSamples = sampcount;
ftheader.nSamplesPre = trigoffsetsamps;
ftheader.nTrials = trialcount;
ftheader.label = chanlabels;


ftconfig = struct();

ftconfig.trl = trialdefs;


ftdata = struct();

ftdata.label = chanlabels;
ftdata.time = timecells;
ftdata.trial = trialdata;

ftdata.sampleinfo = trialdefs(:,1:2);
ftdata.fsample = samprate;

ftdata.hdr = ftheader;
ftdata.cfg = ftconfig;



% Done.
end


%
% This is the end of the file.
