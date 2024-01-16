% Example Code - Data synthesis using the Robinson/Freyer/Hindriks model.
% Written by Christopher Thomas.


%
% Libraries and folders.


% Generated files get put in this folder.

outdir = 'output';


% Library paths.

% NOTE - These are specific to my test system; as long as you have the
% SynthRobinson library path and the Field Trip library path set up, you can
% remove these.
addpath('../../library');
addpath('../../../../neurolab/fieldtrip/fieldtrip-latest');


% Field Trip setup.

evalc('ft_defaults');
ft_notice('off');
ft_info('off');
ft_warning('off');



%
% Configuration.


% Number of trials per test configuration.
trialcount = 10;

% Trial duration before and after "stimulation".
% Signal perturbations take about 1 second to settle after parameter changes.
trial_secs_before = 5;
trial_secs_after = 10;


% Sampling rates and filter parameters.
% Remember that we're generating firing rates, not an LFP. This pretty much
% _is_ rectified multi-unit activity.
samprate_sim = 10000;
samprate_mua = 2000;
mua_lowpass = 50; % Hz.


% Simulation parameters.

chancount = 4;

sim_startup_secs = 2;

% Explicity set this, rather than keeping the Hindriks baseline value.
% 0.07 gets us strong power spectrum peaks at harmonics, not just the
% fundamental mode, but can also go unstable and stay at the maximum rate.
mixture_coupling_coeff = 0.06;

% Factor by which to increase v_es when boosting activity. This is very
% sensitive; raising 1% or decreasing 2% has a very visible effect.
mua_enhance_coeff_factor = 1.01;
mua_suppress_coeff_factor = 0.98;

% Mixing is already gamma-delayed, so an extra delay isn't necessary.
mixing_delays_ms = zeros(chancount, chancount);


% Mixing matrices for various test cases.

% We have to keep each population's input consistent between cases, since
% behavior is very sensitive to it. Do this by having rows sum to consistent
% values but take contributions from different population distributions.

% Everything couples to everything. All links bidirectional.
mixing_matrix_uniform = ...
[ 0     0.333 0.333 0.333 ; ...
  0.333 0     0.333 0.333 ; ...
  0.333 0.333 0     0.333 ; ...
  0.333 0.333 0.333 0     ];

% 1-2 and 3-4, bidirectional.
mixing_matrix_routing_A = ...
[ 0 1 0 0 ; ...
  1 0 0 0 ; ...
  0 0 0 1 ; ...
  0 0 1 0 ];

% 1-4 and 2-3, bidirectional.
mixing_matrix_routing_B = ...
[ 0 0 0 1 ; ...
  0 0 1 0 ; ...
  0 1 0 0 ; ...
  1 0 0 0 ];

% 1-2-3-4-1, bidirectional.
mixing_matrix_loop_bidirectional = ...
[ 0   0.5 0   0.5 ; ...
  0.5 0   0.5 0   ; ...
  0   0.5 0   0.5 ; ...
  0.5 0   0.5 0   ];

% 1->2->3->4->1, directional.
mixing_matrix_loop_ascending = ...
[ 0 0 0 1 ; ...
  1 0 0 0 ; ...
  0 1 0 0 ; ...
  0 0 1 0 ];

% 1<-2<-3<-4<-1, directional.
mixing_matrix_loop_descending = ...
[ 0 1 0 0 ; ...
  0 0 1 0 ; ...
  0 0 0 1 ; ...
  1 0 0 0 ];



%
% Generate and save the synthetic trials.


% Get baseline simulation parameters.

[ modelparams_baseline intcouplings_baseline ] = ...
  synthRFH_getModelParamsHindriks;

% Override the default mixture coupling.
modelparams_baseline.mixturecoupling = mixture_coupling_coeff;


% Set up two epochs (before and after stimulation).

durationlist = [ trial_secs_before, trial_secs_after ];

modelparams = modelparams_baseline;
modelparams(2) = modelparams_baseline;

intcouplings = intcouplings_baseline;
intcouplings(:,:,2) = intcouplings_baseline;

% We never change this, so just make the original dual-epoch.
mixing_delays_ms(:,:,2) = mixing_delays_ms(:,:,1);



% Build test cases using two different baselines: uniform, and bidirectional
% loop.

% Uniform isn't actually that useful as a baseline, but include it anyways.
% The bidirectional loop is equal to the two routing states superimposed
% and to the two directed loop cases superimposed, so it's a better baseline.

baselinematrices = ...
  { mixing_matrix_uniform, mixing_matrix_loop_bidirectional };
baselinelabels = { 'unibase', 'loopbase' };
baselinetitles = { 'uniform', 'bidirected loop' };

muafactors = [ 1.0, mua_enhance_coeff_factor, mua_suppress_coeff_factor ];
mualabels = { 'baseline', 'strong', 'weak' };
muatitles = { 'baseline', 'stronger', 'weaker' };

mixcasematrices = { mixing_matrix_routing_A, mixing_matrix_routing_B, ...
  mixing_matrix_loop_ascending, mixing_matrix_loop_descending };
mixcaselabels = { 'routeA', 'routeB', 'loopup', 'loopdown' };
mixcasetitles = { 'routing state A', 'routing state B', ...
  'ascending loop', 'descending loop' };

[ idxlut namelut ] = synthRFH_getRegionInfo();

totaltime = 0;

for baseidx = 1:length(baselinematrices)

  disp([ '== Generating cases with ' baselinetitles{baseidx} ' baseline.' ]);

  mixing_matrix = baselinematrices{baseidx};
  mixing_matrix(:,:,2) = mixing_matrix(:,:,1);


  % Generate different strengths of MUA response. This includes the baseline.

  ves_baseline = intcouplings(idxlut.cortex_excitatory, ...
    idxlut.thalamus_specific, 2);

  for caseidx = 1:length(muafactors)
    disp([ '-- Generating ' muatitles{caseidx} ' MUA response.' ]);
    tic;

    intcouplings(idxlut.cortex_excitatory, idxlut.thalamus_specific, 2) = ...
      ves_baseline * muafactors(caseidx);

    ftdata = synthRFH_ftWrapper_simulateNetwork( ...
      trialcount, trial_secs_before, {}, true, ...
      durationlist, sim_startup_secs, 1 / samprate_sim, modelparams, ...
      intcouplings, chancount, mixing_matrix, mixing_delays_ms );

    fname = [ outdir filesep 'mua-' baselinelabels{baseidx} ...
      '-' mualabels{caseidx} '.mat' ];
    helper_writeFTData( fname, ftdata, modelparams, intcouplings, ...
      mixing_matrix, mixing_delays_ms, mua_lowpass, samprate_mua );

    totaltime = totaltime + toc;
    durstring = helper_makePrettyTime(toc);
    disp([ '.. Elapsed time: ', durstring ]);
  end

  % Restore the original couplings.
  intcouplings(:,:,2) = intcouplings(:,:,1);


  % Generate different routing states. Baseline is already covered.

  for caseidx = 1:length(mixcasematrices)
    disp([ '-- Generating response with ' mixcasetitles{caseidx} ...
      ' mixing.' ]);

    mixing_matrix(:,:,2) = mixcasematrices{caseidx};

    ftdata = synthRFH_ftWrapper_simulateNetwork( ...
      trialcount, trial_secs_before, {}, true, ...
      durationlist, sim_startup_secs, 1 / samprate_sim, modelparams, ...
      intcouplings, chancount, mixing_matrix, mixing_delays_ms );

    fname = [ outdir filesep 'mua-' baselinelabels{baseidx} ...
      '-' mixcaselabels{caseidx} '.mat' ];
    helper_writeFTData( fname, ftdata, modelparams, intcouplings, ...
      mixing_matrix, mixing_delays_ms, mua_lowpass, samprate_mua );

    totaltime = totaltime + toc;
    durstring = helper_makePrettyTime(toc);
    disp([ '.. Elapsed time: ', durstring ]);
  end

  % Restore the original mixing matrix.
  mixing_matrix(:,:,2) = mixing_matrix(:,:,1);

end

disp('== Finished generating cases.');

durstring = helper_makePrettyTime(totaltime);
disp([ '== Total elapsed time: ', durstring ]);



%
% Helper Functions


% This writes a Field Trip dataset and auxiliary metadata to a matlab file.
% The Field Trip dataset is filtered and downsampled before writing.

function helper_writeFTData( fname, ftdata_mua, modelparams, intcouplings, ...
  mixing_matrix, mixing_delays_ms, lowpass_corner, resample_rate )

  % Field Trip configuration structures.
  ftconfig_filt = struct( 'lpfilter', 'yes', 'lpfilttype', 'but', ...
    'lpfreq', lowpass_corner, 'feedback', 'no' );
  ftconfig_resample = struct( 'resamplefs', resample_rate, ...
    'detrend', 'no', 'feedback', 'no' );

  % Extract native-rate information that we want to keep.
  % NOTE - We know that synthRFH_makeFTDataFromMatrices() provides these.
  origheader = ftdata_mua.hdr;
  origconfigtrl = ftdata_mua.cfg.trl;

  % Filter and downsample.
  ftdata_mua = ft_preprocessing( ftconfig_filt, ftdata_mua );
  ftdata_mua = ft_resampledata( ftconfig_resample, ftdata_mua );

  % Save data and metadata.
  save( fname, 'ftdata_mua', 'origheader', 'origconfigtrl', ...
    'modelparams', 'intcouplings', 'mixing_matrix', 'mixing_delays_ms', ...
    '-v7.3' );

end



%
% This is the end of the file.
