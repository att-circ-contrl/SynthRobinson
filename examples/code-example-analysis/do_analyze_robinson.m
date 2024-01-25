% Example Code - Analysis of behavior of the Robinson/Freyer/Hindriks model.
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


% Pick a model to test.

model_name = 'hindriks';
%model_name = 'freyer';


% Switch for measuring operating points by simulation. This takes time.

want_sim = true;

sim_duration_secs = 30;
sim_startup_secs = 2;
sim_samprate = 10000;


% Loop detection. Any coupling weight smaller than this is treated as zero.

minweight = 0.01;


% Joint optimization of loop gains by adjusting coupling weights.

want_adjust_weights = false;

adjust_taulimit = 1.0;
adjust_bestfactor = 1.2;

% Default is 200 * paramcount, or about 1400.
% Set NaN or [] to keep the default.
%adjust_maxprobes = 100;
adjust_maxprobes = 10000;

% Joint optimization goal depends on the model.
if strcmp('hindriks', model_name)
  % Test perturbing Hindriks.
  adjust_goals = { 'ES', 'grow' };
else
  % Test perturbing Freyer.
  adjust_goals = { 'ES', 'biggest' };
% adjust_goals = { 'ES', 'decay' };
end


% Flags for reporting gradient measurements to console.

tattle_rate_gradients = false;
tattle_edge_gradients = false;
tattle_loop_gradients = false;



%
% Set up models and get metadata.

if strcmp('hindriks', model_name)
  disp('-- Using Hindriks 2023 model parameters.');
  [ modelparams intcouplings ] = synthRFH_getModelParamsHindriks();
else
  disp('-- Using Freyer 2011 model parameters.');
  [ modelparams intcouplings ] = synthRFH_getModelParamsFreyer();
end

regioncount = size(intcouplings,1);

[ indices_lut names_lut ] = synthRFH_getRegionInfo();
labels_from_indices = { names_lut.title };
letters_from_indices = { names_lut.letter };



%
% Get operating points by various methods.


% Simulation.

if want_sim

  disp('-- Measuring operating point by simulation.');
  tic;

  % Only simulate one population with no cross-coupling.
  popcount = 1;

  allrates = synthRFH_simulateNetwork( ...
    sim_duration_secs, sim_startup_secs, 1 / sim_samprate, ...
    modelparams, intcouplings, popcount, [], [] );

  allpotentials = synthRFH_getSigmoidInverse( allrates, ...
    modelparams.qmax, modelparams.threshlevel, modelparams.threshsigma );


  % NOTE - This assumes unimodal signals!
  % Signals that switch between several behaviors may give bad estimates.

  oprates = [];
  oppotentials = [];
  for ridx = 1:regioncount
    oprates(ridx,1) = mean(mean( allrates(ridx,:,:) ));
    oppotentials(ridx,1) = mean(mean( allpotentials(ridx,:,:) ));
  end


  durstring = helper_makePrettyTime(toc);
  disp([ '-- Simulation finished in ' durstring '.' ]);

  ratemsg = 'Real rates:              ';
  potmsg =  'Real potentials:        ';
  for ridx = 1:regioncount
    ratemsg = [ ratemsg sprintf('  %s: %.1f /sec', ...
      letters_from_indices{ridx}, oprates(ridx)) ];
    potmsg = [ potmsg sprintf('   %s: %.2f mV', ...
      letters_from_indices{ridx}, oppotentials(ridx)) ];
  end

  disp(ratemsg);
  disp(potmsg);

end


% Linear approximation.
% NOTE - This might wander outside the range of valid solutions!

[ oprates oppotentials ] = ...
  synthRFH_estimateOperatingPointLinear( modelparams, intcouplings );

ratemsg = 'Linear rates:            ';
potmsg =  'Linear potentials:      ';
for ridx = 1:regioncount
  ratemsg = [ ratemsg sprintf('  %s: %.1f /sec', ...
    letters_from_indices{ridx}, oprates(ridx)) ];
  potmsg = [ potmsg sprintf('   %s: %.2f mV', ...
    letters_from_indices{ridx}, oppotentials(ridx)) ];
end

disp(ratemsg);
disp(potmsg);


% Exponential approximation.
% NOTE - Fsolve works by black magic, so there are no guarantees about
% which solution it'll find out of the valid ones.

[ oprates oppotentials ] = ...
  synthRFH_estimateOperatingPointExponential( ...
    modelparams, intcouplings, [] );

ratemsg = 'Exponential rates:       ';
potmsg =  'Exponential potentials: ';
for ridx = 1:regioncount
  ratemsg = [ ratemsg sprintf('  %s: %.1f /sec', ...
    letters_from_indices{ridx}, oprates(ridx)) ];
  potmsg = [ potmsg sprintf('   %s: %.2f mV', ...
    letters_from_indices{ridx}, oppotentials(ridx)) ];
end

disp(ratemsg);
disp(potmsg);


% NOTE - We're keeping the exponential operating point for later use.



%
% Get loop information.


disp('-- Loop analysis:');

loopinfo = synthRFH_findLoops( modelparams, intcouplings, minweight );
edgegains = synthRFH_getEdgeGains( modelparams, intcouplings, oprates );
loopinfo = synthRFH_addLoopGainInfo( loopinfo, edgegains, {} );

[ looptext looptable ] = helper_makePrettyLoopTable(loopinfo);

disp(looptext);
save( [ outdir filesep 'loop-info-' model_name '.mat' ], ...
  'looptable', '-v7.3' );

disp('-- Finished loop analysis.');



%
% Gradient analysis.


disp('-- Computing gradients of potentials and firing rates.');

gradstep = 0.01;
zerohandling = 'nonzero';

[ rategrad potgrad ] = ...
  synthRFH_getOperatingPointGradient( ...
    modelparams, intcouplings, oppotentials, gradstep, zerohandling );

save( [ outdir filesep 'rate-gradients-' model_name '.mat' ], ...
  'rategrad', 'potgrad', '-v7.3' );

if tattle_rate_gradients

  disp('Firing rate gradients around exponential operating point:');

  for ridx = 1:regioncount
    disp([ '.. ' labels_from_indices{ridx} ':' ]);
    disp(rategrad{ridx});
  end

  disp('Potential gradients around exponential operating point:');

  for ridx = 1:regioncount
    disp([ '.. ' labels_from_indices{ridx} ':' ]);
    disp(potgrad{ridx});
  end

end


disp('-- Computing gradients of edge and loop gains.');

% We already have edge gains themselves from the loop analysis.
% Get the edge gain gradients.

edgegaingradients = synthRFH_getEdgeGainGradients( ...
  modelparams, intcouplings, oprates, rategrad );

% Add edge gain gradient information to the loop table.
% This computes loop gain gradients, storing in 'cyclegaingradient'.
loopinfo = synthRFH_addLoopGainInfo( loopinfo, edgegains, edgegaingradients );

edgesvalid = (abs(intcouplings) >= minweight);

save( [ outdir filesep 'edge-loop-gradients-' model_name '.mat' ], ...
  'edgegaingradients', 'edgesvalid', 'loopinfo', '-v7.3' );


if tattle_edge_gradients

  disp('-- Edge gains around exponential operating point:');

  disp(edgegains);

  disp('-- Edge gain gradients around exponential operating point:');

  for dstidx = 1:regioncount
    for srcidx = 1:regioncount
      if edgesvalid(dstidx, srcidx)

        disp([ '.. ' labels_from_indices{dstidx} ' <-- ' ...
          labels_from_indices{srcidx} ...
          sprintf( ' (gain %.3f):', edgegains(dstidx,srcidx) ) ]);
        disp(edgegaingradients{dstidx,srcidx});

      end
    end
  end

end


% Report loop gain gradients.

if tattle_loop_gradients

  disp('-- Loop gain gradients around exponential operating point:');

  for lidx = 1:length(loopinfo)
    thisrec = loopinfo(lidx);

    disp([ '.. Loop "' thisrec.label '" gain gradient:' ]);
    disp(thisrec.cyclegaingradient);
  end

end


disp('-- Finished gradient calculations.');



%
% Tune the model if requested.

if want_adjust_weights

  disp('-- Adjusting coupling weights.');
  tic;

  disp('.. Initial couplings:');
  disp(intcouplings);

  [ intcouplings opterr ] = synthRFH_optimizeCouplings( ...
    modelparams, intcouplings, loopinfo, adjust_goals, ...
    adjust_taulimit, adjust_bestfactor, adjust_maxprobes );

  disp('.. Final couplings:');
  disp(intcouplings);

  if (opterr > 0.1)
    disp(sprintf( '###  Optimization failed! Final error:  %.3f', opterr ));
  else
    % FIXME - Diagnostics.
    disp(sprintf( '.. Optimization succeeded. Final error:  %.3f', opterr ));
  end

  durstring = helper_makePrettyTime(toc);
  disp([ '-- Finished adjusting coupling weights (' durstring ').' ]);


  disp('-- Revised loop analysis:');

  loopinfo = synthRFH_findLoops( modelparams, intcouplings, minweight );
  edgegains = synthRFH_getEdgeGains( modelparams, intcouplings, oprates );
  loopinfo = synthRFH_addLoopGainInfo( loopinfo, edgegains, {} );

  [ looptext looptable ] = helper_makePrettyLoopTable(loopinfo);

  disp(looptext);

  disp('-- End of revised loop analysis.');
end



%
% This is the end of the file.
