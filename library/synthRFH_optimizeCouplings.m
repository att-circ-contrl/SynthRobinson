function [ bestcouplings besterr ] = synthRFH_optimizeCouplings( ...
  modelparams, startcouplings, loopinfo, loopgoals, ...
  taulimit, bestfactor, maxprobes )

% function [ bestcouplings besterr ] = synthRFH_optimizeCouplings( ...
%   modelparams, startcouplings, loopinfo, loopgoals, ...
%   taulimit, bestfactor, maxprobes )
%
% This attempts to optimize the internal coupling weights of a Robinson
% neural model to modify the loop gains to meet a specified set of
% constraints.
%
% This works via brute-force gradient descent (using "fsolve").
%
% Note that the output coupling matrix follows the constraints described
% in Robinson 2002, Freyer 2011, and Hindriks 2023: The nu_ix weights are
% set to the same values as the nu_ex weights, and the weights that are
% zero in the references are set to be zero here.
%
% "modelparams" is a model parameter structure with the fields described in
%   MODELPARAMSROBINSON.txt.
% "startcouplings" is a 4x4 matrix indexed by (destination,source) that
%   provides the coupling weights (in mV*s) between excitatory, inhibitory,
%   specific nucleus, and reticular nucleus neural populations.
% "loopinfo" is a loop metadata structure returned by synthRFH_findLoops().
% "loopgoals" is a Nx2 cell array. Each row contains a loop label in the
%   first column and a goal keyword in the second column. Goal keywords are
%   'biggest', 'grow', 'decay', or 'dontcare'. Loops with positive tau less
%   than the specified limit are growing, loops with negative tau with
%   absolute value less than the specified limit are decaying, and a growing
%   loop is "biggest" if its tau is smaller than all other growing tau
%   values by the specified factor. Loops not listed default to 'dontcare'.
% "taulimit" is the longest time constant (in seconds) that an envelope may
%   have to be counted as "growing" or "decaying". Typically 1.0 or less.
% "bestfactor" is the factor by which the "best" loop's tau must be smaller
%   than all other growing tau values. Typically 1.2-1.5.
% "maxprobes" is the maximum number of probes to make, or [] or NaN for the
%   default number (about 1400).
%
% "bestcouplings" is a perturbed version of "startcouplings" that produces
%   the desired loop behaviors.
% "besterr" is the error value associated with "bestcouplings". This is in
%   the range 0..1, with 0 being perfect and 1 being terrible.


% Initialize output.

bestcouplings = startcouplings;



%
% Augment the loop information structure with goals.

goallabels = {};
goalkeywords = {};

if ~isempty(loopgoals)
  goallabels = loopgoals(:,1);
  goalkeywords = loopgoals(:,2);
end

for lidx = 1:length(loopinfo)
  thislabel = loopinfo(lidx).label;
  thisgoal = 'dontcare';

  thisidx = min(find(strcmp( thislabel, goallabels )));
  if ~isempty(thisidx)
    thisgoal = goalkeywords{thisidx};
  end

  loopinfo(lidx).goal = thisgoal;
end

% FIXME - Diagnostics.
%disp(transpose( [ { loopinfo(:).label } ; { loopinfo(:).goal } ] ));



%
% Set up fsolve() and see what it gives us.

reducedparams = helper_getReducedFromCouplings(startcouplings);

anonfunc = @(testvec) helper_calcSolutionError( ...
  testvec, modelparams, loopinfo, taulimit, bestfactor );

if isempty(maxprobes)
  options = optimoptions('fsolve', 'Display', 'off');
elseif isnan(maxprobes)
  options = optimoptions('fsolve', 'Display', 'off');
else
  options = optimoptions('fsolve', 'Display', 'off', ...
    'MaxFunctionEvaluations', maxprobes);
end

[ reducedparams, besterr ] = fsolve( anonfunc, reducedparams, options );

bestcouplings = helper_getCouplingsFromReduced(reducedparams);


% Done.
end



%
% Helper Functions


function solutionerror = helper_calcSolutionError( ...
  reducedparams, modelparams, loopinfo, taulimit, bestfactor )

  %
  % Figure out the operating point, loop gains, and loop tau values.

  intcouplings = helper_getCouplingsFromReduced( reducedparams );

  % This can take non-negligible time.
  [ rates potentials ] = ...
    synthRFH_estimateOperatingPointExponential( ...
      modelparams, intcouplings, [] );

  edgegains = ...
    synthRFH_getEdgeGains( modelparams, intcouplings, rates );

  loopinfo = synthRFH_addLoopGainInfo( loopinfo, edgegains, {} );


  %
  % Get "growness", "decayness", and "bestness" for each loop.
  % These should be continuous values in the range 0 to 1.
  % Acceptance is at 0.2, clearly good is at 0.8 or so.

  looptau = [ loopinfo(:).envelopetau ];

  growmask = looptau > 0;
  decaymask = looptau < 0;


  % NOTE - To converge, "growness" and "decayness" both have to be defined
  % even when the loop is doing the opposite.
  % We're mapping correct-sense values to 0.2..1.0, and incorrect to 0..0.2.

  % This gives long-period at 0, accepted +1, good +5, short-period +inf.
  % Negative values are the wrong sense.

  growness = taulimit ./ looptau;
  decayness = -growness;

  % This maps (correct sense) values to 0.2 (0 growness), 0.4 (accepted),
  % 0.83 (growness +3), 1.0 (inf growness).

  correctgrowness = 1 ./ ( 1 + 4 * exp(-growness) );
  correctdecayness = 1 ./ ( 1 + 4 * exp(-decayness) );

  % This maps (incorrect sense) values to 0 (inf decayness), 0.1 (accepted
  % decayness), 0.2 (0 decayness).
  % It's linear rather than exponential so that we can still get a gradient
  % even with strong growth/decay in the direction we don't want.

  incorrectgrowness = 0.2 ./ (1 + growness);
  incorrectdecayness = 0.2 ./ (1 + decayness);

  % Multiplex these two together to get growness/decayness.

  growness = incorrectdecayness;
  growness(growmask) = correctgrowness(growmask);

  decayness = incorrectgrowness;
  decayness(decaymask) = correctdecayness(decaymask);


  % If we only have one growing loop, it's the best.
  % If we have at least two, find the second-best and compute FOMs.

  bestness = zeros(size(looptau));

  if 1 == sum(growmask)
    bestness(growmask) = 1;
  elseif any(growmask)
    % Find the second-shortest tau value.
    sortedvals = sort(looptau(growmask));
    secondbesttau = sortedvals(2);

    % 0 to inf, accepted at 1, goal at bestfactor (greater than 1).
    bestness = secondbesttau ./ looptau;

    % 0 to 1, accepted at 0.2, goal at 0.8.
    powerval = log(16) / log(bestfactor);
    bestness = 1 - 1 ./ (1 + 0.25 * (bestness .^ powerval));
    bestness(~growmask) = 0;
  end


  %
  % Derive a single figure of merit from the goals and tau values.


  % First, find per-loop goal-ness.
  % A "don't care" goal doesn't contribute to goal-ness.
  % We don't have to remember which loop each FOM came from.

  loopgoals = { loopinfo(:).goal };

  goalmaskgrow = strcmp('grow', loopgoals);
  goalmaskdecay = strcmp('decay', loopgoals);
  goalmaskbiggest = strcmp('biggest', loopgoals);

  % Have "biggest" imply "grow".
  goalmaskgrow = goalmaskgrow | goalmaskbiggest;

  % Concatenate these into one joint vector.
  % We have to do this rather than MUXing because "biggest" and "grow"
  % can make the same loop contribute twice.

  loopfoms = [ growness(goalmaskgrow), decayness(goalmaskdecay), ...
    bestness(goalmaskbiggest) ];


  % Now, turn the vector of "okay at 0.2, good at 0.8" values into an error.
  % We can just take the sum of (1-err), or the mean. This guarantees that
  % any given error's gradient contribution remains relevant.

  % The mean will have the range 0..1, with 0 being optimal.
  solutionerror = mean(1 - loopfoms);

end


function reducedparams = helper_getReducedFromCouplings( intcouplings )

  % FIXME - Magic matrix index values (per synthRFH_getRegionInfo).
  % FIXME - fsolve() wants a vector, not a struct, so magic vector indices.

  reducedparams = [];

  reducedparams(1) = intcouplings(1,1);    % ee
  reducedparams(2) = - intcouplings(1,2);  % ei - negative!
  reducedparams(3) = intcouplings(1,3);    % es

  % ie, ii, and is are equal to ee, ei, and es.

  reducedparams(4) = intcouplings(3,1);    % se
  reducedparams(5) = - intcouplings(3,4);  % sr - negative!

  reducedparams(6) = intcouplings(4,1);    % re
  reducedparams(7) = intcouplings(4,3);    % rs

  % All other weights are zero.


  % Map to -inf..+inf.
  for ridx = 1:length(reducedparams)
    reducedparams(ridx) = helper_getUncompressed( reducedparams(ridx) );
  end

end


function intcouplings = helper_getCouplingsFromReduced( reducedparams )

  % Map -inf..+inf to a biologically plausible range.
  for ridx = 1:length(reducedparams)
    reducedparams(ridx) = helper_getCompressed( reducedparams(ridx) );
  end


  % FIXME - Magic matrix index values (per synthRFH_getRegionInfo).
  % FIXME - fsolve() wants a vector, not a struct, so magic vector indices.

  intcouplings = zeros(4,4);

  intcouplings(1,1) = reducedparams(1);    % ee
  intcouplings(1,2) = - reducedparams(2);  % ei - negative!
  intcouplings(1,3) = reducedparams(3);    % es

  % ie, ii, and is are equal to ee, ei, and es.
  intcouplings(2,:) = intcouplings(1,:);

  intcouplings(3,1) = reducedparams(4);    % se
  intcouplings(3,4) = - reducedparams(5);  % sr - negative!

  intcouplings(4,1) = reducedparams(6);    % re
  intcouplings(4,3) = reducedparams(7);    % rs

  % All other weights are zero.

end


function comprange = helper_getCompressed( infrange )

  % This turns a -inf..+inf scalar into a biologically plausible coupling
  % weight (0.05..10.0, per Robinson 2002).

  % Map to 0..1.
  comprange = 1 / (1 + exp(-infrange));

  % Map _that_ to 0.05..10, on a log scale.
  % This is about -3..+2.3.
  comprange = -3 + (5.3 * comprange);
  comprange = exp(comprange);

end


function infrange = helper_getUncompressed( comprange )

  % This turns a biologically plausible coupling weight (0.05..10.0, per
  % Robinson 2002) back into a -inf..+inf scalar.

  % Convert to log scale. This is about -3..+2.3.
  comprange = log(comprange);

  % Map to 0..1 and clamp.
  comprange = (comprange + 3) / 5.3;
  comprange = max(0.000001,comprange);
  comprange = min(0.999999,comprange);

  % Map to -inf..inf.
  infrange = - log( (1 / comprange) - 1 );

end



%
% This is the end of the file.
