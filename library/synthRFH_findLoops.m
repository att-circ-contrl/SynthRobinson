function loopinfo = ...
  synthRFH_findLoops( modelparams, intcouplings, minweight )

% function loopinfo = ...
%   synthRFH_findLoops( modelparams, intcouplings, minweight )
%
% This examines a Robinson model coupling matrix and identifies loops.
% Loop metadata is extracted.
%
% This does not extract loop gain, since that varies with operating point.
%
% "modelparams" is a model parameter structure with the fields described in
%   MODELPARAMSROBINSON.txt.
% "intcouplings" is a 4x4 matrix indexed by (destination, source) that
%   provides the coupling weights (in mV*s) between excitatory, inhibitory,
%   specific nucleus, and reticular nucleus neural populations.
% "minweight" is the threshold to use when evaluating whether a coupling
%   weight is nonzero (the absolute value must be at least "minweight").
%
% "loopinfo" is a structure array with one element per identified loop, and
%   the following fields (per LOOPINFO.txt):
%   "label" is a unique identifier for the loop (the concatenation of the
%     letters associated with each region in the loop's path).
%   "regionsvisited" is a vector containing the indices of each region in the
%     loop's path. The first index is repeated as the last index.
%   "delay" is the propagation time for one cycle around the loop, in seconds.
%   "isinverting" is true if the product of the loop's edge couplings is
%     negative, and false if the product is positive.
%   "frequency" is the loop's fundamental mode frequency in Hz (1/delay if
%     non-inverting, half that if inverting).
%   "attenuation" is the loop's attenuation at its fundamental mode
%     frequency from the alpha, beta, and gamma model parameters. This will
%     be 1 if the signal is passed perfectly and less than 1 if not.


%
% Get metadata.

% Low-frequency approximation of  neuron dyamics delays.
% Gamma has two poles at the same frequency, so it gets applied twice.
delay_albet = (1 / modelparams.alpha) + (1 / modelparams.beta);
delay_gamma = (1 / modelparams.gamma) * 2;

% Convert the cortex/thalamus delay to seconds.
delay_ct = modelparams.halfdelay_ms * 0.001;

% Convenience copies.
alpha = modelparams.alpha;
beta = modelparams.beta;
gamma = modelparams.gamma;


% Be flexible about the number of regions, in case we extend the model.

regioncount = size(intcouplings,1);

[ indices_lut names_lut ] = synthRFH_getRegionInfo();

regionletters = { names_lut.letter };

% Identify cortex and non-cortex regions by black magic.
scratch = { names_lut.lutfield };
iscortex = contains( scratch, 'cortex' );

% Likewise with regions where the gamma delay applies (excitatory cortex).
isgamma = contains( scratch, 'excitatory' ) & iscortex;

% Build an "edge exists" mask.
edgesvalid = abs(intcouplings) >= minweight;



%
% Search for loops.

loopinfo = struct( 'label', {}, 'regionsvisited', {}, 'delay', {}', ...
  'isinverting', {}, 'frequency', {}, 'attenuation', {} );

% Do this by brute force permutation searching.
% The number of possible loops is small enough for this to be practical.

for looplength = 1:regioncount

  regionsleft = 1:regioncount;

  pathlist = helper_getPathFragments( ...
    regionsleft, looplength, edgesvalid, true);

  for pidx = 1:length(pathlist)
    % Only proceed if the final edge (endpoint to starting point) is valid,
    % since the helper only considers open-path edges.

    thispath = pathlist{pidx};

    dstidx = thispath(1);
    srcidx = thispath(length(thispath));

    if edgesvalid(dstidx,srcidx)

      % Get the closed-loop path (appending a copy of the starting point).
      thispath = [ thispath thispath(1) ];

      % Build path metadata.

      thislabel = '';
      thisdelay = 0;
      thisproduct = 1;

      hadgamma = false;

      for nidx = 2:length(thispath)
        srcidx = thispath(nidx-1);
        dstidx = thispath(nidx);

        thisweight = intcouplings(dstidx,srcidx);

        thislabel = [ thislabel regionletters{srcidx} ];
        thisproduct = thisproduct * thisweight;

        thisdelay = thisdelay + delay_albet;
        if isgamma(srcidx)
          thisdelay = thisdelay + delay_gamma;
          hadgamma = true;
        end
        if iscortex(srcidx) ~= iscortex(dstidx)
          thisdelay = thisdelay + delay_ct;
        end
      end

      thisinverting = (thisproduct < 0);
      thisfreq = 1 / thisdelay;
      if thisinverting
        thisfreq = 0.5 * thisfreq;
      end


      % Calculate attenuation due to filter effects.
      % It's (a b) / (jw + a) (jw + b) for all nodes, and g^2 / (jw + g)^2
      % if we have gamma.

      thisjomega = i * thisfreq * 2 * pi;
      thiscount = length(thispath) - 1;

      thisatten = (alpha * beta) ...
        / ( ( thisjomega + alpha ) * (thisjomega + beta) );
      thisatten = thisatten ^ thiscount;

      if hadgamma
        scratch = gamma / ( thisjomega + gamma );
        thisatten = thisatten * scratch * scratch;
      end

      thisatten = abs(thisatten);


      thisrec = struct( 'label', thislabel, 'regionsvisited', thispath, ...
        'delay', thisdelay, 'isinverting', thisinverting, ...
        'frequency', thisfreq, 'attenuation', thisatten );

      loopinfo = [ loopinfo thisrec ];

    end
  end

end


% Done.
end


%
% Helper Functions

function pathlist = helper_getPathFragments( ...
  nodes, pathlength, edgesvalid, extcall )

  pathlist = {};

  if (pathlength > 0) && (length(nodes) >= pathlength)
    for nidx = 1:length(nodes)

      thisprefix = nodes(nidx);

      if pathlength > 1

        % If this is an external call, prune everything before the current
        % node as well as the current node.
        % That way we don't get every possible starting point for each loop.
        % If this is a recursive call, allow all permutations.

        if extcall
          % If this is shorter than the path length, recursion will return
          % an empty fragment list, which is fine.
          newnodes = nodes(nidx+1:length(nodes));
        else
          newnodes = nodes(nodes ~= thisprefix);
        end

        suffixlist = helper_getPathFragments( ...
          newnodes, pathlength - 1, edgesvalid, false );

        for sidx = 1:length(suffixlist)
          thissuffix = suffixlist{sidx};

          % Only append the prefix and suffix if the connection is valid.
          thissrc = thisprefix;
          thisdst = thissuffix(1);

          if edgesvalid(thisdst,thissrc)
            pathlist = [ pathlist { [ thisprefix thissuffix ] } ];
          end
        end

      else
        pathlist = [ pathlist { thisprefix } ];
      end

    end
  end

end


%
% This is the end of the file.
