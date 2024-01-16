function newloopinfo = synthRFH_addLoopGainInfo( ...
  oldloopinfo, edgegains, edgegaingradients )

% function newloopinfo = synthRFH_addLoopGainInfo( ...
%   oldloopinfo, edgegains, edgegaingradients )
%
% This accepts a loop metadata structure array and augments each record
% with gain-related information.
%
% An empty cell array may be supplied for the gain gradients, to omit
% gradient information.
%
% "oldloopinfo" is a structure returned by synthRFH_findLoops().
% "edgegains" is a 4x4 matrix indexed by (destination, source) that contains
%   the small-signal firing rate gains between each source and destination
%   for the excitatory, inhibitory, specific nucleus, and reticular nucleus
%   neural populations.
% "edgegaingradients" is a 4x4 cell array indexed by (destination, source)
%   that contains the gradient with respect to "intcouplings" of the
%   small-signal firing rate gains in "edgegains". Passing an empty cell
%   array skips calculation of loop gain gradients.
%
% "newloopinfo" is a copy of "oldloopinfo" with the following fields added
%   to each record:
%   "cyclegainraw" is the small-signal gain from traversing once around the
%     loop, without taking into account filter attenuation.
%   "cyclegain" is the small-signal gain from traversing once around the
%     loop with filter attenuation taken into account.
%   "envelopetau" is the time constant for the growth (positive) or decay
%     (negative) of the oscillation envelope. The envelope is exp(t/tau).
%   "cyclegaingradient" is the gradient with respect to "intcouplings" of
%     "cyclegain". This is a 4x4 matrix (per "intcouplings"). If an empty
%     cell array is passed as "edgegaingradients", this field is omitted.


wantgrad = ~isempty(edgegaingradients);

regioncount = size(edgegains,1);
loopcount = length(oldloopinfo);


newloopinfo = oldloopinfo;

for lidx = 1:loopcount
  thispath = oldloopinfo(lidx).regionsvisited;
  thispathlength = length(thispath) - 1;

  thisdelay = oldloopinfo(lidx).delay;
  thisinvert = oldloopinfo(lidx).isinverting;

  thiscyclegain = 1;

  for pidx = 1:thispathlength
    srcidx = thispath(pidx);
    dstidx = thispath(pidx+1);

    thiscyclegain = thiscyclegain * edgegains(dstidx, srcidx);
  end

  newloopinfo(lidx).cyclegainraw = thiscyclegain;

  thiscyclegain = thiscyclegain * oldloopinfo(lidx).attenuation;
  thistau = thisdelay / log(abs(thiscyclegain));

  newloopinfo(lidx).cyclegain = thiscyclegain;
  newloopinfo(lidx).envelopetau = thistau;

  if wantgrad

    % Compute gradients using the product rule.

    thiscyclegrad = zeros(regioncount, regioncount);

    for pidx1 = 1:thispathlength
      thisterm = 1;

      for pidx2 = 1:thispathlength
        srcidx = thispath(pidx2);
        dstidx = thispath(pidx2+1);

        if pidx1 == pidx2
          thisterm = thisterm * edgegaingradients{dstidx,srcidx};
        else
          thisterm = thisterm * edgegains(dstidx,srcidx);
        end
      end

      thiscyclegrad = thiscyclegrad + thisterm;
    end

    % This gets attenuated just like gain did.
    thiscyclegrad = thiscyclegrad * oldloopinfo(lidx).attenuation;

    newloopinfo(lidx).cyclegaingradient = thiscyclegrad;

  end
end



% Done.
end


%
% This is the end of the file.
