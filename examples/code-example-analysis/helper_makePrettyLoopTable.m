function [ tabletext tabledata ] = doOp_makePrettyLoopTable( loopinfo )

% function [ tabletext tabledata ] = doOp_makePrettyLoopTable( loopinfo )
%
% This makes a human-readable version of Robinson model loop metadata.
%
% "loopinfo" is a structure returned by synthRFH_findLoops(),
%   optionally with gain information added by synthRFH_addLoopGainInfo.
%
% "tabletext" is a character vector containing a human-readable table
%   describing the loops.
% "tabledata" is a Matlab table containing the same information.


tabletext = '';
tabledata = table();

loopcount = length(loopinfo);

if isfield(loopinfo, 'cyclegain')

  % Augmented information with cyclegain, envelopetau, cyclegradient.

  labels = { loopinfo.label };
  delays_ms = [ loopinfo.delay ] * 1000;
  cycgainsraw = [ loopinfo.cyclegainraw ];
  cycgains = [ loopinfo.cyclegain ];
  freqs = [ loopinfo.frequency ];
  envtaus_ms = [ loopinfo.envelopetau ] * 1000;

  % Transposing the structure array doesn't help, so do this manually.
  if ~iscolumn(labels) ; labels = transpose(labels); end
  if ~iscolumn(delays_ms) ; delays_ms = transpose(delays_ms); end
  if ~iscolumn(cycgainsraw) ; cycgainsraw = transpose(cycgainsraw); end
  if ~iscolumn(cycgains) ; cycgains = transpose(cycgains); end
  if ~iscolumn(freqs) ; freqs = transpose(freqs); end
  if ~iscolumn(envtaus_ms) ; envtaus_ms = transpose(envtaus_ms); end

  tabledata.('Label') = labels;
  tabledata.('Delay') = delays_ms;
  tabledata.('Gain') = cycgainsraw;
  tabledata.('Gain') = cycgains;
  tabledata.('Freq') = freqs;
  tabledata.('Tau') = envtaus_ms;

  tabletext = sprintf('  %-5s  %5s  %6s  %6s  %5s  %5s\n', ...
    'Label', 'Delay', 'Raw G', 'Gain', 'Freq', 'Tau' );
  tabletext = [ tabletext sprintf('  %-5s  %5s  %6s  %6s  %5s  %5s\n', ...
    '', '(ms)', '', '', '(Hz)', '(ms)' ) ];

  for lidx = 1:loopcount
    tabletext = [ tabletext ...
      sprintf('  %-5s  %5d  %6.2f  %6.2f  %5.1f  %5d\n', ...
        labels{lidx}, round(delays_ms(lidx)), cycgainsraw(lidx), ...
        cycgains(lidx), freqs(lidx), round(envtaus_ms(lidx)) ) ];
  end

else

  % Baseline information with label, regionsvisited, delay, isinverting,
  % frequency, and attenuation.

  labels = { loopinfo.label };
  delays_ms = [ loopinfo.delay ] * 1000;
  invflags = [ loopinfo.isinverting ];
  freqs = [ loopinfo.frequency ];
  attens = [ loopinfo.attenuation ];

  % Transposing the structure array doesn't help, so do this manually.
  if ~iscolumn(labels) ; labels = transpose(labels); end
  if ~iscolumn(delays_ms) ; delays_ms = transpose(delays_ms); end
  if ~iscolumn(invflags) ; invflags = transpose(invflags); end
  if ~iscolumn(freqs) ; freqs = transpose(freqs); end
  if ~iscolumn(attens) ; attens = transpose(attens); end

  tabledata.('Label') = labels;
  tabledata.('Delay') = delays_ms;
  tabledata.('Inv') = invflags;
  tabledata.('Freq') = freqs;
  tabledata.('Atten') = attens;

  invletters = cell(size(invflags));
  invletters(:) = { 'n' };
  invletters(invflags) = { 'Y' };

  tabletext = sprintf('  %-5s  %5s  %5s  %5s  %5s\n', ...
    'Label', 'Delay', 'Inv', 'Freq', 'Atten' );
  tabletext = [ tabletext sprintf('  %-5s  %5s  %5s  %5s  %5s\n', ...
    '', '(ms)', '', '(Hz)', '' ) ];

  for lidx = 1:loopcount
    tabletext = [ tabletext sprintf('  %-5s  %5d  %5s  %5.1f  %5.3f\n', ...
      labels{lidx}, round(delays_ms(lidx)), invletters{lidx}, ...
      freqs(lidx), attens(lidx) ) ];
  end

end


% Done.
end


%
% This is the end of the file.
