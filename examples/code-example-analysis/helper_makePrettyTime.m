function durstring = helper_makePrettyTime(dursecs)

% function durstring = helper_makePrettyTime(dursecs)
%
% This formats a duration (in seconds) in a meaningful human-readable way.
% Examples would be "5.0 ms" or "5d12h".
%
% "dursecs" is a duration in seconds to format. This may be fractional.
%
% "durstring" is a character array containing a terse human-readable
%   summary of the duration.


% Initialize.
durstring = '-bogus-';


% Check for durations of seconds or less.
if dursecs < 1e-4
  durstring = sprintf('%.1e s', dursecs);
elseif dursecs < 2e-3
  durstring = sprintf('%.2f ms', dursecs * 1e3);
elseif dursecs < 2e-2
  durstring = sprintf('%.1f ms', dursecs * 1e3);
elseif dursecs < 2e-1
  durstring = sprintf('%d ms', round(dursecs * 1e3));
elseif dursecs < 2
  durstring = sprintf('%.2f s', dursecs);
elseif dursecs < 20
  durstring = sprintf('%.1f s', dursecs);
else
  % We're in days/hours/minutes/seconds territory.

  dursecs = round(dursecs);
  scratch = dursecs;
  dursecs = mod(scratch, 60);
  scratch = round((scratch - dursecs) / 60);
  durmins = mod(scratch, 60);
  scratch = round((scratch - durmins) / 60);
  durhours = mod(scratch, 24);
  durdays = round((scratch - durhours) / 24);

  if durdays > 0
    durstring = ...
      sprintf('%dd%dh%dm%ds', durdays, durhours, durmins, dursecs);
  elseif durhours > 0
    durstring = sprintf('%dh%dm%ds', durhours, durmins, dursecs);
  elseif durmins > 0
    durstring = sprintf('%dm%ds', durmins, dursecs);
  else
    durstring = sprintf('%d s', dursecs);
  end
end


% Done.

end


%
% This is the end of the file.
