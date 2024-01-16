function [ indices_lut names_lut ]  = synthRFH_getRegionInfo()

% function [ indices_lut names_lut ]  = synthRFH_getRegionInfo()
%
% This returns metadata associating each region simulated by the Robinson
% 2002 model with a row/column index, a pretty name, and abbreviated names.
%
% No arguments.
%
% "indices_lut" is a structure with the following fields:
%   "cortex_excitatory" is the row/column index corresponding to excitatory
%     neurons in the cortex.
%   "cortex_inhibitory" is the row/column index corresponding to inhibitory
%     neurons in the cortex.
%   "thalamus_specific" is the row/column index corresponding to "specific
%     nucleus" neurons in the thalamus (also called the relay population).
%   "thalamus_reticular" is the row/column index corresponding to "reticular
%     nucleus" neurons in the thalamus.
%
% "names_lut" is a structure array indexed by region number, with the
%   following fields (all character vectors):
%   "title" is a verbose plot-safe name.
%   "label" is a terse filename-safe and plot-safe name.
%   "letter" is a capital letter.
%   "lutfield" is the name of the corresponding field in "indices_lut".


% Follow the naming conventions in Hindriks 2023 and Robinson 2002.

regiontitles = { 'Excitatory Cortex', 'Inhibitory Cortex', ...
  'Specific Nucleus', 'Reticular Nucleus' };

regionlabels = { 'excitatory', 'inhibitory', 'specific', 'reticular' };

regionletters = { 'E', 'I', 'S', 'R' };

regionfields = { 'cortex_excitatory', 'cortex_inhibitory', ...
  'thalamus_specific', 'thalamus_reticular' };


% Build the lookup tables.

names_lut = struct( 'title', regiontitles, 'label', regionlabels, ...
  'letter', regionletters, 'lutfield', regionfields );

indices_lut = struct();
for lidx = 1:length(names_lut)
  indices_lut.( names_lut(lidx).lutfield ) = lidx;
end


% Done.
end


%
% This is the end of the file.
