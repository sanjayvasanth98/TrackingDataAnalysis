function filtered = filter_breakup_by_ar(allBreakup, arThreshold)
%FILTER_BREAKUP_BY_AR  Keep only breakup events where parentAR >= arThreshold.
%
%   filtered = filter_breakup_by_ar(allBreakup, arThreshold)
%
%   allBreakup:   struct array with fields .caseName, optional .Re, .kD, .events
%   arThreshold:  scalar, minimum parent aspect ratio to retain
%
%   Returns the same struct array with events filtered per case.

filtered = allBreakup;
for ci = 1:numel(filtered)
    ev = filtered(ci).events;
    if isempty(ev), continue; end
    keep = [ev.parentAR] >= arThreshold;
    filtered(ci).events = ev(keep);
end
end
