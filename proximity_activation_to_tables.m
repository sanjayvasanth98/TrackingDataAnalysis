function [pairTable, eventTable, binnedStats, summaryTable] = proximity_activation_to_tables(allProximityActivation)
%PROXIMITY_ACTIVATION_TO_TABLES  Concatenate proximity-analysis outputs.

pairTable = table();
eventTable = table();
binnedStats = table();
summaryTable = table();

if isempty(allProximityActivation)
    return;
end

for i = 1:numel(allProximityActivation)
    d = allProximityActivation(i);
    if isfield(d, 'pairTable') && istable(d.pairTable) && height(d.pairTable) > 0
        pairTable = append_table_compat_local(pairTable, d.pairTable);
    end
    if isfield(d, 'eventTable') && istable(d.eventTable) && height(d.eventTable) > 0
        eventTable = append_table_compat_local(eventTable, d.eventTable);
    end
    if isfield(d, 'binnedStats') && istable(d.binnedStats) && height(d.binnedStats) > 0
        binnedStats = append_table_compat_local(binnedStats, d.binnedStats);
    end
    if isfield(d, 'summaryTable') && istable(d.summaryTable) && height(d.summaryTable) > 0
        summaryTable = append_table_compat_local(summaryTable, d.summaryTable);
    end
end
end


% =========================================================================
function outTbl = append_table_compat_local(outTbl, newRows)
if isempty(newRows) || height(newRows) == 0
    return;
end
if isempty(outTbl) || width(outTbl) == 0
    outTbl = newRows;
else
    outTbl = [outTbl; newRows]; %#ok<AGROW>
end
end
