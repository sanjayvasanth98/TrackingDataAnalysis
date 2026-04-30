function [matDir, matFile] = test_plotmat_location(matFileName)
%TEST_PLOTMAT_LOCATION Central plot-data .mat folder for test scripts.
%
% Edit manualMatDir once here, then all Testing/test_*.m plotting scripts
% can use the same saved main-run .mat files.
%
% Usage:
%   matDir = test_plotmat_location();
%   [matDir, matFile] = test_plotmat_location("activation_summary_by_case.mat");

if nargin < 1
    matFileName = "";
end
matFileName = string(matFileName);

% ---- EDIT THIS PATH WHEN YOU WANT ALL TEST PLOTS TO USE A NEW RUN -------
manualMatDir = "C:\Users\kbsanjayvasanth\Downloads\plot_data_mat (1)";

candidateDirs = strings(0, 1);
candidateDirs(end+1, 1) = manualMatDir;
candidateDirs(end+1, 1) = "C:\Users\kbsanjayvasanth\Downloads\plot_data_mat";

userProfile = string(getenv('USERPROFILE'));
if strlength(userProfile) > 0
    candidateDirs(end+1, 1) = string(fullfile(userProfile, 'Downloads', 'plot_data_mat (1)'));
    candidateDirs(end+1, 1) = string(fullfile(userProfile, 'Downloads', 'plot_data_mat'));
end

homeDir = string(getenv('HOME'));
if strlength(homeDir) > 0
    candidateDirs(end+1, 1) = string(fullfile(homeDir, 'Downloads', 'plot_data_mat (1)'));
    candidateDirs(end+1, 1) = string(fullfile(homeDir, 'Downloads', 'plot_data_mat'));
end

candidateDirs = unique(candidateDirs(strlength(strtrim(candidateDirs)) > 0), 'stable');

matDir = "";
for i = 1:numel(candidateDirs)
    if isfolder(candidateDirs(i))
        if strlength(strtrim(matFileName)) == 0 || isfile(fullfile(candidateDirs(i), matFileName))
            matDir = candidateDirs(i);
            break;
        end
    end
end

if strlength(matDir) == 0
    msg = sprintf('Could not find a plot_data_mat folder containing "%s".\nChecked:\n', char(matFileName));
    for i = 1:numel(candidateDirs)
        msg = sprintf('%s  %s\n', msg, char(candidateDirs(i)));
    end
    msg = sprintf('%s\nEdit manualMatDir in Testing/test_plotmat_location.m.', msg);
    error('%s', msg);
end

if strlength(strtrim(matFileName)) > 0
    matFile = string(fullfile(matDir, matFileName));
else
    matFile = "";
end
end
