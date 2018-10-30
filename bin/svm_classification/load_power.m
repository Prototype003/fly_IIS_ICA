function [powers, frequencies, chronux_params] = load_power(bin_location)
% Loads power values
% Note power values are averaged across trials if across flies is specified
%
% Inputs:
%   across: 0 = within flies, 1 = across flies
%       Not used, now just loads power for within flies, across flies is obtained
%       by averaging across trials (the second dimension in the powers matrix
%
% Outputs:
%

results_directory = [bin_location 'workspace_results/'];
%results_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_power_classification_across' num2str(0) '.mat'];
results_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_power_classification.mat'];

load([results_directory results_filename]);


end

