function [ output_args ] = build_tpm(fly_data, tau, n_values)
% NOTE: Requires implementation of pyphi.convert.state2loli_index() from pyphi
%
% Builds a tpm for one fly and one condition
% Based off of build_tpm in fly_phi.py
%
% Inputs:
%   fly_data = matix (of discretised data) with dimensions (samples x channels x epoch-trials)
%       Holds data for one fly, one condition
%   tau = integer - the lag between current and future states
%   n_values = number of states each *node* can be in (e.g. 2 for ON and OFF)
%
% Outputs:
%   tpm = matrix with dimensions (n_values^channels x n_values^channels)
%       Each row holds the probabilities of a past state transitioning into future states (columns)

% Determine number of possible system states
n_states = n_values ^ size(fly_data, 2);

% Declare TPM
tpm = zeros(n_states, n_states);

transition_counter = zeros(n_states, 1); % transition counter for each state

for trial = 1 : size(fly_data, 3)
    for sample = 1 : size(fly_data, 1) - tau
        sample_current = fly_data(sample, :, trial);
        sample_future = fly_data(sample+tau, :, trial);
    end
end

    function [index] = state2loli_index(state)
        % Placeholder function to mimic pyphi.convert.state2loli_index() for the 2 channel scenario
        % Inputs:
        %   state = vector of 1s and 0s
        % Outputs:
        %   index = corresponding loli index (1-indexed, not 0-indexed as in Python)
        
        if state == [0 0]
            index = 1;
        elseif state == [1 0]
            index = 2;
        elseif state == [0 1]
            index = 3;
        else % state == [1 1]
            index = 4;
        end
    end

end

