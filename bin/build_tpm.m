function [ tpm ] = build_tpm(fly_data, tau, n_values)
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
        
        % Identify current state
        state_current = state2loli_index(sample_current);
        
        % Identify future state
        state_future = state2loli_index(sample_future);
        
        % Increment TPM transition by 1
        tpm(state_current, state_future) = tpm(state_current, state_future) + 1;
        
        % Increment transition counter
        transition_counter(state_current) = transition_counter(state_current) + 1;
    end
end

% Divide elements in TPM by transition counter
% If counter is 0, then transition never occurred - to avoid dividing by 0
for counter = 1 : length(transition_counter)
    if transition_counter(counter) == 0
        transition_counter(counter) = 1;
    end
end
for future_state = 1 : size(tpm, 2)
    tpm(:, future_state) = tpm(:, future_state) ./ transition_counter;
end

    function [index] = state2loli_index(state)
        % Placeholder function to mimic pyphi.convert.state2loli_index() for the 2 channel scenario
        % Inputs:
        %   state = vector of 1s and 0s
        % Outputs:
        %   index = corresponding loli index (1-indexed, not 0-indexed as in Python)
        
        if length(state) == 1
            index = state + 1; % Assumes that states start from 0, not 1 (e.g. 0,1 and not 1,2)
            return
        end
        
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

