function [moded, moded_frequencies] = mip_mode(mips, frequencies)
% Given a set of mips and their frequencies within each trial, finds the
% most frequently occurring state and its corresponding mip, as well as its
% frequency of occurrence within the trial
%
% Inputs:
%   mips = phithrees{n}.mips
%   frequencies = phithrees{n}.state_counters
%
% Outputs:
%   moded = mips, except holds the most frequent mip across states
%   moded_frequencies = frequency of the corresponding most frequent mip


moded_dimensions = size(frequencies);

moded = cell(moded_dimensions(2:end));
moded_frequencies = cell(moded_dimensions(2:end));

for tau = 1 : size(frequencies, 6)
    for condition = 1 : size(frequencies, 5)
        for fly = 1 : size(frequencies, 4)
            for set = 1 : size(frequencies, 2)
                for trial = 1 : size(frequencies, 3)
                    % Determine most frequent state - corresponding MIP is the most frequent MIP
                    [frequency, state] = max(frequencies(:, set, trial, fly, condition, tau));
                    moded{set, trial, fly, condition, tau} = mips{state, set, fly, condition, tau};
                    moded_frequencies{set, trial, fly, condition, tau} = frequency;
                end
            end
        end
    end
end

end