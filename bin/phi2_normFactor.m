function [normal_factor] = phi2_normFactor(partition, cov_past_past)
% Calculates normalisation factor for phi for a given partition
%
% Taken from phi_comp.m from Haun's phi_toolbox_Feb2014
%
% Inputs:
%   partition
%   cov_past_past

nParts = max(partition);

entropies_past = zeros(nParts, 1);

%grab the different quantities for the parts, as desribed by the partition
M_cell = cell(nParts, 1);
for i = 1: nParts
    M_cell{i} = find(partition == i);
end

for part_counter = 1 : nParts
    M = M_cell{part_counter};
    CovXtau_p = cov_past_past(M, M);
    entropies_past(part_counter) = H_gauss(CovXtau_p);
end

normal_factor = (nParts - 1) * min(entropies_past);

end

