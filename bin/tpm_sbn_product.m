function [tpm_prod] = tpm_sbn_product(tpm_ind)
% Combines two independent state-by-node TPMs (using LOLI/fortran indexing)
%
% Inputs:
%   tpm_ind = independent TPMs for single binary nodes; matrix with
%       dimensions (states x nodes) = (2 x nodes); each column gives the
%       state-by-node TPM for a single node
%
% Outputs:
%   tpm_prod = combined state-by-node TPM; matrix with dimensions 
%       (2^nodes x number of nodes)

n_nodes = size(tpm_ind, 2);

n_states = 2^n_nodes;

tpm_prod = zeros(n_states, n_nodes);

% Build list of possible states
% Iterate through each state and convert decimal index into binary string
state_list = zeros(size(tpm_prod));
for state = 1 : n_states
    state_list(state, :) = dec2bi(state-1, n_nodes);
end

% Iterate through each node, copy probabilities from independent TPMs
for node = 1 : n_nodes
    tpm_prod(:, node) = tpm_ind(state_list(:, node)+1, node);
end

    function [bit_string] = dec2bi(dec, nBits)
        % Converts decimal number into binary string with LOLI/fortran
        % order; lowest order bit is on the left
        %
        % Inputs:
        %   dec = decimal number; within range which can be represented by
        %       nBits bits
        %   nBits = number of bits to consider
        %
        % Outputs:
        %   bin = binary string of nBits bits
        
        bit_string = zeros(1, nBits);
        for bit = 1 : nBits
            fits = floor(dec / (2^(nBits-1))) > 0;
            if fits == 1
                bit_string(bit) = 1;
                dec = dec - 2^(nBits-1);
            end
            nBits = nBits - 1;
        end
        
        bit_string = fliplr(bit_string);
    end

end

