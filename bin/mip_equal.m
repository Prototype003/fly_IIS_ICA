function [equal] = mip_equal(a, b)
% Compares 2 MIPs for equality
%
% Inputs:
%   a = MIP; cell array, each cell holds a vector of nodes
%   b = MIP; cell array, each cell holds a vector of nodes
%
% Outputs:
%   equal = 1 if a==b, 0 otherwise

equal = 0;

% Compare each group of MIPA to each group of MIPB
% If there is a match, we can remove the matching group from the test partition so that we don't need to look at it again
% If the B becomes empty (all groups removed), and matches is the size of the MIPA, then all groups matched
match_count = 0;
for group_a = 1 : length(a)
    for group_b = 1 : length(b)
        if length(a{group_a}) == length(b{group_b}) % Check if groups are same size
            match_group = sum(sort(a{group_a}) == sort(b{group_b})); % Check if groups have the same nodes (after sorting)
            if match_group > 0
                b{group_b} = -1; % Mark group for deletion
                match_count = match_count + 1;
            end
        end
    end
end

% Delete matched groups in the test partition (delete in reverse to avoid out-of-index errors)
for group = length(b) : -1 : 1
    if b{group} == -1
        b(group) = [];
    end
end

if match_count == length(a) && isempty(b)
    equal = 1;
end

end