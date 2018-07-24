function [ mip_new ] = partition_z2cell(mip_old)
% Converts Phi-toolbox style Z partition format into cell array
%
% Inputs:
%   mip_old:
%       partition in vector format:
%           Each index corresponds to a node
%           Each value corresponds to the group to which the node is placed
%
% Outputs:
%   mip_new:
%       partition in cell format:
%           Each cell corresponds to a group, and holds a vector holding all
%           the nodes which are part of that group

mip_new = cell(max(mip_old), 1);

for group = 1 : max(mip_old)
    mip_new{group} = mip_old(mip_old == group);
end

end

