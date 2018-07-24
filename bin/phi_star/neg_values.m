%% Description

%{
Issue: some phi-star values for 2 channels are negative

Get networks/flies/conditions/trials which give negative phi_star
%}

%% Setup and Load

tau_counter = 1;

bin_dir = '../';

addpath(genpath([bin_dir 'PhiToolbox-master/']));

data_dir = [bin_dir 'workspace_results/'];
data_file = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim';

%% Get locations of negative values

values = phis{1}.phis;
channel_sets = nchoosek((1:15), 2);

neg_info = cell(0);
neg_counter = 0;
for fly = 1 : size(values, 3)
    for condition = 1 : size(values, 4)
        for network = 1 : size(values, 1)
            for trial = 1 : size(values, 2)
                if values(network, trial, fly, condition, tau_counter) < 0
                    neg_counter = neg_counter + 1;
                    tmp = struct();
                    tmp.fly = fly;
                    tmp.condition = condition;
                    tmp.network = network;
                    tmp.channel_set = channel_sets(network, :);
                    tmp.trial = trial;
                    neg_info{neg_counter} = tmp;
                end
            end
        end
    end
end

%% Try recomputing phi-star value to confirm

options.type_of_dist = 'discrete';
options.type_of_phi = 'star'; % 'star'
options.type_of_MIPsearch = 'Exhaustive';

% Parameters
params.tau = 4;
params.number_of_states = 2;

data = load([data_dir data_file]);
fly_data = data.fly_data;

for neg = 1 : length(neg_info)
    tmp = neg_info{neg};
    neg_info{neg}.data = fly_data(:, tmp.channel_set, tmp.trial, tmp.fly, tmp.condition)';
    medians = median(neg_info{neg}.data, 2);
    neg_info{neg}.medians = repmat(medians, [1, size(neg_info{neg}.data, 2)]);
    neg_info{neg}.data_binary = (neg_info{neg}.data > neg_info{neg}.medians) + 1;
    
    [mip_old, phi] = MIP_search(neg_info{neg}.data_binary, params, options);
    % mip_old format: each index corresponds to a channel, the element gives the partition which the channel belongs to, in the MIP
    
    % Reformat mip
    mip = partition_z2cell(mip_old);
    
    neg_info{neg}.phi = phi;
    neg_info{neg}.mip = mip_old;
    neg_info{neg}.mip_format2 = mip;
    neg_info{neg}.options = options;
    neg_info{neg}.params = params;

end
