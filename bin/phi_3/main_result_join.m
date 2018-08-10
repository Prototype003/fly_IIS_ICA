%% Description

%{

For joining multiple result files

%}

%% Setup

nChannels = 4;
flies = (1:13);
conditions = (1:2);
%taus = 2.^(0:9);
taus = 2.^(0:6);

nchoosek(15, nChannels);
sets = floor(linspace(1, nchoosek(15, nChannels), 10));

trials = (1);
tau_offset = 0;

global_tpm = 1;

source_dir = 'results_split/';
source_prefix = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim'...
    '_nChannels' num2str(nChannels)...
    '_globalTPM' num2str(global_tpm)...
    ];

tau_string = 'tau'; % 'tauBin' or 'tau'

%% Load phis

% sets x flies x conditions x taus
phi_values = zeros(length(sets), length(flies), length(conditions), length(taus));

for set_counter = 1 : length(sets)
    network = sets(set_counter);
    disp(network);
    for fly = flies
        for condition = conditions
            for tau_counter = 1 : length(taus)
                for trial = trials
                    
                    % File name
                    source_file = [source_prefix...
                        '_f' sprintf('%02d', fly)...
                        'c' num2str(condition)...
                        tau_string num2str(taus(tau_counter))...
                        'tauOffset' num2str(tau_offset)...
                        's' sprintf('%04d', network)...
                        't' num2str(trial)...
                        '.mat'...
                        ];
                    
                    try
                        % Load file
                        tmp = load([source_dir source_file]);
                    catch ME
                        disp(ME.message);
                        disp(['Failed file: ' source_file]);
                        continue; % Skip to the next file, leave the data entry structure for this entry empty
                    end
                    
                    phi_values(set_counter, fly, condition, tau_counter) = tmp.phi.phi;
                    
                end
            end
        end
    end
end

%% Save (so we don't need to keep loading)

save(['tmp/' tau_string '.mat'], 'phi_values');