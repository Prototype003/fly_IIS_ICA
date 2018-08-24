%% Description

%{

Classification across trials for a single channel set, single fly

Need to classify:
    Using 2nd order concepts only
    Using 3rd order concepts only
    Using 4th order concepts only
For
    State dependent (all states x concepts)
    Unweighted mean (sum of concept phis across possible states / number of possible states)
    Weighted mean (weighted mean of concept phis across states)

%}

%% Setup

addpath('C:\Users\this_\Documents\MATLAB\Toolboxes\liblinear-2.20\windows');

addpath('../svm_classification/');

nChannels = 4;
fly = 1;
conditions = (1:2);
channel_set = 2;
tau = 4;
trials = (1:8);
global_tpm = 0;
tau_bin = 0;
sample_offset = 0;

% split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels4_globalTPM0_f03c2tau4tauOffset0s0939t6.mat
file_prefix = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels' num2str(nChannels)];

nStates = 2^nChannels;
nConcepts = 15;

% (state-dependence)x(which concepts to use)
classification_features = cell(3, 5);
for concept_column = 1 : size(classification_features, 2)
    classification_features{1, concept_column}.states = (1:nStates); classification_features{1, concept_column}.average_states = 0; % All states
    classification_features{2, concept_column}.states = (1:nStates); classification_features{2, concept_column}.average_states = 1; % Unweighted averaged states
    classification_features{3, concept_column}.states = nStates+1; classification_features{3, concept_column}.average_states = 0; % Weighted averaged states
end
for state_row = 1 : size(classification_features, 1)
    classification_features{state_row, 1}.concepts = (1:4); % first order
    classification_features{state_row, 2}.concepts = (5:10); % second order
    classification_features{state_row, 3}.concepts = (11:14); % third order
    classification_features{state_row, 4}.concepts = (15); % fourth order
    classification_features{state_row, 5}.concepts = (5:15); % 2nd-4th orders
end

%% Load compositions for each trial

tables = zeros(nStates + 1, 15, length(conditions), length(trials));

phi_states = zeros(nStates, length(conditions), length(trials));
phi_unweighted = zeros(length(conditions), length(trials));
phi_weighted = zeros(length(conditions), length(trials));

for condition = conditions
    for trial = trials
        
        file_infix = [...
            '_globalTPM' num2str(global_tpm)...
            '_f' sprintf('%02d', fly)...
            'c' num2str(condition)...
            'tau' num2str(tau)...
            'tauOffset' num2str(sample_offset)...
            's' sprintf('%04d', channel_set)...
            't' num2str(trial)...
            '.mat'
            ];
        
        [phi, diff_table, unparted, parted] = composition_table(nChannels, [file_prefix file_infix]);
        
        tables(:, :, condition, trial) = diff_table; % Small phi values
        
        phi_weighted(condition, trial) = phi.phi; % Big phi values
        phi_unweighted(condition, trial) = mean(phi.state_phis);
        phi_states(:, condition, trial) = phi.state_phis;
        
    end
end

%% Classify using each set of features

classifications = cell(size(classification_features));

cost_powers = (-20 : 20);
costs = 2.^cost_powers;

for state_features = 1 : size(classification_features, 1)
    for concept_features = 1 : size(classification_features, 2)
        
        % Get relevant features from full composition tables
        feature_table = tables(classification_features{state_features, concept_features}.states, classification_features{state_features, concept_features}.concepts, :, :);
        
        % Calculate average across states if required
        if classification_features{state_features, concept_features}.average_states == 1
            feature_table = mean(feature_table, 1);
        end
        
        % Collapse features to the first dimension as required
        feature_table_collapsed = zeros(size(feature_table, 1) * size(feature_table, 2), size(feature_table, 3), size(feature_table, 4));
        for row = 1 : size(feature_table, 1)
            feature_table_collapsed((1:size(feature_table, 2))+((row-1)*size(feature_table, 2)), :, :) = permute(feature_table(row, :, :, :), [2 3 4 1]);
        end
        
        % Reformat for input to SVM classification
        feature_table = permute(feature_table_collapsed, [3 1 2]); % trials x features x conditions (observations x features x classes)
        
        cost_results = cell(size(costs));
        for cost_counter = 1 : length(costs)
            % Submit for classification
            cost_results{cost_counter} = svm_lol_liblinear_manual(feature_table, costs(cost_counter));
            classifications{state_features, concept_features} = cost_results;
        end
        
    end
end

% Classification using just big-phi values (1 feature)
phi_classifications = cell(3, length(costs));
for cost_counter = 1 : length(costs)
    phi_classifications{1, cost_counter} = svm_lol_liblinear_manual(permute(phi_states, [3 1 2]), costs(cost_counter));
    phi_classifications{2, cost_counter} = svm_lol_liblinear_manual(permute(phi_unweighted, [2 3 1]), costs(cost_counter));
    phi_classifications{3, cost_counter} = svm_lol_liblinear_manual(permute(phi_weighted, [2 3 1]), costs(cost_counter));
end

%% Extract classification accuracy for each feature-type

accuracies = zeros([size(classification_features) length(costs)]);

for state_features = 1 : size(classification_features, 1)
    for concept_features = 1 : size(classification_features, 2)
        for cost = 1 : length(costs)
            accuracies(state_features, concept_features, cost) = classifications{state_features, concept_features}{cost}.accuracy;
        end
    end
end

phi_accuracies = zeros(3, length(costs));
for state_features = 1 : size(phi_classifications, 1)
    for cost = 1 : length(costs)
        phi_accuracies(state_features, cost) = phi_classifications{state_features, cost}.accuracy;
    end
end

% Join into one matrix
big_phi_mat = zeros(3, 1, length(costs));
big_phi_mat(:, :, :) = NaN;
big_phi_mat(:, 1, :) = phi_accuracies; % This 3rd row corresponds to weighted average across states
accuracies = cat(2, accuracies, big_phi_mat);

%% Plot accuracies

% Plot format:
%   1 row per feature type
%   1 colour/line pattern per state dependence type

figure;

dependence_colours = 'rgb';

feature_titles = {'1st', '2nd', '3rd', '4th', '2/3/4th', '\Phi'};

for feature_type = 1 : size(classification_features, 2)
    subplot(size(classification_features, 2), 1, feature_type);
    for state_dependence = 1 : size(classification_features, 1)
        line = permute(accuracies(state_dependence, feature_type, :), [3 1 2]);
        plot(cost_powers, line, dependence_colours(state_dependence)); hold on;
        title(feature_titles{feature_type});
        ylim([0.5 1]);
    end
    
    if feature_type == 1
        legend('x16 states' ,'unweighted mean', 'weighted mean');
    end
end

xlabel('cost (2^n)');
ylabel('accuracy');

%% Video frames for cost search
figure;
clim = [min(accuracies(:)) max(accuracies(:))];

frames = cell(size(costs));
for cost = 1 : length(costs)
    imagesc(accuracies(:, :, cost), clim); c=colorbar; ylabel(c, '% accuracy');
    title(['Cost = 2^{' num2str(cost_powers(cost)) '}']);
    set(gca, 'XTick', (1:size(accuracies, 2)), 'XTickLabel', {'1st', '2nd', '3rd', '4th', '2/3/4th', '\Phi'}); xlabel('concept order');
    set(gca, 'YTick', (1:size(accuracies, 1)), 'YTickLabel', {'x16 states', 'unweighted avg', 'weighted avg'}); ylabel('state dependence');
    drawnow;
    frame = getframe(gcf);
    frames{cost} = frame2im(frame);
end

%% Write frames into gif

video_duration = 20; % in seconds
frame_duration = video_duration / length(costs);

output_file = 'cost_search.gif';

for frame = 1 : length(frames)
    [mapped_frame, map] = rgb2ind(frames{frame}, 256);
    if frame == 1
        imwrite(mapped_frame, map, output_file, 'gif', 'LoopCount', Inf, 'DelayTime', frame_duration);
    else
        imwrite(mapped_frame, map, output_file, 'gif', 'WriteMode', 'append', 'DelayTime', frame_duration);
    end
end

%% Write frames into video

 % create the video writer with 1 fps
 writerObj = VideoWriter('cost_search.avi');
 writerObj.FrameRate = 1 / (video_duration / length(frames));

 % open the video writer
 open(writerObj);

 % write the frames to the video
 for u=1:length(frames)
     % convert the image to a frame
     frame = im2frame(frames{u});
     
     % write to video
     writeVideo(writerObj, frame);
 end

 % close the writer object
 close(writerObj);


%% Extract classification accuracy for each feature-type

% accuracies = zeros(size(classification_features));
% 
% for state_features = 1 : size(classification_features, 1)
%     for concept_features = 1 : size(classification_features, 2)
%         accuracies(state_features, concept_features) = classifications{state_features, concept_features}.accuracy;
%     end
% end
% 
% % Plot
% figure;
% imagesc(accuracies); colorbar;