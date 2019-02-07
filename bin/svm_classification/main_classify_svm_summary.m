%% Description

%{

Summarise power and phi-3 SVM classification results

%}

%% Setup

results_location = 'results/';
measures = {'power', 'coherence', 'phi3'};

%% Within-fly classification

results = cell(1, 2);

results{1} = load([results_location 'power_svm_within.mat']); % Power
results{2} = load([results_location 'coherence_svm_within.mat']); % Coherence
results{3} = load([results_location 'phi3_svm_within.mat']); % Phi-3

% Average accuracies across flies
for measure = 1 : length(results)
    nFlies = length(results{measure}.classifications);
    accuracies = zeros(nFlies, 1);
    
    for fly = 1 : nFlies
        accuracies(fly) = results{measure}.classifications{fly}.accuracy;
    end
    
    results{measure}.accuracies = accuracies;
    results{measure}.accuracy = mean(accuracies);
    
    disp([measures{measure} ' - mean within-fly accuracy across N=13 flies: ' num2str(results{measure}.accuracy)]);
end

%% Across-fly classification

results = cell(1, 2);

results{1} = load([results_location 'power_svm_across.mat']); % Power
results{2} = load([results_location 'coherence_svm_across.mat']); % Coherence
results{3} = load([results_location 'phi3_svm_across.mat']); % Phi-3

for measure = 1 : length(results)
    disp([measures{measure} ' - across-fly accuracy: ' num2str(results{measure}.classifications.accuracy)]);
end

