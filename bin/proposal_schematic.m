%{
Figure creation for research proposal

Experimental setup > data > binarised data > transition probability matrix
> phi
%}

%{
dimensions = [3 10]; % channels x samples
sample_range = [-30 60];

sample_data = rand(dimensions);
sample_data = sample_range(1) + (sample_range(2) - sample_range(1)) .* sample_data;

middle = median(sample_data(:));
middle_mat = repmat(middle, dimensions);

binarised_data = sample_data > middle_mat;

figure;
colormap('jet');
subplot(4, 1, [1 2]);
imagesc(sample_data);

subplot(4, 1, [3 4]);
imagesc(binarised_data); colorbar;
%}

% Load
load('../../flies/fly_data/trials_anesthDescAdded11092014_bPlrRerefTyp2/Analyzed_WalkingDrosDror115610072014/trials.mat');

channels = (1:3);
samples = (81787:81787+19);

sample_data = [trials.LFP];
sample_data = sample_data(channels, samples);

% Each channel is binarised based on its median value
middle = median(sample_data, 2);
middle_mat = repmat(middle, [1 length(samples)]);
%middle_mat = repmat(median(sample_data(:)), [length(channels) length(samples)]);

binarised_data = sample_data > middle_mat;

% TPM
nStates = length(channels)^2;
tpm = zeros(nStates);

transition_counter = zeros(nStates, 1);



figure;
sub1 = subplot(4, 1, [1 2]);
imagesc(sample_data); cbar = colorbar;
set(gca, 'YTick', [1 2 3], 'XTickLabel', '');
xlabel(cbar, '\muV');
colormap(sub1, 'jet');

sub2 = subplot(4, 1, [3 4]);
imagesc(binarised_data);
set(gca, 'YTick', [1 2 3]);
cbar = colorbar; caxis([0 1]);
set(cbar, 'YTick', [0.25 0.75], 'YTickLabel', {'off', 'on'});
xlabel('time sample'); ylabel('channel');
colormap(sub2, [0 0 0; 1 1 1]);