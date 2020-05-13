%% Load data

load('split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim.mat');

%%

fly = 4;
condition = 2;
trial = 1;
channels = [1 3 5];

%%

data_shape = size(fly_data);

data = reshape(permute(fly_data, [2 4 1 3 5]), [data_shape(2) data_shape(4) data_shape(1)*data_shape(3)*data_shape(5)]);
data = permute(data, [3 1 2]);

%%

addpath('CCM_L_M/');

% parameters
tau = 1; % time step 
E   = 2; % dimension of reconstruction
LMN = 5; % number of neigborhoods for L and M methods
[ SugiC , SugiR , LM , SugiY , SugiX , origY , origX ]=SugiLM(fly_data(:,14,1,1,1),fly_data(:,15,1,1,1),tau,E,LMN);        
% results
disp('Sugiharas CMM correlation for estimate of X and original X in coupled case')
SugiC(1)
disp('Sugiharas CMM correlation for estimate of Y and original Y in coupled case')
SugiC(2)
plot(SugiY,origY,'ro',SugiX,origX,'b*')
title('Estimated vs. original data in coupled case')
xlabel('Estimated data') 
ylabel('Original data') 
legend('Y','X')

%% Project 3 channel time-series into 3D manifold

figure;

plot3(...
    data((1:2247), channels(1), trial, fly, condition),...
    data((2:2248), channels(1), trial, fly, condition),...
    data((4:2250), channels(1), trial, fly, condition)...
    );

%%

figure;

plot3(...
    data(:, channels(1), fly),...
    data(:, channels(2), fly),...
    data(:, channels(3), fly)...
    );

%%
hold on;
for t = 1 : 2247%size(fly_data)
    scatter3(...
        fly_data(t, channels(1), trial, fly, condition),...
        fly_data(t+1, channels(1), trial, fly, condition),...
        fly_data(t+2, channels(1), trial, fly, condition),...
        '.'...
        );
    hold on;
    drawnow;
end

%%
hold on;
for t = 1 : 2246%size(fly_data)
    plot3(...
        [fly_data(t, channels(1), trial, fly, condition) fly_data(t+1, channels(1), trial, fly, condition)],...
        [fly_data(t, channels(2), trial, fly, condition) fly_data(t+1, channels(2), trial, fly, condition)],...
        [fly_data(t, channels(3), trial, fly, condition) fly_data(t+1, channels(3), trial, fly, condition)],...
        'Color', [0 0 0 0.1]...
        );
    hold on;
    drawnow;
end

%%
hold on;
for t = 1 : 2246%size(fly_data)
    plot3(...
        [data(t, channels(1), fly) data(t+1, channels(1), fly)],...
        [data(t, channels(2), fly) data(t+1, channels(2), fly)],...
        [data(t, channels(3), fly) data(t+1, channels(3), fly)],...
        'Color', [0 0 0 0.1]...
        );
    hold on;
    drawnow;
end

%%

hold on;
for t = 1 : 2246%size(fly_data)
    plot3(...
        [fly_data(t, channels(1), trial, fly, condition) fly_data(t+1, channels(1), trial, fly, condition)],...
        [fly_data(t+1, channels(1), trial, fly, condition) fly_data(t+1+1, channels(1), trial, fly, condition)],...
        [fly_data(t+2, channels(1), trial, fly, condition) fly_data(t+2+1, channels(1), trial, fly, condition)],...
        'Color', [0 0 0 0.1]...
        );
    hold on;
    drawnow;
end
