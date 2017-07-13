% Jamin's phi-star tutorial: https://github.com/wjmn/phi_demo/blob/master/phi_star/DEMO_Phi_Calculations.ipynb

T = 406;
tau = 6;

data_tmp = repmat(linspace(1, 3054, 3054), 118, 1);
current_data_chn_1 = data_tmp;
current_data_chn_2 = data_tmp;

current_bin = [current_data_chn_1(1, 1:T); current_data_chn_2(1, 1:T)];

state_present = current_bin(:, (1 + tau): T);
state_past = current_bin(:, 1:(T - tau));

mean_present = mean(state_present, 2);
mean_past = mean(state_past, 2);

mean_present_array = repmat(mean_present, 1, (T - tau));
mean_past_array = repmat(mean_past, 1, (T - tau));

normaliser = (T - tau) - 1;

dev_present = state_present - mean_present_array;
dev_past = state_past - mean_past_array;

cov_present = (dev_present * dev_present') / normaliser;
cov_cross = (dev_present * dev_past') / normaliser;

% Jamin's calculation and Cov_comp_sample.m give different results

disp('Cov_X');
disp(cov_present);
disp('Cov_XY');
disp(cov_cross);
disp('Cov_comp_sample');
[tut_cov_X, tut_cov_XY] = Cov_comp_sample(current_bin, tau);
disp('cov_X');
disp(tut_cov_X);
disp('cov_XY');
disp(tut_cov_XY);