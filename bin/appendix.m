%% PHI LME

load('analysis_results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_phithree_lmeStats.mat')

model_spec
model_null_specs

for i = 1 : 3
    compare(model_nulls{i}, model_full)
end

%% PHI-STAR LME

load('analysis_results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_phistar_lmeStats.mat')

model_spec
model_null_specs

for i = 1 : 3
    compare(model_nulls{i}, model_full)
end

%% FEEDBACK ANOVA

anova_feedback

%% PHI PHI-STAR CORRELATIONS ANOVA

anova_correlations

%% PHI PHI-STAR MIP MATCH ANOVA

for i = 1 : 3
    anova_mips{i}
end