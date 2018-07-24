%%

taus = [1 2 3 4 5 6 7 8 16];

phis = zeros(size(taus));
for tau_counter = 1 : length(taus)
    tau = taus(tau_counter);
    
    loaded = load(['results_split/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels2_globalTPM0_f01c1tauBin' num2str(tau) 's0001t1.mat']);
    
    phis(tau_counter) = loaded.phi.phi;
    
end

