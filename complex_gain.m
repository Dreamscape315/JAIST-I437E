clear;
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end
parpool('local', 128);
num_trials = 100000; %100 times of target ber

%------with knowledge------%
target = 1e-3;

ratio_db = 1:15;
ratio = (10.^(ratio_db ./ 10));

snr_range = -3:18;
[cfgEnc, cfgDec] = getProtoMatrix(648,324);  %depend on your code block length
gain = zeros(length(ratio), 1);
for ratio_index = 1:length(ratio)
    r = ratio(ratio_index);
    gain(ratio_index) = simulate_ldpc(num_trials, snr_range, r, cfgEnc, cfgDec, 50, target);
    fprintf('ratio %.2f complete\n', r);
end

save('gain.mat','gain');
