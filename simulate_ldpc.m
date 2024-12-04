function [gain] = simulate_ldpc(num_trials, snr_range, ratio, cfgLDPCEnc, cfgLDPCDec, num_iter, target)
    %{
    input arguments:
    num_trials: number of simulation trials per SNR point
    snr_range: SNR range
    ratio: variance ratio in dB
    cfgLDPCEnc: LDPC encoder configuration
    cfgLDPCDec: LDPC decoder configuration
    num_iter: number of iterations
    target: target BER value
    output arguments:
    gain: gain in dB
    %}

    BER_known = zeros(1, length(snr_range));        % BER with known variance
    BER_unknow = zeros(1, length(snr_range));       % BER with unknown variance
    info_length = cfgLDPCEnc.NumInformationBits;    % number of information bits
    block_length = cfgLDPCDec.BlockLength;          % block length
    num_bits = num_trials * block_length;           % number of bits
    rate = info_length / block_length;              % code rate
    num_tasks = 125;                                % number of tasks
    trials_per_task = num_trials / num_tasks;       % trials per task

    for snr_index = 1:length(snr_range)             % loop over SNR points
        tic;
        snr = snr_range(snr_index);                 % current SNR value
        variance_vector = zeros(block_length, 1);   % variance vector
        variance1 = 1 / (2 * rate * 10^(snr / 10)); % variance at odd bits
        variance2 = variance1 * ratio;              % variance at even bits
        variance_vector(1:2:end) = variance1;       % interleaved variance
        variance_vector(2:2:end) = variance2;       % interleaved variance
        variance_for_decoding_known = variance_vector;
        variance_for_decoding_unknow = (variance1 + variance2) / 2;
        task_errors_known = zeros(1, num_tasks);    % errors with known variance per task
        task_errors_unknow = zeros(1, num_tasks);   % errors with unknown variance per task
        parfor task_id = 1:num_tasks                % loop over tasks
            local_errors_known = 0;
            local_errors_unknow = 0;
            for trial = 1:trials_per_task
                info_bits = randi([0, 1], info_length, 1);          % generate random information bits
                codeword = ldpcEncode(info_bits, cfgLDPCEnc);       % encode information bits
                codeword_m = 1 - 2 * codeword;                      % BPSK modulation
                noise_vector = sqrt(variance_vector) .* randn(block_length, 1); % generate noise vector
                recieve = codeword_m + noise_vector;                % received signal
                % decoding with known variance
                llr_known = (2 * recieve) ./ variance_for_decoding_known;
                decoded_known = ldpcDecode(llr_known, cfgLDPCDec, num_iter, "OutputFormat", "whole");
                local_errors_known = local_errors_known + biterr(codeword, decoded_known);
                % decoding with unknown variance
                llr_unknow = (2 * recieve) ./ variance_for_decoding_unknow;
                decoded_unknow = ldpcDecode(llr_unknow, cfgLDPCDec, num_iter, "OutputFormat", "whole");
                local_errors_unknow = local_errors_unknow + biterr(codeword, decoded_unknow);
            end
            task_errors_known(task_id) = local_errors_known;
            task_errors_unknow(task_id) = local_errors_unknow;
        end
        errors_known = sum(task_errors_known);
        errors_unknow = sum(task_errors_unknow);
        BER_known(1, snr_index) = errors_known / num_bits;
        BER_unknow(1, snr_index) = errors_unknow / num_bits;
        toc;
    end
    [known_unique, known_idx] = unique(BER_known, "stable");
    [unknow_unique, unknow_idx] = unique(BER_unknow, "stable");
    snr_known = snr_range(known_idx);
    snr_unknow = snr_range(unknow_idx);
    known = interp1(known_unique, snr_known, target, "pchip", "extrap");
    unknow = interp1(unknow_unique, snr_unknow, target, "pchip", "extrap");
    gain = unknow - known;
end
