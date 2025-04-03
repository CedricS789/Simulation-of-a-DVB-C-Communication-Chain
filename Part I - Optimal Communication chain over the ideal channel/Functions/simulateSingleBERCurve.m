function [EbN0_domain_dB, ber_values, actual_bits_simulated, actual_errors_accumulated] = simulateSingleBERCurve(params)
    %SIMULATESINGLEBERCURVE Runs BER simulation for ONE specific parameter set.
    %   Generates data, runs through Tx -> Channel -> Rx chain over a range
    %   of Eb/N0 values for a given modulation and system configuration.
    %
    %   Inputs:
    %       params - Struct containing all parameters for this run (from coreParameters).
    %
    %   Outputs:
    %       EbN0_domain_dB - Vector of Eb/N0 values used (dB).
    %       ber_values     - Vector of corresponding simulated BER values.
    %       actual_bits_simulated - Total bits simulated per Eb/N0 point.
    %       actual_errors_accumulated - Total errors accumulated per Eb/N0 point.

    % --- Extract parameters needed ---
    Nbps        = params.modulation.Nbps;
    NumBits     = params.timing.NumBits;
    ModType     = params.modulation.ModulationType;
    OSF         = params.sampling.OversamplingFactor;
    SymRate     = params.timing.SymbolRate;
    BitRate     = params.timing.BitRate;
    Beta        = params.filter.RolloffFactor;
    NumTaps     = params.filter.NumFilterTaps;
    EbN0_min_dB = params.simulation.EbN0_min_dB;
    EbN0_max_dB = params.simulation.EbN0_max_dB;
    EbN0_step_dB = params.simulation.EbN0_step_dB;
    iterations_per_EbN0 = params.simulation.iterations_per_EbN0;

    % --- Simulation setup ---
    EbN0_domain_dB = EbN0_min_dB:EbN0_step_dB:EbN0_max_dB;
    num_EbN0_points = length(EbN0_domain_dB);
    ber_values = zeros(1, num_EbN0_points);
    actual_bits_simulated = zeros(1, num_EbN0_points);
    actual_errors_accumulated = zeros(1, num_EbN0_points);

    fprintf('\n\nStarting Simulation for %d-%s...', params.modulation.ModulationOrder, upper(ModType));
    fprintf('\n  Eb/N0 Range: %.1f dB to %.1f dB (Step: %.1f dB)', EbN0_min_dB, EbN0_max_dB, EbN0_step_dB);
    fprintf('\n  Bits per Tx Block (fixed per Eb/N0 point): %d', NumBits);
    fprintf('\n  Iterations per Eb/N0 point: %d', iterations_per_EbN0);
    fprintf('\n  Total Bits per Eb/N0 point: %d', NumBits * iterations_per_EbN0);

    % -------- Generate RRC filter (once per simulation run) --------
    h_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps);
    fprintf('\n  Generated RRC filter (beta=%.2f, %d taps).', Beta, NumTaps);

    % -------- Simulation Loop over Eb/N0 --------
    fprintf('\n  Simulating Eb/N0 points:');
    for idx_EbN0 = 1:num_EbN0_points
        EbN0dB = EbN0_domain_dB(idx_EbN0);

        % -------- 0. Transmitter Processing (Data Generation Once Per Eb/N0 Point) --------
        bit_tx = randi([0, 1], 1, NumBits);
        symb_tx = mapping(bit_tx, Nbps, ModType);
        symb_tx_up = upSampler(symb_tx, OSF).';
        signal_tx = applyFilter(symb_tx_up, h_rrc, NumTaps);
        signalPower = mean(abs(signal_tx).^2);
        Eb = signalPower / BitRate; % Energy per bit based on AVERAGE power after shaping

        total_bit_errors_point = 0;
        total_bits_sim_point = 0;

        % --- Inner Loop for Averaging ---
        for iter = 1:iterations_per_EbN0
            % -------- 1. AWGN Channel --------
             % Ensure addAWGN uses the calculated Eb for THIS signal
            signal_tx_noisy = addAWGN(signal_tx, Eb, EbN0dB, OSF, SymRate);

            % -------- 2. Receiver Chain --------
            signal_rx_filtered = applyFilter(signal_tx_noisy, h_rrc, NumTaps);
            symb_rx = downSampler(signal_rx_filtered, OSF).';
            bit_rx = demapping(symb_rx, Nbps, ModType);
            bit_rx = bit_rx(:).';

            % -------- 3. Calculate Errors for this iteration --------
            num_errors_iter = sum(bit_tx ~= bit_rx(1:length(bit_tx))); % Ensure comparison length matches Tx
            bits_iter = length(bit_tx); % Should equal NumBits

            total_bit_errors_point = total_bit_errors_point + num_errors_iter;
            total_bits_sim_point = total_bits_sim_point + bits_iter;
        end

        % -------- Calculate Average BER for this Eb/N0 point --------
        if total_bits_sim_point > 0
            ber_values(idx_EbN0) = total_bit_errors_point / total_bits_sim_point;
        else
            ber_values(idx_EbN0) = NaN; % Avoid division by zero if no bits were simulated
        end

        actual_bits_simulated(idx_EbN0) = total_bits_sim_point;
        actual_errors_accumulated(idx_EbN0) = total_bit_errors_point;

        % Print progress marker
        if mod(idx_EbN0, 5) == 0 || idx_EbN0 == num_EbN0_points
            fprintf('.');
        end
%         fprintf('\n  Eb/N0 = %5.1f dB :  Eb = %.2e,   Total Bits = %8d,   Total Errors = %6d,  Avg BER = %.3e', ...
%             EbN0dB, Eb, total_bits_sim_point, total_bit_errors_point, ber_values(idx_EbN0));

    end
    fprintf(' Done.\n');

end % function simulateSingleBERCurve

% --- Add required helper functions below or ensure they are on the MATLAB path ---
% function out = mapping( ... )
% function out = demapping( ... )
% function out = upSampler( ... )
% function out = downSampler( ... )
% function y = applyFilter( ... )
% function coeffs = rrcFilter( ... )
% function noisySignal = addAWGN( ... )
% (You should already have these defined)