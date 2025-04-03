function ber_averages_datas = generateBERData(params)
    %   Runs BER simulation for ONE specific parameter set.
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
    EbN0_domain_dB      = params.simulation.EbN0_domain_dB;
    iterations_per_EbN0 = params.simulation.iterations_per_EbN0;
    num_EbN0_points     = length(EbN0_domain_dB);

    % ---- Setup ----
    ber_averages_datas = zeros(1, num_EbN0_points);

    fprintf('========================================');
    fprintf('\n    Simulation Initialization         ');
    fprintf('\n========================================');
    fprintf('\n Modulation       : %d-%s', params.modulation.ModulationOrder, upper(ModType));
    fprintf('\n Eb/N0 Range      : %.1f dB to %.1f dB (Step: %.1f dB)', EbN0_min_dB, EbN0_max_dB, EbN0_step_dB);
    fprintf('\n Bits per Tx Block: %d', NumBits);
    fprintf('\n Iterations       : %d', iterations_per_EbN0);
    fprintf('\n Total Bits       : %d', NumBits * iterations_per_EbN0);
    fprintf('\n========================================');

    % -------- Generate RRC filter (once per simulation run) --------
    h_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps);
    fprintf('\n\n========================================');
    fprintf('\n      RRC Filter Generation           ');
    fprintf('\n========================================');
    fprintf('\n Generated RRC filter (beta=%.2f, %d taps).', Beta, NumTaps);
    fprintf('\n========================================');

    % -------- Simulation Loop over Eb/N0 --------
    fprintf('\n\n================================================================================');
    fprintf('\n       Eb/N0 Simulation Loop          ');
    fprintf('\n================================================================================');
    fprintf('\n Simulating Eb/N0 points:');
    for idx_EbN0 = 1:num_EbN0_points
        EbN0dB = EbN0_domain_dB(idx_EbN0);

        % -------- 0. Transmitter Processing (Data Generation Once Per Eb/N0 Point) --------
        bit_tx = randi([0, 1], 1, NumBits);
        symb_tx = mapping(bit_tx, Nbps, ModType);
        symb_tx_up = upSampler(symb_tx, OSF).';
        signal_tx = applyFilter(symb_tx_up, h_rrc, NumTaps);
        signalPower = mean(abs(signal_tx).^2);
        Eb = signalPower / BitRate;             % Energy per bit based on AVERAGE power after shaping

        total_bit_errors_point  = 0;
        total_bits_sim_point    = 0;

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
            num_errors_iter = sum(bit_tx ~= bit_rx);                    % Ensure comparison length matches Tx
            bits_iter = length(bit_tx);                                 % Should equal NumBits

            total_bit_errors_point = total_bit_errors_point + num_errors_iter;
            total_bits_sim_point = total_bits_sim_point + bits_iter;
        end

        ber_averages_datas(idx_EbN0) = total_bit_errors_point / total_bits_sim_point; % Calculate and store the average BER
            
        fprintf('\n  Eb/N0 = %5.1f dB :  Eb = %.2e,   Total Bits = %8d,   Total Errors = %6d,  Avg BER = %.3e', ...
             EbN0dB, Eb, total_bits_sim_point, total_bit_errors_point, ber_averages_datas(idx_EbN0));

    end
    fprintf('\n================================================================================');
    fprintf('\n\n========================================');
    fprintf('\n       Generation Complete            ');
    fprintf('\n========================================\n');

end