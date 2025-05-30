%% ======= Step 2_2 - Assessing the Impact of Synchronization Errors - Plotting BER Curve - Using a Function =========
%
% Use the function addSyncErrors() to apply Synchronization errors to the signal. And Plot of the BER curve.
%
% ==============================================================================================================================================


clear; close all; clc;
addpath('../Part I - Optimal Communication chain over the ideal channel/p1_functions');
addpath('p2_functions');

%% ========================================== Load Simulation Parameters  ==========================================
Nbps = 2;                                                           % Number of bits per symbol (2^Nbps = ModOrder)
params = initParameters_2(Nbps);                                    % Initialize fixed parameters from external function

% --- Extract parameters needed ---
NumBits     = params.timing.NumBits;                                % Bits per Tx block (frame)
ModType     = params.modulation.ModulationType;                     % Modulation type ('pam' or 'qam')
ModOrder    = params.modulation.ModulationOrder;                    % Modulation order (M = 2^Nbps)
OSF         = params.sampling.OversamplingFactor;                   % Oversampling factor
SymRate     = params.timing.SymbolRate;                             % Symbol rate (symbols/sec)
Tsymb       = params.timing.SymbolPeriod;
BitRate     = params.timing.BitRate;                                % Bit rate (bits/sec)
Fs          = params.sampling.SamplingFrequency;                    % Sampling frequency (samples/sec)
Ts          = params.sampling.SamplePeriod;                         % Sample period (seconds/sample)
Beta        = params.filter.RolloffFactor;                          % RRC filter roll-off factor
NumTaps     = params.filter.NumFilterTaps;                          % Number of RRC filter taps

% --- BER Curve Parameters ---
EbN0_min_dB         = params.simulation.EbN0_min_dB;                % Minimum Eb/N0 value in dB
EbN0_max_dB         = params.simulation.EbN0_max_dB;                % Maximum Eb/N0 value in dB
EbN0_step_dB        = params.simulation.EbN0_step_dB;               % Step size for Eb/N0 sweep in dB
iterations_per_EbN0 = params.simulation.iterations_per_EbN0;        % Iterations for averaging BER at each Eb/N0 point
EbN0_domain_dB      = params.simulation.EbN0_domain_dB;             % Range of Eb/N0 values to simulate (dB)
num_EbN0_points     = length(EbN0_domain_dB);                       % Number of points on the BER curve

displayParameters(params);

% ---- CFO and sample time offset Parameters ----
Fc = 600e6;                                     % Carrier frequency in Hz
ppm = 1;
delta_cfo = ppm * 1e-6 * Fc;           
phi_0  = 0;                         % Phase offset in rad
sample_time_offset = 0;                      % Sample offset (0 samples)


% ---- Data Segmentation Parameters ----
pilot_position = 50;
pilot_size = 20;
averaging_window = 8; % Number of symbols to average over for CFO estimation


% --- Pre-allocate results array ---
ber_data = zeros(num_EbN0_points, 1);                     % Stores simulated BER for each Eb/N0 point
fprintf('\n\n========================================');
fprintf('\n    BER Curve Simulation Setup         ');
fprintf('\n========================================');
fprintf('\n Modulation       : %d-%s', ModOrder, upper(ModType));
fprintf('\n Eb/N0 Range      : %.1f dB to %.1f dB (Step: %.1f dB).', EbN0_min_dB, EbN0_max_dB, EbN0_step_dB);
fprintf('\n Bits per Tx Block: %d', NumBits);
fprintf('\n Iterations       : %d', iterations_per_EbN0);
fprintf('\n Total Bits       : %d', NumBits * iterations_per_EbN0);
fprintf('\n========================================');


%% ========================================== Transmitter Chain (Filter Only - Executed Once) ==========================================
fprintf('\n\n========================================');
fprintf('\nCommon Processing');
fprintf('\n========================================');

% -------- Generate RRC filter --------
% Generate Root Raised Cosine filter coefficients
h_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps);
fprintf('\n  Generated RRC filter (beta=%.2f, %d taps).', Beta, NumTaps);
fprintf('\n========================================');


%% ========================================== Simulation Loop over Eb/N0 ==========================================
fprintf('\n\n================================================================================');
fprintf('\nStarting Eb/N0 Loop (One Tx block per Eb/N0 point)...');
fprintf('\n================================================================================');


% Outer Loop: Iterate through each specified Eb/N0 value
for idx_EbN0 = 1:num_EbN0_points
    EbN0dB = EbN0_domain_dB(idx_EbN0);                              % Current Eb/N0 value in dB for this outer loop iteration

    % -------- 0. Transmitter Processing (Data Generation Once Per Eb/N0 Point) --------
    bit_tx      = randi([0, 1], 1, NumBits).';                        % Generate random source bits for this Eb/N0 point
    symb_tx     = mapping(bit_tx, Nbps, ModType);                   % Map bits to complex symbols
    unuseful    = symb_tx(1 : pilot_position - 1);
    pilot       = symb_tx(pilot_position : pilot_position + pilot_size - 1);
    signal_tx  = upSampler(symb_tx, OSF).';                        % Upsample by inserting zeros, transpose for filter function
    signal_tx_filtered   = applyFilter(signal_tx, h_rrc, NumTaps);          % Apply RRC pulse shaping filter
    signalPower = mean(abs(signal_tx_filtered).^2);                          % Average Power of baseband signal after pulse shaping
    Eb = signalPower / BitRate;                                     % Energy per bit Eb = P_avg / R_bit
    total_bit_errors_point = 0;                                     % Reset error counter for this Eb/N0 point
    total_bits_sim_point = 0;                                       % Reset bit counter for this Eb/N0 point
    time_vector = (0 : length(signal_tx_filtered) - 1).' * Ts;
    time_vector_symb = (0 : length(symb_tx) - 1).' * Tsymb;

    % --- Inner Loop for Averaging ---
    % Run multiple iterations using the SAME transmitted signal (signal_tx)
    % but adding a DIFFERENT random noise instance each time.
    % This allows for averaging the BER over multiple noise realizations to have a more accurate estimate of the BER.
    
    for iter = 1:iterations_per_EbN0

        % -------- 1. AWGN Channel --------
        signal_tx_noisy = addAWGN(signal_tx_filtered, Eb, EbN0dB, OSF, SymRate);     % Add Additive White Gaussian Noise based on defined Eb and EbN0dB
        signal_tx_offset = signal_tx_noisy .* exp(1j * (2 * pi * delta_cfo * time_vector + phi_0));


        % -------- 2. Receiver Chain --------
        signal_rx_matched_filtered = applyFilter(signal_tx_offset, h_rrc, NumTaps);   % Apply matched filter (same RRC filter)
        symb_rx = downSampler(signal_rx_matched_filtered, OSF);                       % Downsample to symbol rate, transpose back
        [toa, cfo_hat] = frameFreqAcquisition(pilot, symb_rx, averaging_window, Tsymb);
        symb_rx = symb_rx .* exp(-1j * (2 * pi * cfo_hat * time_vector_symb));    % Compensate for CFO and phase offset
        bit_rx = demapping(symb_rx, Nbps, ModType);                                             % Demap received symbols to bits

        % -------- 3. Calculate Errors for this iteration --------
        num_errors_iter = sum(bit_tx ~= bit_rx);                            % Compare transmitted and received bits to count errors
        bits_iter = length(bit_tx);                                         % Number of bits transmitted in this block (should equal NumBits)

        total_bit_errors_point = total_bit_errors_point + num_errors_iter;  % Accumulate errors over iterations for this EbN0
        total_bits_sim_point = total_bits_sim_point + bits_iter;            % Accumulate total bits simulated for this EbN0

    end

    % -------- Calculate Average BER for this EbN0 point --------
    ber_data(idx_EbN0) = total_bit_errors_point / total_bits_sim_point;     % Calculate and store the average BER

    % Print results for the each Eb/N0 point
    fprintf('\n  Eb/N0 = %5.1f dB :  Eb = %.2e,   Total Bits = %8d,   Total Errors = %6d,  Avg BER = %.3e', ...
        EbN0dB, Eb, total_bits_sim_point, total_bit_errors_point, ber_data(idx_EbN0));

end
fprintf('\n================================================================================');


fprintf('\n\n========================================');
fprintf('\nEb/N0 Loop Complete.');
fprintf('\n========================================');


%% ====================== Generate Plots  =======================
plotBERCurve(ber_data, params);
bits_to_plot = min(params.timing.NumBits, 100 * Nbps); 
plotConstellation_Tx_Rx(ModOrder, ModType, symb_tx, symb_rx);
plotBitstream_Tx_Rx(bit_tx, bit_rx, bits_to_plot);
plotFilterCharacteristics(h_rrc, Beta, Fs, OSF);
plotPSD_Tx_Rx(signal_tx_filtered, signal_rx_matched_filtered, Fs);
plotBasebandFrequencyResponse(signal_tx_filtered, signal_rx_matched_filtered, Fs);

fprintf('\n\n========================================');
fprintf('\nPlotting complete.');
fprintf('\n========================================');