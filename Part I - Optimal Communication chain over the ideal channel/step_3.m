%%   =================== Step 3 - Simulation for BER Curve Generation ===================
%
%       Orchestrates the simulation over a range of Eb/N0 values to generate
%       a BER curve by simulating Tx -> AWGN Channel -> Rx repeatedly.
%       Uses standardized Eb definition based on symbol energy Es=1.
%       Generates one random bit sequence per Eb/N0 point, used for all averaging iterations.
%
%   =======================================================================================

clear; close all; clc;

%% ========================================== Simulation Parameters  ==========================================
params = coreParameters();                                          % Initialize fixed parameters from external function

% --- Extract parameters needed ---
Nbps        = params.modulation.Nbps;                               % Number of bits per symbol
NumBits     = params.timing.NumBits;                                % Bits per Tx block (frame)
ModType     = params.modulation.ModulationType;                     % Modulation type ('pam' or 'qam')
ModOrder    = params.modulation.ModulationOrder;                    % Modulation order (M = 2^Nbps)
OSF         = params.sampling.OversamplingFactor;                   % Oversampling factor
SymRate     = params.timing.SymbolRate;                             % Symbol rate (symbols/sec)
BitRate     = params.timing.BitRate;                                % Bit rate (bits/sec)
Fs          = params.sampling.SamplingFrequency;                    % Sampling frequency (samples/sec)
Beta        = params.filter.RolloffFactor;                          % RRC filter roll-off factor
NumTaps     = params.filter.NumFilterTaps;                          % Number of RRC filter taps

% --- BER Curve Parameters ---
EbN0_min_dB         = params.simulation.EbN0_min_dB;                % Minimum Eb/N0 value in dB
EbN0_max_dB         = params.simulation.EbN0_max_dB;                % Maximum Eb/N0 value in dB
EbN0_step_dB        = params.simulation.EbN0_step_dB;               % Step size for Eb/N0 sweep in dB
iterations_per_EbN0 = params.simulation.iterations_per_EbN0;        % Iterations for averaging BER at each Eb/N0 point
EbN0_domain_dB      = EbN0_min_dB:EbN0_step_dB:EbN0_max_dB;         % Range of Eb/N0 values to simulate (dB)
num_EbN0_points     = length(EbN0_domain_dB);                       % Number of points on the BER curve

% --- Pre-allocate results array ---
ber_values = zeros(1, num_EbN0_points);                             % Stores simulated BER for each Eb/N0 point
actual_bits_simulated = zeros(1, num_EbN0_points);                  % Stores total bits simulated per Eb/N0 point
actual_errors_accumulated = zeros(1, num_EbN0_points);              % Stores total errors accumulated per Eb/N0 point

fprintf('\n\n------------------------------------');
fprintf('\nBER Curve Simulation Setup:');
fprintf('\n------------------------------------');
fprintf('\n  Modulation: %d-%s', ModOrder, upper(ModType));
fprintf('\n  Eb/N0 Range: %.1f dB to %.1f dB (Step: %.1f dB)', EbN0_min_dB, EbN0_max_dB, EbN0_step_dB);
fprintf('\n  Bits per Tx Block (fixed per Eb/N0 point): %d', NumBits);
fprintf('\n  Iterations per Eb/N0 point: %d', iterations_per_EbN0);
fprintf('\n  Total Bits per Eb/N0 point: %d', NumBits * iterations_per_EbN0);


%% ========================================== Transmitter Chain (Filter Only - Executed Once) ==========================================
fprintf('\n\n------------------------------------');
fprintf('\nCommon Processing');
fprintf('\n------------------------------------');

% -------- Generate RRC filter --------
h_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps);                     % Generate Root Raised Cosine filter coefficients
fprintf('\n  Generated RRC filter (beta=%.2f, %d taps).', Beta, NumTaps);


%% ========================================== Simulation Loop over Eb/N0 ==========================================
fprintf('\n\n------------------------------------');
fprintf('\nStarting Eb/N0 Loop (One Tx block per Eb/N0 point)...');
fprintf('\n------------------------------------');


% Outer Loop: Iterate through each specified Eb/N0 value
for idx_EbN0 = 1:num_EbN0_points
    EbN0dB = EbN0_domain_dB(idx_EbN0);                              % Current Eb/N0 value in dB for this outer loop iteration

    % -------- 0. Transmitter Processing (Data Generation Once Per Eb/N0 Point) --------
    bit_tx = randi([0, 1], 1, NumBits);                             % Generate random source bits for this Eb/N0 point
    symb_tx = mapping(bit_tx, Nbps, ModType);                       % Map bits to complex symbols
    symb_tx_up = upSampler(symb_tx, OSF).';                         % Upsample by inserting zeros, transpose for filter function
    signal_tx = applyFilter(symb_tx_up, h_rrc, NumTaps);            % Apply RRC pulse shaping filter
    signalPower = mean(abs(signal_tx).^2);                          % Average Power of baseband signal after pulse shaping
    Eb = signalPower / BitRate;                                     % Energy per bit Eb = P_avg / R_bit

    total_bit_errors_point = 0;                                     % Reset error counter for this Eb/N0 point
    total_bits_sim_point = 0;                                       % Reset bit counter for this Eb/N0 point

    % --- Inner Loop for Averaging ---
    % Run multiple iterations using the SAME transmitted signal (signal_tx)
    % but adding a DIFFERENT random noise instance each time.
    % This allows for averaging the BER over multiple noise realizations to have a more accurate estimate of the BER.
    for iter = 1:iterations_per_EbN0

        % -------- 1. AWGN Channel --------
        signal_tx_noisy = addAWGN(signal_tx, Eb, EbN0dB, OSF, SymRate); % Add Additive White Gaussian Noise based on defined Eb and EbN0dB

        % -------- 2. Receiver Chain --------
        signal_rx_filtered = applyFilter(signal_tx_noisy, h_rrc, NumTaps); % Apply matched filter (same RRC filter)
        symb_rx = downSampler(signal_rx_filtered, OSF).';           % Downsample to symbol rate, transpose back
        bit_rx = demapping(symb_rx, Nbps, ModType);                 % Demap received symbols to bits
        bit_rx = bit_rx(:).';                                       % Ensure bit_rx is a row vector

        % -------- 3. Calculate Errors for this iteration --------
        num_errors_iter = sum(bit_tx ~= bit_rx);                    % Compare transmitted and received bits to count errors
        bits_iter = length(bit_tx);                                 % Number of bits transmitted in this block (should equal NumBits)

        total_bit_errors_point = total_bit_errors_point + num_errors_iter;  % Accumulate errors over iterations for this EbN0
        total_bits_sim_point = total_bits_sim_point + bits_iter;            % Accumulate total bits simulated for this EbN0

    end

    % -------- Calculate Average BER for this EbN0 point --------
    ber_values(idx_EbN0) = total_bit_errors_point / total_bits_sim_point; % Calculate and store the average BER



    % Print results for the each Eb/N0 point
    fprintf('\n  Eb/N0 = %5.1f dB :  Eb = %.2e,   Total Bits = %8d,   Total Errors = %6d,  Avg BER = %.3e', ...
        EbN0dB, Eb, total_bits_sim_point, total_bit_errors_point, ber_values(idx_EbN0));

end

fprintf('\n\n------------------------------------');
fprintf('\nEb/N0 Loop Complete.');
fprintf('\n------------------------------------');


%% ========================================== Plotting Results ==========================================
fprintf('\n\n------------------------------------');
fprintf('\nPlotting BER Curve...');
fprintf('\n------------------------------------');

plotBERCurve(EbN0_domain_dB, ber_values, ModType, ModOrder, OSF, Beta);

fprintf('\n\n------------------------------------');
fprintf('\nPlotting complete.');
fprintf('\n------------------------------------');