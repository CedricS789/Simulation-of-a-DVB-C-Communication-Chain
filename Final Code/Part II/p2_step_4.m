%% =================== Step 4 - Frame and Frequency Acquisition ===================
clear; close all; clc;
addpath('../Part I/p1_functions');
addpath('p2_functions')


%% ========================================== Load Simulation Parameters  ==========================================
Nbps = 4;
params = initParameters(Nbps);
displayParameters(params);
NumBits = params.timing.NumBits;
ModType = params.modulation.ModulationType;
ModOrder= params.modulation.ModulationOrder;
OSF = params.sampling.OversamplingFactor;
SymRate = params.timing.SymbolRate;
BitRate = params.timing.BitRate;
Fs = params.sampling.SamplingFrequency;
Ts = params.sampling.SamplePeriod;
Beta = params.filter.RolloffFactor;
NumTaps = params.filter.NumFilterTaps;
iterations = params.simulation.iterations_per_EbN0;

% ---- CFO and sample time offset Parameters ----
Fc = 600e6;                     % Carrier frequency in Hz
ppm = 2;
delta_cfo = ppm * 1e-6 * Fc;    % Frequency offset in Hz (1 ppm)
phi_0  = 0;                     % Phase offset in rad
sample_time_offset = 0;         % Sample time offset

% ---- Data Segmentation Parameters ----
pilot_position = 50;
pilot_size = 20;
averaging_window = 8;           % Number of symbols to average over for CFO estimation


%% ========================================== Communication Chain ==========================================
% --- Transmitter  ---
bit_tx = randi([0, 1], 1, NumBits).';
symb_tx = mapping(bit_tx, Nbps, ModType);
unuseful = symb_tx(1 : pilot_position - 1);
pilot = symb_tx(pilot_position : pilot_position + pilot_size - 1);
signal_tx = upSampler(symb_tx, OSF).';
h_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps);
signal_tx_filtered = applyFilter(signal_tx, h_rrc, NumTaps);
signalPower_tx = mean(abs(signal_tx_filtered).^2);
Eb = signalPower_tx / BitRate;

% -- Introduce Noise --
EbN0dB = 10;
signal_tx_noisy = addAWGN(signal_tx_filtered, Eb, EbN0dB, OSF, SymRate);

% --- Introduce CFO, phase offset and sample time offset ---
time_vector = (0 : length(signal_tx) - 1).' * Ts;
signal_tx_distorted  = signal_tx_noisy .* exp(1j * (2 * pi * delta_cfo * time_vector + phi_0));    % Apply CFO to the transmitted signal


% --- Receiver Chain ---
signal_rx_matched_filtered  = applyFilter(signal_tx_distorted, h_rrc, NumTaps);                    % Apply Matched Filter
symb_rx    = downSampler(signal_rx_matched_filtered, OSF);
[toa, cfo_hat] = frameFreqAcquisition(pilot, symb_rx, averaging_window, Tsymb);
time_vector_symb = (0 : length(symb_rx) - 1).' * Tsymb;
symb_rx = symb_rx .* exp(-1j * (2 * pi * cfo_hat * time_vector_symb));                          % Compensate for CFO and phase offset
bit_rx = demapping(symb_rx, Nbps, ModType); 



%% ====================== BER  =======================
% --- BER Curve Parameters ---
EbN0_min_dB         = params.simulation.EbN0_min_dB;                % Minimum Eb/N0 value in dB
EbN0_max_dB         = params.simulation.EbN0_max_dB;                % Maximum Eb/N0 value in dB
EbN0_step_dB        = params.simulation.EbN0_step_dB;               % Step size for Eb/N0 sweep in dB
iterations_per_EbN0 = params.simulation.iterations_per_EbN0;        % Iterations for averaging BER at each Eb/N0 point
EbN0_domain_dB      = params.simulation.EbN0_domain_dB;             % Range of Eb/N0 values to simulate (dB)
num_EbN0_points     = length(EbN0_domain_dB);                       % Number of points on the BER curve

% --- Pre-allocate results array ---
ber_data = zeros(num_EbN0_points, 1);                               % Stores simulated BER for each Eb/N0 point
fprintf('\n\n========================================');
fprintf('\n    BER Curve Simulation Setup         ');
fprintf('\n========================================');
fprintf('\n Modulation       : %d-%s', ModOrder, upper(ModType));
fprintf('\n Eb/N0 Range      : %.1f dB to %.1f dB (Step: %.1f dB)', EbN0_min_dB, EbN0_max_dB, EbN0_step_dB);
fprintf('\n Bits per Tx Block: %d', NumBits);
fprintf('\n Iterations       : %d', iterations_per_EbN0);
fprintf('\n Total Bits       : %d', NumBits * iterations_per_EbN0);
fprintf('\n========================================');

% Outer Loop: Iterate through each specified Eb/N0 value
for idx_EbN0 = 1:num_EbN0_points
    EbN0dB = EbN0_domain_dB(idx_EbN0);                                      % Current Eb/N0 value in dB for this outer loop iteration

    % -------- 0. Transmitter Processing (Data Generation Once Per Eb/N0 Point) --------
    bit_tx  = randi([0, 1], 1, NumBits).';                              % Generate random source bits for this Eb/N0 point
    symb_tx  = mapping(bit_tx, Nbps, ModType);                           % Map bits to complex symbols
    unuseful = symb_tx(1 : pilot_position - 1);
    pilot = symb_tx(pilot_position : pilot_position + pilot_size - 1);
    signal_tx = upSampler(symb_tx, OSF).';                                 % Upsample by inserting zeros, transpose for filter function
    signal_tx_filtered   = applyFilter(signal_tx, h_rrc, NumTaps);          % Apply RRC pulse shaping filter
    signalPower = mean(abs(signal_tx_filtered).^2);                         % Average Power of baseband signal after pulse shaping
    Eb = signalPower / BitRate;                                             % Energy per bit Eb = P_avg / R_bit
    total_bit_errors_point = 0;                                             % Reset error counter for this Eb/N0 point
    total_bits_sim_point = 0;                                               % Reset bit counter for this Eb/N0 point
    time_vector = (0 : length(signal_tx_filtered) - 1).' * Ts;
    time_vector_symb = (0 : length(symb_tx) - 1).' * Tsymb;

    % --- Inner Loop for Averaging ---
    % Run multiple iterations using the SAME transmitted signal (signal_tx)
    % but adding a DIFFERENT random noise instance each time.
    % This allows for averaging the BER over multiple noise realizations to have a more accurate estimate of the BER.
    
    for iter = 1:iterations_per_EbN0

        % -------- 1. AWGN Channel --------
        signal_tx_noisy = addAWGN(signal_tx_filtered, Eb, EbN0dB, OSF, SymRate);     % Add Additive White Gaussian Noise based on defined Eb and EbN0dB
        signal_tx_distorted = signal_tx_noisy .* exp(1j * (2 * pi * delta_cfo * time_vector + phi_0));


        % -------- 2. Receiver Chain --------
        signal_rx_matched_filtered = applyFilter(signal_tx_distorted, h_rrc, NumTaps);     % Apply matched filter (same RRC filter)
        symb_rx = downSampler(signal_rx_matched_filtered, OSF);                         % Downsample to symbol rate, transpose back
        [toa, cfo_hat] = frameFreqAcquisition(pilot, symb_rx, averaging_window, Tsymb);
        symb_rx = symb_rx .* exp(-1j * (2 * pi * cfo_hat * time_vector_symb));          % Compensate for CFO and phase offset
        bit_rx = demapping(symb_rx, Nbps, ModType);                                     % Demap received symbols to bits

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


%% ====================== Generate Plots  =======================
plotBERCurve(ber_data, params);
bits_to_plot = min(params.timing.NumBits, 100 * Nbps); 
plotConstellation_Tx_Rx(ModOrder, ModType, symb_tx, symb_rx);
plotBitstream_Tx_Rx(bit_tx, bit_rx, bits_to_plot);
plotFilterCharacteristics(h_rrc, Beta, Fs, OSF);
plotPSD_Tx_Rx(signal_tx_filtered, signal_rx_matched_filtered, Fs);
plotBasebandFrequencyResponse(signal_tx_filtered, signal_rx_matched_filtered, Fs);



%% ===================== Robustness Analysis ====================
averaging_window_vector = (1:1:20).';
toa_hat_vector = zeros(length(averaging_window_vector), 1);
cfo_hat_vector = zeros(length(averaging_window_vector), 1);
for i = 1:length(averaging_window_vector)
    [toa_hat_vector(i), cfo_hat_vector(i)] = frameFreqAcquisition(pilot, symb_rx, averaging_window_vector(i), Tsymb);
end

figure;
plot(averaging_window_vector, abs(toa_hat_vector - pilot_position));
xlabel('Averaging Window Size');
ylabel('abs(toa/_hat - pilot/_position)');
title('TOA Estimation Error vs. Averaging Window');
grid on;

figure;
plot(averaging_window_vector, abs(cfo_hat_vector - delta_cfo));
xlabel('Averaging Window Size');
ylabel('abs(.)');
title('TOA Estimation Error vs. Averaging Window');
grid on;


pilot_size_vector = (1:1:60).';
toa_hat_vector = zeros(length(pilot_size_vector), 1);
cfo_hat_vector = zeros(length(averaging_window_vector), 1);

K = 8;

for j = 1:length(pilot_size_vector)
    pilot = symb_tx(pilot_position : pilot_position + pilot_size_vector(j) - 1);
    [toa_hat_vector(j), cfo_hat_vector(j)] = frameFreqAcquisition(pilot, symb_rx, K, Tsymb);
end
figure;
plot(pilot_size_vector, abs(toa_hat_vector - pilot_position));
xlabel('Pilot Size (Number of Symbols)');
ylabel('abs(toa/_hat - pilot/_position)');
title('TOA Estimation Error vs. Pilot Size');
grid on;

figure;
plot(pilot_size_vector, abs(cfo_hat_vector - delta_cfo));
xlabel('Pilot Size (Number of Symbols)');
ylabel('abs(toa/_hat - pilot/_position)');
title('TOA Estimation Error vs. Pilot Size');
grid on;