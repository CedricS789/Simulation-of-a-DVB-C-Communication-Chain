%%   =================== Step 3 - Semi-Modular Simulation for BER Curve ===================
%
%       This script orchestrates the simulation of a communication chain over
%       an AWGN channel to generate a Bit Error Rate (BER) curve.
%       Transmitter and Receiver logic are inlined. Noise addition and plotting
%       are done using separate functions.
%
%   =======================================================================================

clear; close all; clc;

%% ========================================== Simulation Parameters  ==========================================
params = initializeParameters(); % Initialize fixed parameters

% --- Extract parameters needed ---
Nbps        = params.modulation.Nbps;
NumBits     = params.timing.NumBits;
ModType     = params.modulation.ModulationType;
ModOrder    = params.modulation.ModulationOrder;
OSF         = params.sampling.OversamplingFactor;
SymRate     = params.timing.SymbolRate;
Fs          = params.sampling.SamplingFrequency;
Beta        = params.filter.RolloffFactor;
NumTaps     = params.filter.NumFilterTaps;
BitRate     = params.timing.BitRate;

% --- BER Curve Parameters ---
EbN0_min_dB         = params.simulation.EbN0_min_dB;
EbN0_max_dB         = params.simulation.EbN0_max_dB;
EbN0_step_dB        = params.simulation.EbN0_step_dB;
iterations_per_EbN0 = params.simulation.iterations_per_EbN0;

EbN0_domain_dB  = EbN0_min_dB:EbN0_step_dB:EbN0_max_dB; % This is all the values of Eb/N0 we’re testing
num_EbN0_points = length(EbN0_domain_dB);               % How many points are in that range?

% --- Pre-allocate results array ---
% BER_values stores the final BER for each Eb/N0 point after averaging
ber_values = zeros(1, num_EbN0_points);

fprintf('\n\n------------------------------------\n');
fprintf('BER Curve Simulation Setup:\n');
fprintf('------------------------------------\n');
fprintf('  Eb/N0 Range: %.1f dB to %.1f dB (Step: %.1f dB)\n', EbN0_min_dB, EbN0_max_dB, EbN0_step_dB);
fprintf('  Iterations per Eb/N0 point: %d\n', iterations_per_EbN0);
fprintf('  Total Bits per Eb/N0 point: %d\n', NumBits * iterations_per_EbN0);


%% ========================================== Transmitter Chain ==========================================
fprintf('\n\n------------------------------------\n');
fprintf('Transmitter Processing\n');
fprintf('------------------------------------\n');

% -------- 1. Generate Source Data Bits --------
bit_tx = randi([0, 1], 1, NumBits);
fprintf('  Generated %d random bits.\n', NumBits);

% -------- 2. Symbol Mapping --------
symb_tx = mapping(bit_tx, Nbps, ModType);
fprintf('  Mapped bits to %d %d-%s symbols.\n', length(symb_tx), ModOrder, upper(ModType));

% -------- 3. Upsampling --------
symb_tx_up = upSampler(symb_tx, OSF).'; % Transpose for filter conv
fprintf('  Upsampled symbols by a factor of %d.\n', OSF);

% -------- 4. Pulse Shaping (Transmit Filtering) --------
h_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps); % Generates filter + plots its characteristics
fprintf('  Generated RRC filter (beta=%.2f, %d taps).\n', Beta, NumTaps);
signal_tx = applyFilter(symb_tx_up, h_rrc, NumTaps);
fprintf('  Applied RRC pulse shaping filter.\n');

% -------- 5. Calculate Signal Metrics --------
signalPower_tx = mean(abs(signal_tx).^2);
Eb_tx = signalPower_tx / BitRate; % energy per bit Eb = Ptx / Rb
fprintf('  Calculated Tx Signal Power: %.2e\n', signalPower_tx);
fprintf('  Calculated Energy per Bit (Eb): %.2e\n', Eb_tx);
fprintf('Transmitter Processing Complete.\n');


%% ========================================== Simulation Loop over Eb/N0 ==========================================
% - Outer Loop (idx_EbN0): We’re stepping through a bunch of Eb/N0 values (like 0 to 15 dB) to see how the system holds up with different amounts of noise.
%   - Eb/N0 is just a fancy way of saying signal-to-noise ratio per bit. Low Eb/N0 means tons of noise; high Eb/N0 means it’s pretty clean.
% - Inner Loop (iter): For each Eb/N0, we run the simulation a few times (iterations_per_EbN0) to smooth out the randomness of the noise.
%   - Noise is all over the place, so one shot isn’t enough. Averaging gives us a solid BER number we can trust.
% - What we’re doing: For every Eb/N0, we toss some AWGN on the signal, run it through the receiver, count the bit errors, and average them to get the BER.

fprintf('\n\n------------------------------------\n');
fprintf('Starting Eb/N0 Loop...\n');
fprintf('------------------------------------\n');

hw = waitbar(0, 'Simulating BER Curve...');

% Loop over each Eb/N0 value to check out different noise levels
% This outer loop goes through all the Eb/N0 points we set up (like 16 steps from 0 to 15 dB)
% - num_EbN0_points is how many Eb/N0 values we’re testing
% - The goal here is to simulate the whole chain with varying noise to plot that BER curve
for idx_EbN0 = 1:num_EbN0_points
    % Grab the current Eb/N0 value in dB from our list (e.g., 0 dB, 1 dB, up to 15 dB)
    % - EbN0dB is what we feed into addAWGN to set the noise level; bigger number = less noise
    EbN0dB = EbN0_domain_dB(idx_EbN0);
    
    % Pop up a little progress bar so we know where we’re at
    % - Shows something like "Simulating Eb/N0 = 5.0 dB..." and updates as we go
    waitbar(idx_EbN0/num_EbN0_points, hw, sprintf('Simulating Eb/N0 = %.1f dB...', EbN0dB));

    % total_bit_errors keeps track of all the mistakes we find at this Eb/N0
    % - Start at 0 for each new Eb/N0 so we’re counting fresh
    % - We’ll add up errors from each run to figure out the average BER later
    total_bit_errors = 0; % Reset the error counter for this Eb/N0

    % --- Inner Loop for Averaging ---
    % We’re running this a few times to average out the noise craziness
    % - iterations_per_EbN0 (say, 5) is how many rounds we do for this Eb/N0
    % - Noise is random, so a single run could be way off. This keeps it steady.
    for iter = 1:iterations_per_EbN0
        % -------- 1. AWGN Channel --------
        % Throw some noise on the signal with addAWGN
        % - signal_tx is the clean signal we sent out
        % - Eb_tx is the energy per bit we calculated earlier (power divided by bit rate)
        % - EbN0dB sets how noisy it gets
        % - Fs is the sampling frequency, messing with the noise bandwidth
        % - signal_tx_noisy is what we get: the signal all messed up by AWGN
        signal_tx_noisy = addAWGN(signal_tx, Eb_tx, EbN0dB, Fs);

        % -------- 2. Receiver Chain (Inlined) --------
        % Matched Filtering
        % - Run the noisy signal through the RRC filter to clean it up a bit
        signal_rx_filtered = applyFilter(signal_tx_noisy, h_rrc, NumTaps);
        
        % Downsampling
        % - Pull it back to one sample per symbol from all those extras
        symb_rx = downSampler(signal_rx_filtered, OSF).'; % Flip it to a row vector
        
        % Symbol Demapping
        % - Turn those symbols back into bits
        bit_rx = demapping(symb_rx, Nbps, ModType)'; % Flip it to match bit_tx

        % -------- 3. Calculate Errors for this iteration --------
        % Count how many bits got flipped in this run
        % - bit_tx is what we sent, bit_rx is what we got back
        % - num_errors_iter is just the number of differences (like 3 if 3 bits are wrong)
        % - Using sum(bit_tx ~= bit_rx) to quickly tally up the errors
        num_errors_iter = sum(bit_tx ~= bit_rx);
        
        % Add this run’s errors to the total
        % - total_bit_errors builds up (e.g., 3 + 2 + 4 = 9 errors after 3 runs)
        % - We’re summing them to average later for a solid BER
        total_bit_errors = total_bit_errors + num_errors_iter; % Pile on the errors
    end % Done with the averaging runs

    % -------- Calculate Average BER for this EbN0 --------
    % Figure out how many bits we tested total
    % - NumBits is bits per run (e.g., 400 if we’ve got 100 symbols at 4 bits each)
    % - iterations_per_EbN0 is how many times we ran it (e.g., 5)
    % - total_bits_simulated is the grand total (e.g., 400 * 5 = 2000 bits)
    % - Gotta know this to turn errors into a rate
    total_bits_simulated = NumBits * iterations_per_EbN0;
    
    % Make sure we’ve got some bits to work with
    if total_bits_simulated > 0
        % Stick the average BER for this Eb/N0 into ber_values
        % - Total errors divided by total bits gives us the error rate
        % - Like, 10 errors out of 2000 bits = 0.005 (BER = 5e-3)
        % - idx_EbN0 is where this goes in our BER array
        % - This is what we’ll plot later (Section 4.3 stuff)
        ber_values(idx_EbN0) = total_bit_errors / total_bits_simulated;
    else
        % If something’s off and we’ve got no bits, just set it to NaN
        % - Keeps the plot from freaking out
        ber_values(idx_EbN0) = NaN; % No division by zero here
    end
    
    % Show what we got for this Eb/N0
    % - Prints the Eb/N0, total errors, and the average BER
    % - Looks like: "Eb/N0 =  5.0 dB: Accumulated Errors =     10, Average BER = 5.000e-03"
    % - Nice to see how it’s going as we run
    fprintf('  Eb/N0 = %5.1f dB: Accumulated Errors = %6d, Average BER = %.3e\n', ...
            EbN0dB, total_bit_errors, ber_values(idx_EbN0));
end % Wrap up the Eb/N0 loop

close(hw); % Close waitbar
fprintf('\n------------------------------------\n');
fprintf('Eb/N0 Loop Complete.\n');
fprintf('------------------------------------\n');


%% ========================================== Plotting Results ==========================================
fprintf('\n\n------------------------------------\n');
fprintf('Plotting BER Curve...\n');
fprintf('------------------------------------\n');

plotBERCurve(EbN0_domain_dB, BER_simulated, params);

fprintf('Plotting complete.\n');
fprintf('------------------------------------\n');