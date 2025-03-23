%% ----------------- Step 3: AWGN Channel and BER Evaluation -----------------
% Define range of Eb/N0 values to simulate (in dB)
Eb_N0_dB = 0:2:20;  % Range from 0 to 20 dB with 2 dB steps
BER = zeros(size(Eb_N0_dB)); % Initialize BER array

% Calculate signal power
sig_power = mean(abs(tx_filtered).^2);

% Loop through each Eb/N0 value
for i = 1:length(Eb_N0_dB)
    % Convert Eb/N0 from dB to linear scale
    Eb_N0 = 10^(Eb_N0_dB(i)/10);
    
    % Calculate required noise power
    % No = Eb/SNR, and we need sigma^2 = No/2 for complex noise
    Eb = sig_power * Tsym / Nbs; % Energy per bit
    No = Eb / Eb_N0;            % Noise spectral density
    noise_var = No/2;           % Noise variance for complex noise (No/2 per dimension)
    
    % Generate complex Gaussian noise
    noise = sqrt(noise_var) * (randn(size(tx_filtered)) + 1j*randn(size(tx_filtered)));
    
    % Add noise to transmitted signal
    rx_signal = tx_filtered + noise;
    
    % Pass through receiver filter (matched filter)
    rx_filtered = filter(rrc_filter, 1, rx_signal);
    
    % Downsample to symbol rate (skip filter transients)
    filter_delay = (length(rrc_filter)-1)/2;
    rx_downsampled = rx_filtered(filter_delay+1:M:end);
    rx_symbols = rx_downsampled(1:length(symb_tx)); % Ensure same length as transmitted symbols
    
    % Demapping (symbol to bit conversion)
    if Nbs > 1
        bits_rx = demapping(rx_symbols, Nbs, 'qam');
    else
        bits_rx = demapping(rx_symbols, Nbs, 'pam');
    end
    
    % Calculate BER
    num_bit_errors = sum(bits_tx ~= bits_rx);
    BER(i) = num_bit_errors / Nbt;
    
    % Display progress
    fprintf('Eb/N0 = %2.1f dB, BER = %e\n', Eb_N0_dB(i), BER(i));
end

%% ----------------- Plot: BER vs Eb/N0 -----------------
figure;
semilogy(Eb_N0_dB, BER, 'bo-', 'LineWidth', 1.5);
hold on;

% Plot theoretical BER curves for comparison
Eb_N0_theory = 0:0.1:20;
if Nbs == 1  % BPSK
    BER_theory = 0.5 * erfc(sqrt(10.^(Eb_N0_theory/10)));
elseif Nbs == 2  % QPSK
    BER_theory = 0.5 * erfc(sqrt(10.^(Eb_N0_theory/10)));
elseif Nbs == 4  % 16-QAM
    BER_theory = 3/8 * erfc(sqrt((10.^(Eb_N0_theory/10))*4/10));
elseif Nbs == 6  % 64-QAM
    BER_theory = 7/24 * erfc(sqrt((10.^(Eb_N0_theory/10))*4/42));
elseif Nbs == 8  % 256-QAM
    BER_theory = 15/64 * erfc(sqrt((10.^(Eb_N0_theory/10))*4/170));
end

semilogy(Eb_N0_theory, BER_theory, 'r-', 'LineWidth', 1);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title(['BER vs E_b/N_0 for ' num2str(2^Nbs) '-' modulation_type]);
legend('Simulation', 'Theoretical');
ylim([10^-5 1]);




% Function to run the simulation for different modulation orders
function runBERSimulation()
    % Array of bits per symbol to test
    Nbs_values = [1, 2, 4, 6, 8];  % BPSK, QPSK, 16-QAM, 64-QAM, 256-QAM
    
    % Plot to save all results
    figure;
    hold on;
    grid on;
    
    % Colors for different modulation orders
    colors = {'b', 'g', 'r', 'm', 'c'};
    
    % Loop through each modulation order
    for idx = 1:length(Nbs_values)
        Nbs = Nbs_values(idx);
        
        % Run simulation with current Nbs
        [Eb_N0_dB, BER, BER_theory, modulation_type] = simulateBER(Nbs);
        
        % Plot simulation results
        semilogy(Eb_N0_dB, BER, [colors{idx} 'o-'], 'LineWidth', 1.5, 'DisplayName', ...
            [num2str(2^Nbs) '-' modulation_type ' (Simulation)']);
        
        % Plot theoretical curve
        semilogy(0:0.1:20, BER_theory, [colors{idx} '-'], 'LineWidth', 1, 'DisplayName', ...
            [num2str(2^Nbs) '-' modulation_type ' (Theory)']);
    end
    
    xlabel('E_b/N_0 (dB)');
    ylabel('Bit Error Rate (BER)');
    title('BER vs E_b/N_0 for Different Modulation Orders');
    legend('Location', 'southwest');
    ylim([10^-5 1]);
end

% Function to simulate BER for a specific modulation order
function [Eb_N0_dB, BER, BER_theory, modulation_type] = simulateBER(Nbs)
    %% ----------------- Parameters -----------------
    Nbt = 10000;                         % Total number of bits to transmit
    bits_tx = randi([0 1], Nbt, 1);      % Generate random bit sequence of 0s and 1s

    M = 4;                               % Oversampling factor (samples per symbol)
    sym_rate = 5e6;                      % Symbol rate in symbols per second [5 Msymb/s]
    Tsym = 1 / sym_rate;                 % Symbol duration in seconds
    Fs = M * sym_rate;                   % Sampling frequency in Hz (samples per second)

    rolloff = 0.2;                       % Roll-off factor for RRC filter
    span = 12;                           % Filter span in number of symbol periods
    
    %% ----------------- Step 1: Mapping -----------------
    if Nbs > 1
        symb_tx = mapping(bits_tx, Nbs, 'qam');  % QAM mapping for Nbs > 1
        modulation_type = 'QAM';
    else
        symb_tx = mapping(bits_tx, Nbs, 'pam');  % PAM mapping for Nbs = 1
        modulation_type = 'PAM';
    end

    %% ----------------- Step 2: RRC Filtering -----------------
    symbols_up = upsample(symb_tx, M);                  % Upsampled symbol sequence
    rrc_filter = rcosdesign(rolloff, span, M, 'sqrt');  % RRC filter with square-root response
    tx_filtered = filter(rrc_filter, 1, symbols_up);    % Filtered transmit signal

    %% ----------------- Step 3: AWGN Channel and BER Evaluation -----------------
    Eb_N0_dB = 0:2:20;  % Range from 0 to 20 dB with 2 dB steps
    BER = zeros(size(Eb_N0_dB)); % Initialize BER array

    % Calculate signal power
    sig_power = mean(abs(tx_filtered).^2);

    % Loop through each Eb/N0 value
    for i = 1:length(Eb_N0_dB)
        % Convert Eb/N0 from dB to linear scale
        Eb_N0 = 10^(Eb_N0_dB(i)/10);
        
        % Calculate required noise power
        Eb = sig_power * Tsym / Nbs; % Energy per bit
        No = Eb / Eb_N0;            % Noise spectral density
        noise_var = No/2;           % Noise variance for complex noise
        
        % Generate complex Gaussian noise
        noise = sqrt(noise_var) * (randn(size(tx_filtered)) + 1j*randn(size(tx_filtered)));
        
        % Add noise to transmitted signal
        rx_signal = tx_filtered + noise;
        
        % Pass through receiver filter
        rx_filtered = filter(rrc_filter, 1, rx_signal);
        
        % Downsample to symbol rate (skip filter transients)
        filter_delay = (length(rrc_filter)-1)/2;
        rx_downsampled = rx_filtered(filter_delay+M:M:end);
        rx_symbols = rx_downsampled(1:length(symb_tx)); % Ensure same length
        
        % Demapping
        if Nbs > 1
            bits_rx = demapping(rx_symbols, Nbs, 'qam');
        else
            bits_rx = demapping(rx_symbols, Nbs, 'pam');
        end
        
        % Truncate to match length if needed
        len = min(length(bits_tx), length(bits_rx));
        
        % Calculate BER
        num_bit_errors = sum(bits_tx(1:len) ~= bits_rx(1:len));
        BER(i) = num_bit_errors / len;
        
        fprintf('Modulation: %d-%s, Eb/N0 = %2.1f dB, BER = %e\n', 2^Nbs, modulation_type, Eb_N0_dB(i), BER(i));
    end

    % Calculate theoretical BER for comparison
    Eb_N0_lin = 10.^(0:0.1:20)/10;
    if Nbs == 1  % BPSK
        BER_theory = 0.5 * erfc(sqrt(Eb_N0_lin));
    elseif Nbs == 2  % QPSK
        BER_theory = 0.5 * erfc(sqrt(Eb_N0_lin));
    elseif Nbs == 4  % 16-QAM
        BER_theory = 3/8 * erfc(sqrt((Eb_N0_lin)*4/10));
    elseif Nbs == 6  % 64-QAM
        BER_theory = 7/24 * erfc(sqrt((Eb_N0_lin)*4/42));
    elseif Nbs == 8  % 256-QAM
        BER_theory = 15/64 * erfc(sqrt((Eb_N0_lin)*4/170));
    end
end