function hFig = plotFilterCharacteristics(h_rrc, Beta, Fs, OSF)
    %   Visualizes RRC and combined RC filter properties.
    %   Includes RRC impulse/frequency response and combined RC pulse ISI check.
    %   Inputs:
    %       h_rrc             - RRC filter coefficients (impulse response).
    %       Beta              - Roll-off factor (for annotation).
    %       Fs                - Sampling frequency (Hz).
    %       OSF               - Oversampling factor.
    %   Output:
    %       hFig              - Handle to the created figure object.
        
    hFig = figure('Name', 'Filter Analysis (RRC & RC)', 'NumberTitle', 'off', 'Position', [100, 100, 800, 700]);
    
    % Calculate basic parameters
    NumFilterTaps = length(h_rrc);
    Ts = 1/Fs;
    NFFT = 2^12;

    % ---- Subplot 1: RRC Impulse Response ----
    subplot(3, 1, 1);
    time_axis_us = (-(NumFilterTaps-1)/2 : (NumFilterTaps-1)/2) * Ts * 1e6; % Convert to μs
    plot(time_axis_us, h_rrc, 'MarkerSize', 4, 'Marker', 'o');
    grid on; box on;
    xlabel('Time (\mus)');
    ylabel('Amplitude');
    title(sprintf('RRC Filter Impulse Response (\\beta=%.2f, Taps=%d)', Beta, NumFilterTaps));
    xlim(time_axis_us([1 end]));
    
    % ---- Subplot 2: RRC Frequency Response ----
    subplot(3, 1, 2);
    [H, ~] = freqz(h_rrc, 1, NFFT, 'whole', Fs);
    H_shifted = fftshift(H);
    frequencies_MHz = (-NFFT/2 : NFFT/2 - 1) * (Fs / NFFT) / 1e6; % Vector of frequencies in MHz
    plot(frequencies_MHz, db(H_shifted), 'LineWidth', 1.5);
    grid on; box on;
    xlabel('Frequency (MHz)');
    ylabel('Magnitude (dB)');
    title('RRC Filter Frequency Response (Magnitude)');
    xlim([-Fs/2, Fs/2] / 1e6);
    ylim([-80, 10]);
    
    % ---- Subplot 3: Combined RC Response & ISI Check ----
    subplot(3, 1, 3);
    h_rc = conv(h_rrc, h_rrc);      % Combined Tx+Rx filter responseµ
    NTaps = length(h_rc);           % Number of Fiter taps
    center = floor(NTaps / 2) + 1;
    t_norm = ((1:length(h_rc)) - center) / OSF;     % Time axis normalized to symbol periods (t / Tsymb = index_offset / OSF) - (in symbol periods)
    plot(t_norm, h_rc); hold on;                    % Plot combined RC response
    
    max_plot_offset_sym = 20;                                   % How many symbol periods more or less to show stems for in the plot
    k_offsets = -max_plot_offset_sym : max_plot_offset_sym;     % Integer symbol offsets
    sample_indices = center + k_offsets * OSF;
    
    % Ensure indices are valid before accessing h_rc
    valid_mask = (sample_indices >= 1) & (sample_indices <= length(h_rc));
    valid_k = k_offsets(valid_mask);
    valid_indices = sample_indices(valid_mask);
    
    % Plot ISI points
    stem(valid_k, h_rc(valid_indices), 'r', 'filled', 'MarkerSize', 4, 'BaseValue', 0);
    grid on; box on;

    title('Combined RC Pulse & ISI Check');
    xlabel('Time (Symbol Periods, t / T_{symb})');
    ylabel('Amplitude');
    legend('RC Pulse', 'Samples at k*T_{symb}');
    xlim([-max_plot_offset_sym-0.5, max_plot_offset_sym+0.5]);
    yline(0, 'k:', 'HandleVisibility', 'off');          % Add horizontal line at y=0, excluded from legend
    hold off;
end