function hFig = plotBasebandFrequencyResponse(signal_tx, signal_rx, Fs)
%   Plots magnitude frequency response of two signals.
%   Uses FFT to calculate the spectrum and plots the magnitude in dB,
%   normalized to the peak (0 dB).
%
%   Inputs:
%       signal_tx - Time-domain transmitted signal vector.
%       signal_rx - Time-domain received signal vector.
%       Fs        - Sampling frequency (Hz).
%
%   Output:
%       hFig      - Handle to the created figure object.

    % --- Configuration ---
    txColor         = [0, 0.4470, 0.7410];      % Blue
    rxColor         = [0.8500, 0.3250, 0.0980]; % Red
    lineWidth       = 1.5;
    freqUnitsFactor = 1e6;                      % Hz -> MHz
    freqUnitsLabel  = 'MHz';
    yLimitMinDb     = -80;                      % Lower Y-axis limit in dB

    % --- Calculate Frequency Response (Tx) ---
    N_tx = length(signal_tx);
    % Calculate FFT, shift zero frequency to center
    fft_tx_shifted = fftshift(fft(signal_tx(:))); % Use column vector
    % Calculate magnitude
    mag_tx = abs(fft_tx_shifted);
    % Convert to dB, normalized to peak (0 dB). Add eps for log(0) stability.
    db_tx = 20 * log10(mag_tx / max(mag_tx + eps) + eps);
    % Create frequency axis
    freq_tx = (-N_tx/2 : N_tx/2 - 1) * (Fs / N_tx);

    % --- Calculate Frequency Response (Rx) ---
    N_rx = length(signal_rx);
    fft_rx_shifted = fftshift(fft(signal_rx(:)));
    mag_rx = abs(fft_rx_shifted);
    db_rx = 20 * log10(mag_rx / max(mag_rx + eps) + eps);
    freq_rx = (-N_rx/2 : N_rx/2 - 1) * (Fs / N_rx); % Recalculate in case lengths differ

    % --- Plotting ---
    hFig = figure;
    ax = axes('Parent', hFig);
    hold(ax, 'on');

    % Plot Tx Signal Spectrum
    plot(ax, freq_tx / freqUnitsFactor, db_tx, ...
         'Color', txColor, 'LineWidth', lineWidth, 'DisplayName', 'Transmitted Signal');

    % Plot Rx Signal Spectrum
    plot(ax, freq_rx / freqUnitsFactor, db_rx, ...
         'Color', rxColor, 'LineWidth', lineWidth, 'DisplayName', 'Received Signal');

    % --- Formatting ---
    grid(ax, 'on');
    xlabel(ax, sprintf('Frequency (%s)', freqUnitsLabel));
    ylabel(ax, 'Normalized Magnitude (dB)');
    title(ax, 'Baseband Signal Frequency Response');
    legend(ax, 'show', 'Location', 'best');
    xlim(ax, [-Fs/2, Fs/2] / freqUnitsFactor); % Set x-axis limits to +/- Fs/2
    ylim(ax, [yLimitMinDb, 5]);               % Set y-axis limits (dB)
    box(ax, 'on');
    hold(ax, 'off');

end