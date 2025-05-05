function hFig = plotPSD_Tx_Rx(signal_tx, signal_rx, SamplingFrequency)
     %PLOTPSD_TX_RX Plots the Power Spectral Density (PSD) of transmitted and received signals.
     %   HFIG = PLOTPSD_TX_RX(SIGNAL_TX, SIGNAL_RX, SAMPLINGFREQUENCY)
     %   calculates and plots the PSD for both the transmitted signal (SIGNAL_TX)
     %   and the received signal (SIGNAL_RX) on the same axes. The PSD shows
     %   how the power of the signal is distributed over frequency.
     %
     %   This visualization helps in verifying the spectral characteristics of the
     %   signals, such as the effectiveness of pulse shaping in limiting bandwidth
     %   (for Tx) and the impact of noise and filtering on the received spectrum (for Rx).
     %
     %   The PSD is estimated using Welch's overlapped segment averaging method
     %   (via MATLAB's 'pwelch' function). Plotting appearance and pwelch
     %   parameters are controlled internally.
     %
     %   Inputs:
     %       signal_tx         - A vector (row or column) containing the time-domain
     %                           transmitted signal waveform (e.g., after pulse shaping).
     %       signal_rx         - A vector containing the time-domain received signal
     %                           waveform (e.g., after channel and matched filter).
     %                           Should operate at the same sampling frequency as signal_tx.
     %       SamplingFrequency - The sampling frequency (Fs) in Hz common to both
     %                           signal_tx and signal_rx. Required for correct PSD calculation.
     %
     %   Output:
     %       hFig              - A handle to the created figure object.

     % =====================================================================
     % == Plotting Parameters (Locally Defined) ==
     % =====================================================================
     % --- Colors ---
     txColor = [0, 0.4470, 0.7410];      % Standard MATLAB blue for transmitted signal PSD
     rxColor = [0.8500, 0.3250, 0.0980]; % Standard MATLAB red for received signal PSD
     % --- Line ---
     lineWidth = 1.5;                    % Line width for PSD plots
     % --- General Figure ---
     figureName          = 'Power Spectral Density Comparison'; % Figure window title
     figureNumberTitle   = 'off';        % Hide 'Figure X:' prefix ('on' or 'off')
     % --- Axes Properties ---
     plotTitle           = 'PSD Comparison: Transmitted vs. Received Signals'; % Main title for the plot
     yAxisLabel          = 'Normalized Power Spectral Density (dB/Hz)'; % Y-axis label
     yLimitDb            = [-100, 5];    % Y-axis limits in dB for better visualization of spectral shape/sidelobes
     freqUnitsFactor     = 1e6;          % Factor to convert Hz to MHz for frequency axis
     freqUnitsLabel      = 'MHz';        % Label for Megahertz units
     xAxisLabelFormat    = 'Frequency (%s)'; % Format string for x-axis label
     % --- Legend ---
     legendTextTx        = 'Transmitted Signal (Post-Pulse Shaping)'; % Legend entry for Tx PSD
     legendTextRx        = 'Received Signal (Post-Matched Filter)'; % Legend entry for Rx PSD
     legendLocation      = 'best';       % Location of the legend box
     % --- pwelch Calculation Settings (Welch's method parameters) ---
     pwelchWindowType    = @hamming;     % Window function handle (e.g., @hamming, @hann, @rectwin)
     pwelchDefaultWindowLen = 512;       % Default window length (segment size) for pwelch
     pwelchOverlapPercent= 50;           % Percentage overlap between segments (e.g., 50%)
     pwelchMinNfft       = 1024;         % Minimum number of FFT points (for frequency resolution)
     pwelchSpectrumType  = 'psd';        % Type of estimate ('psd' or 'power')
     pwelchRange         = 'centered';   % Frequency range ('centered' for [-Fs/2, Fs/2], 'onesided' for [0, Fs/2])
     % =====================================================================

     % --- Figure and Axes Creation ---
     % Create the main figure window.
     hFig = figure('Name', figureName, 'NumberTitle', figureNumberTitle);
     % Create axes within the figure.
     axPSD = axes('Parent', hFig);
     % Hold the axes to plot multiple lines.
     hold(axPSD, 'on');

     % --- Prepare pwelch Parameters ---
     % Determine appropriate window length, overlap, and NFFT for Tx signal.
     % Use default length or signal length if shorter.
     pwelch_window_len_tx = min(length(signal_tx), pwelchDefaultWindowLen);
     % Calculate number of overlapping samples based on percentage.
     pwelch_noverlap_tx = floor(pwelchOverlapPercent / 100 * pwelch_window_len_tx);
     % Generate the window vector.
     pwelch_window_tx = pwelchWindowType(pwelch_window_len_tx);
     % Determine NFFT: use minimum setting or next power of 2 of window length, whichever is larger.
     nfft_tx = max(pwelchMinNfft, 2^nextpow2(pwelch_window_len_tx));

     % Determine appropriate window length, overlap, and NFFT for Rx signal (may differ if lengths differ).
     pwelch_window_len_rx = min(length(signal_rx), pwelchDefaultWindowLen);
     pwelch_noverlap_rx = floor(pwelchOverlapPercent / 100 * pwelch_window_len_rx);
     pwelch_window_rx = pwelchWindowType(pwelch_window_len_rx);
     nfft_rx = max(pwelchMinNfft, 2^nextpow2(pwelch_window_len_rx));

     % --- Calculate and Plot Transmitted Signal PSD ---
     % Calculate PSD using Welch's method for the transmitted signal.
     % Ensure signal_tx is a column vector for pwelch.
     [psd_tx, freq_tx] = pwelch(signal_tx(:), pwelch_window_tx, pwelch_noverlap_tx, nfft_tx, ...
                                SamplingFrequency, pwelchRange, pwelchSpectrumType);
     % Plot the PSD in dB scale, normalized to 0 dB peak.
     % Convert frequency axis to desired units (e.g., MHz).
     % Add 'eps' for numerical stability to avoid log10(0) or division by zero.
     plot(axPSD, freq_tx / freqUnitsFactor, 10*log10(psd_tx / max(psd_tx + eps) + eps), ...
          'Color', txColor, 'LineWidth', lineWidth, 'DisplayName', legendTextTx); % Assign DisplayName for legend

     % --- Calculate and Plot Received Signal PSD ---
     % Calculate PSD using Welch's method for the received signal.
     % Ensure signal_rx is a column vector.
     [psd_rx, freq_rx] = pwelch(signal_rx(:), pwelch_window_rx, pwelch_noverlap_rx, nfft_rx, ...
                                SamplingFrequency, pwelchRange, pwelchSpectrumType);
     % Plot the normalized PSD in dB scale for the received signal.
     plot(axPSD, freq_rx / freqUnitsFactor, 10*log10(psd_rx / max(psd_rx + eps) + eps), ...
          'Color', rxColor, 'LineWidth', lineWidth, 'DisplayName', legendTextRx); % Assign DisplayName for legend

     % --- Axes Formatting and Final Touches ---
     grid(axPSD, 'on');                     % Enable grid lines
     % Set axis labels, using the specified format for the x-axis label.
     xlabel(axPSD, sprintf(xAxisLabelFormat, freqUnitsLabel));
     ylabel(axPSD, yAxisLabel);
     title(axPSD, plotTitle);               % Set the main plot title
     % Add the legend with specified labels and location.
     legend(axPSD, 'show', 'Location', legendLocation);
     ylim(axPSD, yLimitDb);                 % Set y-axis limits (dB scale)
     % Set x-axis limits to cover the full spectrum [-Fs/2, Fs/2] in target units.
     xlim(axPSD, [-SamplingFrequency/2/freqUnitsFactor, SamplingFrequency/2/freqUnitsFactor]);
     box(axPSD, 'on');                      % Draw a box around the plot area
     hold(axPSD, 'off');                    % Release the plot hold
 end