function hFig = plotPSD_Tx_Rx(signal_tx, signal_rx, SamplingFrequency)
     %PLOTPSD_TX_RX Plots the Power Spectral Density (PSD) of transmitted and received signals.
     %   Uses Welch's method via 'pwelch'.
     %
     %   Inputs:
     %       signal_tx         - Time-domain transmitted signal vector.
     %       signal_rx         - Time-domain received signal vector.
     %       SamplingFrequency - Sampling frequency (Fs) in Hz.
     %
     %   Output:
     %       hFig              - Handle to the created figure object.

     % =====================================================================
     % == Plotting Parameters (Locally Defined) ==
     % =====================================================================
     txColor = [0, 0.4470, 0.7410];      % Blue
     rxColor = [0.8500, 0.3250, 0.0980]; % Red
     lineWidth = 1.5;
     figureName          = 'Power Spectral Density Comparison';
     figureNumberTitle   = 'off';
     plotTitle           = 'PSD Comparison: Transmitted vs. Received Signals';
     yAxisLabel          = 'Normalized Power Spectral Density (dB/Hz)';
     yLimitDb            = [-100, 5];    % Y-axis limits (dB) for viewing spectral shape
     freqUnitsFactor     = 1e6;          % Hz to MHz
     freqUnitsLabel      = 'MHz';
     xAxisLabelFormat    = 'Frequency (%s)';
     legendTextTx        = 'Transmitted Signal (Post-Pulse Shaping)';
     legendTextRx        = 'Received Signal (Post-Matched Filter)';
     legendLocation      = 'best';
     % --- pwelch Calculation Settings ---
     pwelchWindowType    = @hamming;     % Window function for Welch's method
     pwelchDefaultWindowLen = 512;       % Default segment length
     pwelchOverlapPercent= 50;           % Segment overlap percentage
     pwelchMinNfft       = 1024;         % Minimum FFT points
     pwelchSpectrumType  = 'psd';        % Estimate type ('psd' or 'power')
     pwelchRange         = 'centered';   % Frequency range ('centered' or 'onesided')

     % --- Figure and Axes Creation ---
     hFig = figure('Name', figureName, 'NumberTitle', figureNumberTitle);
     axPSD = axes('Parent', hFig);
     hold(axPSD, 'on');

     % --- Prepare pwelch Parameters ---
     % Determine window, overlap, NFFT for Tx signal
     pwelch_window_len_tx = min(length(signal_tx), pwelchDefaultWindowLen);
     pwelch_noverlap_tx = floor(pwelchOverlapPercent / 100 * pwelch_window_len_tx);
     pwelch_window_tx = pwelchWindowType(pwelch_window_len_tx);
     nfft_tx = max(pwelchMinNfft, 2^nextpow2(pwelch_window_len_tx));

     % Determine window, overlap, NFFT for Rx signal
     pwelch_window_len_rx = min(length(signal_rx), pwelchDefaultWindowLen);
     pwelch_noverlap_rx = floor(pwelchOverlapPercent / 100 * pwelch_window_len_rx);
     pwelch_window_rx = pwelchWindowType(pwelch_window_len_rx);
     nfft_rx = max(pwelchMinNfft, 2^nextpow2(pwelch_window_len_rx));

     % --- Calculate and Plot Transmitted Signal PSD ---
     [psd_tx, freq_tx] = pwelch(signal_tx(:), pwelch_window_tx, pwelch_noverlap_tx, nfft_tx, ...
                                SamplingFrequency, pwelchRange, pwelchSpectrumType);
     % Plot normalized PSD in dB
     plot(axPSD, freq_tx / freqUnitsFactor, 10*log10(psd_tx / max(psd_tx + eps) + eps), ...
          'Color', txColor, 'LineWidth', lineWidth, 'DisplayName', legendTextTx);

     % --- Calculate and Plot Received Signal PSD ---
     [psd_rx, freq_rx] = pwelch(signal_rx(:), pwelch_window_rx, pwelch_noverlap_rx, nfft_rx, ...
                                SamplingFrequency, pwelchRange, pwelchSpectrumType);
     % Plot normalized PSD in dB
     plot(axPSD, freq_rx / freqUnitsFactor, 10*log10(psd_rx / max(psd_rx + eps) + eps), ...
          'Color', rxColor, 'LineWidth', lineWidth, 'DisplayName', legendTextRx);

     % --- Axes Formatting and Final Touches ---
     grid(axPSD, 'on');
     xlabel(axPSD, sprintf(xAxisLabelFormat, freqUnitsLabel));
     ylabel(axPSD, yAxisLabel);
     title(axPSD, plotTitle);
     legend(axPSD, 'show', 'Location', legendLocation);
     ylim(axPSD, yLimitDb);
     xlim(axPSD, [-SamplingFrequency/2/freqUnitsFactor, SamplingFrequency/2/freqUnitsFactor]);
     box(axPSD, 'on');
     hold(axPSD, 'off');
 end