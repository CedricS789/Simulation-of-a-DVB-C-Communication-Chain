function impulseResponse = rrcFilter(Beta, SymRate, OSF, NumTaps)
    %   Generates the time-domain impulse response of a Root Raised Cosine (RRC) filter.
    %   Designed in frequency domain using sqrt(RC response) and transformed via IFFT.
    %   Result is normalized such that the convolution of the filter with itself
    %   results in a Raised Cosine pulse with a peak amplitude of 1 at t=0.
    %
    %   Inputs:
    %       Beta               - Roll-off factor (0 <= Beta <= 1).
    %       SymRate            - Symbol Rate (Rs) in Hz.
    %       OSF                - Oversampling Factor (Integer >= 1). Fs = OSF * SymRate.
    %       NumTaps            - Filter Length (odd number recommended for symmetry).
    %
    %   Output:
    %       impulseResponse    - Row vector of normalized RRC filter coefficients.

    % =====================================================================
    % == Parameter Calculations
    % =====================================================================
    samplingFreq = OSF * SymRate; % Sampling frequency (Fs)
    symbolPeriod = 1 / SymRate;   % Symbol period (Tsym) % Changed name for clarity

    % =====================================================================
    % == Frequency Domain Design ==
    % =====================================================================
    % Create frequency grid
    % Use NFFT slightly larger than NumTaps for better frequency resolution if needed,
    % but for direct IFFT approach, matching length is fine.
    NFFT = NumTaps; % Use filter length for FFT/IFFT
    freqGridHz = linspace(-samplingFreq/2, samplingFreq/2, NFFT); % Centered grid

    % Define characteristic frequencies for Raised Cosine (RC)
    f1_hz = (1 - Beta) / (2 * symbolPeriod); % End of passband
    f2_hz = (1 + Beta) / (2 * symbolPeriod); % Start of stopband

    % Construct ideal RC frequency response magnitude |H_RC(f)|
    raisedCosineFreqResponse = zeros(1, NFFT);
    absFreqHz = abs(freqGridHz);

    % Region 1: Passband (|f| <= f1)
    passbandIndices = (absFreqHz <= f1_hz);
    % RC pulse should have energy Tsym. We design sqrt(RC).
    % Let H_RC(f) scale with Tsym. Then H_RRC(f) scales with sqrt(Tsym)
    raisedCosineFreqResponse(passbandIndices) = symbolPeriod; % Target RC gain in passband

    % Region 2: Roll-off band (f1 < |f| <= f2)
    rolloffIndices = (absFreqHz > f1_hz) & (absFreqHz <= f2_hz);
    raisedCosineFreqResponse(rolloffIndices) = (symbolPeriod / 2) * (1 + cos((pi * symbolPeriod / Beta) * (absFreqHz(rolloffIndices) - f1_hz)));

    % Calculate RRC frequency response: H_RRC(f) = sqrt(H_RC(f))
    rrcFreqResponse = sqrt(raisedCosineFreqResponse);

    % =====================================================================
    % == Time Domain Transformation ==
    % =====================================================================
    % IFFT to get time-domain impulse response, centered using fftshift/ifftshift
    impulseResponseUnnormalized = fftshift(ifft(ifftshift(rrcFreqResponse)));

    % Ensure real output (due to potential numerical inaccuracies)
    impulseResponseUnnormalized = real(impulseResponseUnnormalized);

    % == Normalization for Matched Filtering ==
    % Normalize so that the combined pulse h_rc = conv(h_rrc, h_rrc) has peak amplitude 1.
    % This ensures the signal amplitude at the sampling instant is correctly scaled.
    h_rc_temp = conv(impulseResponseUnnormalized, impulseResponseUnnormalized);
    center_idx_rc = floor(length(h_rc_temp)/2) + 1; % Index of the peak for zero delay
    peak_val_rc = h_rc_temp(center_idx_rc);

    % Avoid division by zero or tiny numbers if filter design failed
    if peak_val_rc <= eps
        warning('RRC filter peak is near zero. Check parameters.');
        impulseResponse = impulseResponseUnnormalized; % Return unnormalized
    else
        % Scale factor is sqrt of the combined peak value
        scale_factor = sqrt(peak_val_rc);
        impulseResponse = impulseResponseUnnormalized / scale_factor;
    end


    % =====================================================================
    % == Plotting the Filter Characteristics ==
    % =====================================================================
    % samplePeriod = 1 / samplingFreq;
    % timeVectorSec = (-(NumTaps - 1) / 2 : (NumTaps - 1) / 2) * samplePeriod;

    % filterFig = figure;
    % set(filterFig, 'Name', sprintf('RRC Filter Characteristics (β=%.2f, Rs=%.1f MHz, OSF=%d, Taps=%d)', Beta, SymRate/1e6, OSF, NumTaps), 'NumberTitle', 'off');

    % % Plot 1: Frequency Domain Magnitude Response |H_RRC(f)|
    % axFreq = subplot(1, 2, 1);
    % plot(axFreq, freqGridHz / 1e6, abs(fftshift(fft(ifftshift(impulseResponse)))), 'b-', 'LineWidth', 1.5); % Plot actual H(f) of normalized h(t)
    % hold(axFreq, 'on');
    % plot(axFreq, [f1_hz, f1_hz]/1e6, ylim(axFreq), 'r--', 'LineWidth', 1.0, 'DisplayName', 'f_1'); % Mark f1
    % plot(axFreq, [-f1_hz, -f1_hz]/1e6, ylim(axFreq), 'r--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    % plot(axFreq, [f2_hz, f2_hz]/1e6, ylim(axFreq), 'g--', 'LineWidth', 1.0, 'DisplayName', 'f_2'); % Mark f2
    % plot(axFreq, [-f2_hz, -f2_hz]/1e6, ylim(axFreq), 'g--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    % hold(axFreq, 'off');
    % grid(axFreq, 'on');
    % xlabel(axFreq, 'Frequency (MHz)');
    % ylabel(axFreq, 'Magnitude |H_{RRC}(f)|');
    % title(axFreq, 'RRC Filter Frequency Response');
    % legend(axFreq, 'show', 'Location', 'best');
    % xlim(axFreq, samplingFreq/2*[-1 1] / 1e6); % Use samplingFreq/2
    % box(axFreq, 'on');

    % % Plot 2: Time Domain Impulse Response h_RRC(t)
    % axTime = subplot(1, 2, 2);
    % plot(axTime, timeVectorSec * 1e6, impulseResponse, 'b-', 'LineWidth', 1.5); % Use final normalized response
    % grid(axTime, 'on');
    % xlabel(axTime, 'Time (\mus)');
    % ylabel(axTime, 'Amplitude h_{RRC}(t)');
    % title(axTime, 'RRC Filter Impulse Response (Normalized)');
    % xlim(axTime, [min(timeVectorSec), max(timeVectorSec)] * 1e6);
    % box(axTime, 'on');
end