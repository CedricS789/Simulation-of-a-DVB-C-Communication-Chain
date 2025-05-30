function h_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps)
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
    NFFT = NumTaps; % Use filter length for FFT/IFFT
    freqGridHz = linspace(-samplingFreq/2, samplingFreq/2, NFFT); % Centered grid

    % Define characteristic frequencies for Raised Cosine (RC)
    f1_hz = (1 - Beta) / (2 * symbolPeriod); % End of passband
    f2_hz = (1 + Beta) / (2 * symbolPeriod); % Start of stopband

    % Construct ideal RC frequency response magnitude |H_RC(f)|
    H_RC = zeros(1, NFFT);
    absFreqHz = abs(freqGridHz);

    % Region 1: Passband (|f| <= f1)
    passbandIndices = (absFreqHz <= f1_hz);
    H_RC(passbandIndices) = symbolPeriod; % Target RC gain in passband

    % Region 2: Roll-off band (f1 < |f| <= f2)
    rolloffIndices = (absFreqHz > f1_hz) & (absFreqHz <= f2_hz);
    H_RC(rolloffIndices) = (symbolPeriod / 2) * (1 + cos((pi * symbolPeriod / Beta) * (absFreqHz(rolloffIndices) - f1_hz)));

    % Calculate RRC frequency response: H_RRC(f) = sqrt(H_RC(f))
    H_RRC = sqrt(H_RC);


    % =====================================================================
    % == Time Domain Transformation ==
    % =====================================================================
    % IFFT to get time-domain impulse response, centered using fftshift/ifftshift
    impulseResponseUnnormalized = fftshift(ifft(ifftshift(H_RRC)));

    % Ensure real output (due to potential numerical inaccuracies)
    impulseResponseUnnormalized = real(impulseResponseUnnormalized)';

    % == Normalization for Matched Filtering ==
    % Normalize so that the combined pulse h_rc = conv(h_rrc, h_rrc) has peak amplitude 1.
    % This ensures the signal amplitude at the sampling instant is correctly scaled.
    h_rc_temp = conv(impulseResponseUnnormalized, impulseResponseUnnormalized);
    center_idx_rc = floor(length(h_rc_temp)/2) + 1; % Index of the peak for zero delay
    peak_val_rc = h_rc_temp(center_idx_rc);
    
    scale_factor = sqrt(peak_val_rc);
    h_rrc = impulseResponseUnnormalized / scale_factor;
    h_rc = h_rrc .* h_rrc;
















      % =====================================================================
      % == Plotting the Filter Characteristics ==
      % =====================================================================
      figure;
      set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14);
      plot(freqGridHz / 1e6, abs(fftshift(fft(ifftshift(h_rrc)))), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, 'Marker', '*', 'DisplayName', '$|H(f)|$'); % MATLAB blue
      hold on;
      plot([f1_hz, f1_hz]/1e6, ylim, 'r--', 'LineWidth', 1.0, 'DisplayName', '$f_1$');
      plot([-f1_hz, -f1_hz]/1e6, ylim, 'r--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
      plot([f2_hz, f2_hz]/1e6, ylim, 'g--', 'LineWidth', 1.0, 'DisplayName', '$f_2$');
      plot([-f2_hz, -f2_hz]/1e6, ylim, 'g--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
      hold off;
      grid on;
      xlabel('Frequency (MHz)', 'Interpreter', 'latex', 'FontSize', 16);
      ylabel('Magnitude $|H(f)|$', 'Interpreter', 'latex', 'FontSize', 16);
      title(['Raised Cosine Filter Discrete Frequency Points ($\beta$ = ', num2str(Beta), ', Taps = ', num2str(NumTaps), ')'], 'Interpreter', 'latex', 'FontSize', 18);
      legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 14);
      xlim(samplingFreq/2*[-1 1] / 1e6);
      box on;

      figure;
      set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14);
      samplePeriod = 1 / samplingFreq;
      timeVectorSec = (-(NumTaps - 1) / 2 : (NumTaps - 1) / 2) * samplePeriod;
      plot(timeVectorSec * 1e6, h_rc, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.5); % MATLAB purple
      grid on;
      xlabel('Time (\mu s)', 'FontSize', 16);
      ylabel('Amplitude h_{RRC}(t)', 'FontSize', 16);
      title(sprintf('RRC Filter Impulse Response (Normalized) - Number of Taps: %d', NumTaps), 'FontSize', 18); % Display NumTaps
      xlim([min(timeVectorSec), max(timeVectorSec)] * 1e6);
      box on;
end