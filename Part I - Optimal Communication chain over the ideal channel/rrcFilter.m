function impulseResponse = rrcFilter(Beta, SymRate, OSF, NumTaps)
    %   Generates the time-domain impulse response of a Root Raised Cosine (RRC) filter.
    %   Designed in frequency domain using sqrt(RC response) and transformed via IFFT.
    %   Result is normalized to unit energy.
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
    symbolPeriod = 1 / SymRate;   % Symbol period (Ts)

    % =====================================================================
    % == Frequency Domain Design ==
    % =====================================================================
    % Create frequency grid
    freqStep = samplingFreq / NumTaps;
    maxGridFreq = freqStep * (NumTaps - 1) / 2;
    freqGridHz = linspace(-maxGridFreq, maxGridFreq, NumTaps);

    % Define characteristic frequencies for Raised Cosine (RC)
    f1_hz = (1 - Beta) / (2 * symbolPeriod); % End of passband
    f2_hz = (1 + Beta) / (2 * symbolPeriod); % Start of stopband

    % Construct ideal RC frequency response magnitude |H_RC(f)|
    raisedCosineFreqResponse = zeros(1, NumTaps);
    absFreqHz = abs(freqGridHz);

    % Region 1: Passband (|f| <= f1)
    passbandIndices = (absFreqHz <= f1_hz);
    raisedCosineFreqResponse(passbandIndices) = symbolPeriod;

    % Region 2: Roll-off band (f1 < |f| <= f2)
    rolloffIndices = (absFreqHz > f1_hz) & (absFreqHz <= f2_hz);
    raisedCosineFreqResponse(rolloffIndices) = (symbolPeriod / 2) * (1 + cos((pi * symbolPeriod / Beta) * (absFreqHz(rolloffIndices) - f1_hz)));

    % Calculate RRC frequency response: H_RRC(f) = sqrt(H_RC(f))
    rrcFreqResponse = sqrt(raisedCosineFreqResponse);

    % =====================================================================
    % == Time Domain Transformation ==
    % =====================================================================
    % IFFT to get time-domain impulse response, centered at t=0
    impulseResponseUnnormalized = fftshift(ifft(ifftshift(rrcFreqResponse)));

    % Normalize the impulse response to have unit energy (L2 norm = 1)
    impulseResponse = impulseResponseUnnormalized / norm(impulseResponseUnnormalized);

    % == Optional: Plotting within function (kept as per original structure) ==
    % NOTE: For a minimal *library* function, this plotting might be removed.
    % Keeping it as instructed to only change comments.

    % =====================================================================
    % == Time Vector for Plotting ==
    % =====================================================================
    samplePeriod = 1 / samplingFreq;
    timeVectorSec = (-(NumTaps - 1) / 2 : (NumTaps - 1) / 2) * samplePeriod;

    % =====================================================================
    % == Plotting the Filter Characteristics ==
    % =====================================================================
    filterFig = figure;
    set(filterFig, 'Name', sprintf('RRC Filter Characteristics (β=%.2f, Rs=%.1f MHz, OSF=%d, Taps=%d)', Beta, SymRate/1e6, OSF, NumTaps), 'NumberTitle', 'off');

    % Plot 1: Frequency Domain Magnitude Response |H_RRC(f)|
    axFreq = subplot(1, 2, 1);
    plot(axFreq, freqGridHz / 1e6, abs(rrcFreqResponse), 'b-', 'LineWidth', 1.5);
    hold(axFreq, 'on');
    plot(axFreq, [f1_hz, f1_hz]/1e6, ylim(axFreq), 'r--', 'LineWidth', 1.0, 'DisplayName', 'f_1'); % Mark f1
    plot(axFreq, [-f1_hz, -f1_hz]/1e6, ylim(axFreq), 'r--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    plot(axFreq, [f2_hz, f2_hz]/1e6, ylim(axFreq), 'g--', 'LineWidth', 1.0, 'DisplayName', 'f_2'); % Mark f2
    plot(axFreq, [-f2_hz, -f2_hz]/1e6, ylim(axFreq), 'g--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    hold(axFreq, 'off');
    grid(axFreq, 'on');
    xlabel(axFreq, 'Frequency (MHz)');
    ylabel(axFreq, 'Magnitude |H_{RRC}(f)|');
    title(axFreq, 'RRC Filter Frequency Response');
    legend(axFreq, 'show', 'Location', 'best');
    xlim(axFreq, [min(freqGridHz), max(freqGridHz)] / 1e6);
    box(axFreq, 'on');

    % Plot 2: Time Domain Impulse Response h_RRC(t)
    axTime = subplot(1, 2, 2);
    plot(axTime, timeVectorSec * 1e6, real(impulseResponse), 'b-', 'LineWidth', 1.5); % Plot real part
    grid(axTime, 'on');
    xlabel(axTime, 'Time (\mus)');
    ylabel(axTime, 'Amplitude h_{RRC}(t)');
    title(axTime, 'RRC Filter Impulse Response');
    xlim(axTime, [min(timeVectorSec), max(timeVectorSec)] * 1e6);
    box(axTime, 'on');
end