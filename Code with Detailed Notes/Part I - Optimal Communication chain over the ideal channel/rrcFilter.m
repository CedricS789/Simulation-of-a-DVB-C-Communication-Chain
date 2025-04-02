function impulseResponse = rrcFilter(Beta, SymRate, OSF, NumTaps)
    %   Generates the time-domain impulse response of a Root Raised Cosine (RRC)
    %   filter. This type of filter is commonly used in digital communication
    %   systems for pulse shaping (at the transmitter) and matched filtering
    %   (at the receiver) to minimize Inter-Symbol Interference (ISI) while
    %   limiting the signal bandwidth.
    %
    %   The filter is designed by first defining the ideal Raised Cosine (RC)
    %   frequency response, taking its square root to get the RRC frequency
    %   response, and then transforming this RRC response to the time domain
    %   using the Inverse Fast Fourier Transform (IFFT). The resulting impulse
    %   response is normalized to have unit energy (L2 norm = 1).
    %
    %   Inputs:
    %       Beta               - Roll-off factor (0 <= Beta <= 1). This parameter
    %                            controls the excess bandwidth of the filter beyond
    %                            the Nyquist minimum bandwidth (SymRate / 2).
    %                            Beta = 0 corresponds to an ideal rectangular filter (brick wall).
    %                            Beta = 1 corresponds to a full cosine roll-off, using
    %                            double the minimum bandwidth.
    %       SymRate            - Symbol Rate (Rs or 1/Ts) in symbols per second (Hz).
    %                            This determines the rate at which symbols are transmitted.
    %       OSF                - Oversampling Factor (Integer >= 1). The number of samples
    %                            used to represent each symbol period in the discrete-time
    %                            signal. It determines the sampling frequency: Fs = OSF * SymRate.
    %       NumTaps            - Number of Taps (Filter Length). The desired length of the
    %                            filter's impulse response in samples. Should ideally be
    %                            an odd number to ensure a symmetric impulse response (linear phase).
    %                            A larger number of taps generally leads to a better approximation
    %                            of the ideal filter characteristics but increases computational complexity.
    %
    %   Output:
    %       impulseResponse    - A row vector containing the normalized time-domain
    %                            impulse response coefficients of the RRC filter.
    %                            The length of the vector is NumTaps.
    %
    %   Note: The design approach here uses frequency sampling. It defines the
    %         target RRC frequency response on a grid and then uses IFFT.
    %         The function also generates plots visualizing the designed filter's
    %         frequency response and time-domain impulse response for verification.

    % =====================================================================
    % == Parameter Calculations
    % =====================================================================
    % Calculate the sampling frequency (Fs) based on the symbol rate and oversampling factor.
    % Fs represents the rate at which the continuous-time signal is sampled.
    samplingFreq = OSF * SymRate; % Sampling frequency (Fs) in Hz

    % Calculate the symbol period (Ts) which is the duration of one symbol.
    symbolPeriod = 1 / SymRate;   % Symbol period (Ts) in seconds

    % =====================================================================
    % == Frequency Domain Design ==
    % =====================================================================
    % Create a frequency grid (vector) for designing the filter's frequency response.
    % The grid is centered around 0 Hz and spans symmetrically up to frequencies
    % relevant for the given sampling frequency and number of taps.
    % The frequency resolution (step size) depends on Fs and NumTaps.
    freqStep = samplingFreq / NumTaps; % Frequency resolution (Hz)
    % Calculate the maximum frequency extent of the grid.
    maxGridFreq = freqStep * (NumTaps - 1) / 2; % Max frequency in the symmetric grid (Hz)
    % Generate the linearly spaced frequency grid vector.
    freqGridHz = linspace(-maxGridFreq, maxGridFreq, NumTaps); % Symmetrical frequency grid (Hz)

    % Define the characteristic frequencies that determine the shape of the Raised Cosine spectrum.
    % f1: End of the flat passband region.
    f1_hz = (1 - Beta) / (2 * symbolPeriod); % Lower transition frequency (Hz)
    % f2: Start of the zero stopband region (end of the roll-off region).
    f2_hz = (1 + Beta) / (2 * symbolPeriod); % Upper transition frequency (Hz)

    % Initialize the frequency response array for the full Raised Cosine (RC) filter.
    % We design the RC response first, then take the square root for the RRC response.
    raisedCosineFreqResponse = zeros(1, NumTaps);

    % Construct the ideal Raised Cosine frequency response magnitude |H_RC(f)|.
    % This is done piecewise based on the frequency regions defined by f1 and f2.
    % Using vectorized operations for efficiency.
    absFreqHz = abs(freqGridHz); % Calculate absolute frequencies for easier comparison

    % Region 1: Passband (|f| <= f1)
    % The response is flat and equal to the symbol period (Ts).
    passbandIndices = (absFreqHz <= f1_hz);
    raisedCosineFreqResponse(passbandIndices) = symbolPeriod;

    % Region 2: Roll-off band (f1 < |f| <= f2)
    % The response follows a cosine shape, transitioning smoothly from Ts to 0.
    rolloffIndices = (absFreqHz > f1_hz) & (absFreqHz <= f2_hz);
    % Apply the standard raised cosine formula for the roll-off region.
    raisedCosineFreqResponse(rolloffIndices) = (symbolPeriod / 2) * (1 + cos((pi * symbolPeriod / Beta) * (absFreqHz(rolloffIndices) - f1_hz)));

    % Region 3: Stopband (|f| > f2)
    % The response is zero in the stopband. This is already handled by the
    % initialization with zeros.

    % Calculate the Root Raised Cosine (RRC) frequency response |H_RRC(f)|.
    % The RRC filter's frequency response is the square root of the RC filter's response.
    % H_RRC(f) = sqrt(H_RC(f)). When two RRC filters are cascaded (Tx and Rx),
    % their combined response is H_RRC(f) * H_RRC(f) = H_RC(f), which satisfies
    % Nyquist's criterion for zero ISI (at symbol sampling instants).
    rrcFreqResponse = sqrt(raisedCosineFreqResponse);

    % =====================================================================
    % == Time Domain Transformation ==
    % =====================================================================
    % Obtain the time-domain impulse response h_RRC(t) by taking the Inverse
    % Fast Fourier Transform (IFFT) of the designed RRC frequency response H_RRC(f).

    % Apply ifftshift before IFFT: The rrcFreqResponse vector was created corresponding
    % to a frequency grid centered at 0 Hz (-maxGridFreq to +maxGridFreq). `ifftshift`
    % rearranges this into the standard FFT order (0 Hz to Fs/2, then -Fs/2 to 0).
    % Apply fftshift after IFFT: The result of `ifft` is typically ordered from time 0
    % onwards. `fftshift` rearranges the time-domain impulse response so that the
    % center tap (corresponding to t=0) is in the middle of the vector.
    impulseResponseUnnormalized = fftshift(ifft(ifftshift(rrcFreqResponse)));

    % Normalize the impulse response to have unit energy (L2 norm = 1).
    % norm(vector) calculates the Euclidean length (sqrt(sum(|elements|^2))).
    % Normalization ensures that the filter doesn't change the signal power and
    % is standard practice for matched filters.
    impulseResponse = impulseResponseUnnormalized / norm(impulseResponseUnnormalized);

    % =====================================================================
    % == Time Vector for Plotting ==
    % =====================================================================
    % Create a time vector corresponding to the filter taps (samples) of the
    % impulse response. The vector should be centered around t=0.
    samplePeriod = 1 / samplingFreq; % Duration of one sample interval (seconds)
    % Generate time instances for each tap, from -(NumTaps-1)/2 * Tsamp to +(NumTaps-1)/2 * Tsamp.
    timeVectorSec = (-(NumTaps - 1) / 2 : (NumTaps - 1) / 2) * samplePeriod;

    % =====================================================================
    % == Plotting the Filter Characteristics ==
    % =====================================================================
    % Create a new figure to display the filter properties.
    filterFig = figure;
    % Set the figure name for clarity, including key parameters.
    set(filterFig, 'Name', sprintf('RRC Filter Characteristics (β=%.2f, Rs=%.1f MHz, OSF=%d, Taps=%d)', Beta, SymRate/1e6, OSF, NumTaps), 'NumberTitle', 'off');

    % Plot 1: Frequency Domain Magnitude Response |H_RRC(f)|
    axFreq = subplot(1, 2, 1); % Create subplot for frequency response
    % Plot the magnitude of the designed RRC frequency response vs frequency in MHz.
    plot(axFreq, freqGridHz / 1e6, abs(rrcFreqResponse), 'b-', 'LineWidth', 1.5);
    hold(axFreq, 'on'); % Hold the plot to add reference lines
    % Add vertical dashed lines indicating the key transition frequencies f1 and f2 (in MHz).
    plot(axFreq, [f1_hz, f1_hz]/1e6, ylim(axFreq), 'r--', 'LineWidth', 1.0, 'DisplayName', 'f_1 = (1-β)/2T_s');
    plot(axFreq, [-f1_hz, -f1_hz]/1e6, ylim(axFreq), 'r--', 'LineWidth', 1.0, 'HandleVisibility', 'off'); % Negative f1
    plot(axFreq, [f2_hz, f2_hz]/1e6, ylim(axFreq), 'g--', 'LineWidth', 1.0, 'DisplayName', 'f_2 = (1+β)/2T_s');
    plot(axFreq, [-f2_hz, -f2_hz]/1e6, ylim(axFreq), 'g--', 'LineWidth', 1.0, 'HandleVisibility', 'off'); % Negative f2
    hold(axFreq, 'off'); % Release the plot hold
    grid(axFreq, 'on'); % Add grid lines
    xlabel(axFreq, 'Frequency (MHz)'); % Set x-axis label
    ylabel(axFreq, 'Magnitude |H_{RRC}(f)|'); % Set y-axis label
    title(axFreq, 'RRC Filter Frequency Response'); % Set title
    legend(axFreq, 'show', 'Location', 'best'); % Show legend
    xlim(axFreq, [min(freqGridHz), max(freqGridHz)] / 1e6); % Set x-axis limits based on the grid
    box(axFreq, 'on'); % Draw a box around the plot

    % Plot 2: Time Domain Impulse Response h_RRC(t)
    axTime = subplot(1, 2, 2); % Create subplot for time response
    % Plot the real part of the normalized impulse response vs time in microseconds.
    % Due to numerical precision, there might be a tiny imaginary component, so we plot the real part.
    plot(axTime, timeVectorSec * 1e6, real(impulseResponse), 'b-', 'LineWidth', 1.5);
    grid(axTime, 'on'); % Add grid lines
    xlabel(axTime, 'Time (\mus)'); % Set x-axis label
    ylabel(axTime, 'Amplitude h_{RRC}(t)'); % Set y-axis label
    title(axTime, 'RRC Filter Impulse Response'); % Set title
    xlim(axTime, [min(timeVectorSec), max(timeVectorSec)] * 1e6); % Set x-axis limits based on the time vector
    box(axTime, 'on'); % Draw a box around the plot
end