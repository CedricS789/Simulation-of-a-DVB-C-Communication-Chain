function signal_out = addAWGN(signal_in, Eb, EbN0dB, Fs)
    %   Adds Additive White Gaussian Noise to a complex baseband signal.
    %   This function simulates the effect of a noisy communication channel
    %   by adding complex AWGN with a specified power level relative to the
    %   signal's energy per bit (Eb).
    %
    %   Inputs:
    %       signal_in   - Input complex baseband signal vector (transmitted signal).
    %                     Assumed to be a row or column vector.
    %       Eb          - Energy per bit of the input signal. This is typically
    %                     calculated as the average symbol energy divided by the
    %                     number of bits per symbol, or average transmit power
    %                     divided by the bit rate (Ptx / BitRate).
    %       EbN0dB      - Desired ratio of energy per bit (Eb) to noise power
    %                     spectral density (N0) in decibels (dB). This value
    %                     determines the signal-to-noise ratio per bit.
    %       Fs          - Sampling frequency of the input signal in Hz. This is
    %                     needed to determine the noise bandwidth and thus the
    %                     total noise power within the signal's bandwidth.
    %
    %   Output:
    %       signal_out  - Output signal vector, which is the input signal plus
    %                     the generated complex AWGN. It has the same dimensions
    %                     as signal_in.
    %
    %   Purpose:
    %       - Simulates the effect of random thermal noise encountered in real-world
    %         wireless or wired communication channels (AWGN = Additive White
    %         Gaussian Noise). AWGN is characterized by a constant power spectral
    %         density across all frequencies ('White') and a Gaussian amplitude
    %         distribution ('Gaussian'). It is added ('Additive') to the signal.
    %       - The strength of the added noise is controlled by the Eb/N0 ratio.
    %         A lower Eb/N0 value corresponds to a higher noise level relative
    %         to the signal energy per bit, simulating a poorer channel condition.
    %         Conversely, a higher Eb/N0 represents a cleaner channel.
    %       - This function is typically used in communication system simulations
    %         (e.g., in Step 3 as mentioned in the original comment) to evaluate
    %         the system's performance (like Bit Error Rate - BER) under various
    %         noise conditions (Section 4.3).

        % --- Calculate Noise Parameters ---
        % Convert the Eb/N0 ratio from decibels (dB) to a linear scale.
        % The formula for converting dB to linear is: linear = 10^(dB / 10).
        EbN0_linear = 10^(EbN0dB / 10);

        % Calculate the one-sided Noise Power Spectral Density (N0).
        % N0 represents the noise power per unit bandwidth (Watts/Hz).
        % It's derived from the definition Eb/N0: N0 = Eb / (Eb/N0)_linear.
        % A smaller N0 indicates less noise power spectral density.
        N0 = Eb / EbN0_linear;          % Noise PSD (W/Hz)

        % Calculate the total noise variance (average power) within the signal bandwidth.
        % For a baseband signal sampled at Fs, the Nyquist bandwidth is Fs/2.
        % However, for complex baseband AWGN, the total noise power across the
        % equivalent double-sided bandwidth Fs is N0 * Fs.
        % This variance represents the average power of the complex noise process n(t).
        % sigma^2 = E[|n(t)|^2] = E[n_real(t)^2] + E[n_imag(t)^2] = N0 * Fs.
        noiseVariance = N0 * Fs;        % Variance (power) of the complex baseband noise

        % --- Generate Complex AWGN ---
        % Generate complex Gaussian noise samples. The noise is complex because
        % the input signal is assumed to be complex baseband.
        % The total noise variance (noiseVariance) is split equally between the
        % real (in-phase) and imaginary (quadrature) components.
        % Therefore, the variance of the real part is noiseVariance / 2, and
        % the variance of the imaginary part is also noiseVariance / 2.
        % The standard deviation for each component is sqrt(variance / 2).
        % randn(size(signal_in)) generates standard Gaussian random numbers
        % (mean 0, variance 1) with the same dimensions as the input signal.
        % We scale these standard Gaussian numbers by the calculated standard deviation.

        % Generate the real part of the noise (In-phase component).
        noise_real = sqrt(noiseVariance / 2) * randn(size(signal_in));

        % Generate the imaginary part of the noise (Quadrature component).
        noise_imag = sqrt(noiseVariance / 2) * randn(size(signal_in));

        % Combine the real and imaginary parts to form the complex AWGN vector.
        noise = noise_real + 1i * noise_imag; % Complex Gaussian noise vector n = n_real + j*n_imag

        % --- Add Noise ---
        % Add the generated complex noise vector sample-by-sample to the
        % input signal vector to simulate the noisy channel effect.
        signal_out = signal_in + noise;
    end