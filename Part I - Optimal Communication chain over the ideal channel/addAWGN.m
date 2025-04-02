function signal_out = addAWGN(signal_in, Eb, EbN0dB, Fs)
    %   Adds Additive White Gaussian Noise to a complex baseband signal.
    %
    %   Inputs:
    %       signal_in   - Input complex baseband signal vector.
    %       Eb          - Energy per bit of the input signal (Ptx / BitRate).
    %       EbN0dB      - Desired Eb/N0 ratio in dB.
    %       Fs          - Sampling frequency of the input signal in Hz.
    %
    %   Output:
    %       signal_out  - Output signal vector (signal_in + complex noise).

        % --- Calculate Noise Parameters ---
        EbN0_linear = 10^(EbN0dB / 10); % Convert dB to linear
        N0 = Eb / EbN0_linear;          % One-sided Noise Power Spectral Density (W/Hz)
        noiseVariance = N0 * Fs;        % Variance (power) of complex baseband noise (sigma^2)

        % --- Generate Complex AWGN ---
        % Variance is split equally between real and imaginary parts (var/2 each)
        % Standard deviation for each component is sqrt(noiseVariance / 2)
        noise_real = sqrt(noiseVariance / 2) * randn(size(signal_in));
        noise_imag = sqrt(noiseVariance / 2) * randn(size(signal_in));
        noise = noise_real + 1i * noise_imag; % Complex noise vector

        % --- Add Noise ---
        signal_out = signal_in + noise;
end