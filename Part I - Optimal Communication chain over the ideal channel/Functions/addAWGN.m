function signal_out = addAWGN(signal_in, Eb, EbN0dB, OSF, SymRate)
    %   Adds Additive White Gaussian Noise to a complex baseband signal.
    %   Uses Eb based on the actual measured power of the input signal stream.
    %
    %   Inputs:
    %       signal_in   - Input complex baseband signal vector (pulse-shaped).
    %       Eb          - Calculated Energy per bit of the input signal stream
    %                     (Avg Signal Power / Bit Rate).
    %       EbN0dB      - Desired Eb/N0 ratio in dB.
    %       OSF         - Oversampling Factor.
    %       SymRate     - Symbol Rate in Hz.
    %
    %   Output:
    %       signal_out  - Output signal vector (signal_in + complex noise).

    % --- Calculate Noise Parameters based on provided Eb ---
    EbN0_linear = 10^(EbN0dB / 10);                         % Convert dB to linear SNR

    % N0 is the one-sided Power Spectral Density (PSD) of the noise in W/Hz.
    % Eb/N0 = (Signal Power / BitRate) / (Noise Power / Bandwidth)
    % For complex baseband AWGN, the total noise power (variance) sigma^2
    % in a bandwidth Fs is N0 * Fs.
    % N0 = Eb / EbN0_linear
    N0 = Eb / EbN0_linear;

    % Calculate the variance (average power) of the complex baseband noise
    % The noise bandwidth is the sampling frequency Fs = OSF * SymRate.
    noiseVariance = N0 * OSF * SymRate;                     % Variance (power) of the complex noise (sigma^2).
                                                            % This is the total variance for N = Nr + j*Ni
                                                            % Var(N) = E[|N|^2] = E[Nr^2] + E[Ni^2] = noiseVariance
                                                            % Assuming equal power in real/imag parts: Var(Nr) = Var(Ni) = noiseVariance / 2

    % --- Generate Complex AWGN ---
    % Generate real and imaginary parts separately, each with variance noiseVariance / 2.
    % Standard deviation for each component is sqrt(noiseVariance / 2).
    noise_real = sqrt(noiseVariance / 2) * randn(size(signal_in));
    noise_imag = sqrt(noiseVariance / 2) * randn(size(signal_in));

    % Combine into complex noise vector
    noise = noise_real + 1i * noise_imag;

    % --- Add Noise ---
    signal_out = signal_in + noise;  % Add generated noise to the input signal

end