function signal_out = addAWGN(signal_in, Eb, EbN0dB, OSF, SymRate)
    EbN0_linear = 10^(EbN0dB / 10); 
    N0 = Eb / EbN0_linear;
    noiseVariance = N0 * OSF * SymRate;
    noise_real = sqrt(noiseVariance / 2) * randn(size(signal_in));
    noise_imag = sqrt(noiseVariance / 2) * randn(size(signal_in));
    noise = noise_real + 1i * noise_imag;
    signal_out = signal_in + noise;
end