% troncate conv at reception, not in filter
function h_t = RRC_filter(beta,BW,OSF,NbTaps) % Fs = OSF*BW, NbTaps is odd

    Fs = OSF*BW; % sampling freq
    step = (1/NbTaps)*Fs; % distance between freq points
    maxFreq = step*(NbTaps-1)/2;
    freqGrid = linspace(-maxFreq, maxFreq, NbTaps);

    T = 1/BW;   % symbol period
    lowerLimit = (1-beta)/(2*T);
    upperLimit = (1+beta)/(2*T);

    filterFreq = zeros(1,NbTaps);

    for i = 1:NbTaps
        if abs(freqGrid(i)) <= lowerLimit
            filterFreq(i) = T;
        elseif abs(freqGrid(i)) <= upperLimit
                filterFreq(i) = T/2 * (1 + cos((pi*T/beta)*(abs(freqGrid(i)) - lowerLimit)));
        end
    end

    h_t = fftshift(ifft(ifftshift(sqrt(filterFreq))));
    h_t = h_t/norm(h_t);

    Delta_t = 1/Fs;
    t = (-(NbTaps-1)/2:(NbTaps-1)/2)*Delta_t;

    figure
    subplot(1,2,1)
    plot(freqGrid, filterFreq)
    subplot(1,2,2)
    plot(t,h_t)
    grid on;
end