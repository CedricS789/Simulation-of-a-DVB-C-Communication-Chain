function [toa_hat, delta_cfo_hat] = frameFreqAcquisition(pilot, symb_rx_down, averaging_window, Tsymb)
    a = pilot;
    K = averaging_window;
    N = length(a);
    D_k = zeros(length(symb_rx_down), 1);
    D = zeros(length(symb_rx_down), K);
    y = symb_rx_down;
    sum_n_hat=zeros(length(symb_rx_down), 1);
    sum_delta_f = 0;
    
    for k = 1:K
        % Compute the differential correlation for each value of k
        for n = 1:(length(symb_rx_down) - (N-1))
            sum_D_k = 0;
            for l = k:N-1
                sum_D_k = sum_D_k + (conj(y(n+l))*a(l+1)) * conj((conj(y(n+l-k))*a(l-k+1)));
            end
            D_k(n) = (1/(N-k))*sum_D_k;
        end
        D(:, k) = D_k;
        sum_n_hat = sum_n_hat + abs(D_k);
    end
    
    [~, toa_hat] = max(sum_n_hat);
    
    for k = 1:K
        sum_delta_f = sum_delta_f + (angle(D(toa_hat, k))/(2*pi*k*Tsymb));
    end
    
    delta_cfo_hat = -(1/K)*sum_delta_f;
end
