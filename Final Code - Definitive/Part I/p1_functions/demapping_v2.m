function [bit_rx] = demapping_v2(symb_rx_down, Nbps, modulation)
    symb_rx_down = symb_rx_down(:);
    Nsymb = size(symb_rx_down, 1); 

    if Nsymb == 0
        bit_rx = [];
        return;
    end

    switch lower(modulation) 

        case 'pam'
            M = 2^Nbps; 

            pam_levels = 0:(M-1);
            mean_level = mean(pam_levels); 
            sigma = sqrt(mean((pam_levels - mean_level).^2));
            
            int_rx = sigma * real(symb_rx_down) + mean_level;

            int_det = round(int_rx);

            int_det(int_det < 0) = 0;
            int_det(int_det > M-1) = M-1;

            try
                mapp_rx = de2bi(int_det, Nbps, 'left-msb');
            catch ME 
                 rethrow(ME); 
            end

            if size(mapp_rx, 1) ~= Nsymb || size(mapp_rx, 2) ~= Nbps
                 if Nsymb == 0 && isempty(mapp_rx)
                     mapp_rx = zeros(0, Nbps); 
                 elseif size(mapp_rx, 2) < Nbps 
                    mapp_rx = [zeros(Nsymb, Nbps - size(mapp_rx, 2)), mapp_rx];
                 elseif size(mapp_rx, 1) ~= Nsymb
                 end
            end

            bit_rx2 = zeros(Nsymb, Nbps); 
            bit_rx2(:, 1) = mapp_rx(:, 1); 
            for ii = 2:Nbps
                bit_rx2(:, ii) = xor(mapp_rx(:, ii - 1), mapp_rx(:, ii));
            end

            bit_rx = reshape(bit_rx2', Nsymb * Nbps, []); 

        case 'qam'
            
            Nbps_per_dim = Nbps / 2; 
            M_per_dim = 2^Nbps_per_dim; 

            pam_levels_1D = 0:(M_per_dim - 1);
            mean_level_1D = mean(pam_levels_1D); 
            sigma_1D = sqrt(mean((pam_levels_1D - mean_level_1D).^2));
            
            symb_rxI = real(symb_rx_down); 

            int_rxI = sigma_1D * sqrt(2) * symb_rxI + mean_level_1D;

            int_detI = round(int_rxI);
            int_detI(int_detI < 0) = 0;
            int_detI(int_detI > M_per_dim - 1) = M_per_dim - 1;

            try
                mapp_rxI = de2bi(int_detI, Nbps_per_dim, 'left-msb');
            catch ME
                 rethrow(ME);
            end

            if size(mapp_rxI, 1) ~= Nsymb || size(mapp_rxI, 2) ~= Nbps_per_dim
                 if Nsymb == 0 && isempty(mapp_rxI)
                     mapp_rxI = zeros(0, Nbps_per_dim);
                 elseif size(mapp_rxI, 2) < Nbps_per_dim
                    mapp_rxI = [zeros(Nsymb, Nbps_per_dim - size(mapp_rxI, 2)), mapp_rxI];
                 elseif size(mapp_rxI, 1) ~= Nsymb
                 end
            end

            bit_rx2I = zeros(Nsymb, Nbps_per_dim);
            bit_rx2I(:, 1) = mapp_rxI(:, 1);
            for ii = 2:Nbps_per_dim
                bit_rx2I(:, ii) = xor(mapp_rxI(:, ii - 1), mapp_rxI(:, ii));
            end

            symb_rxQ = imag(symb_rx_down); 

            int_rxQ = sigma_1D * sqrt(2) * symb_rxQ + mean_level_1D;

            int_detQ = round(int_rxQ);
            int_detQ(int_detQ < 0) = 0;
            int_detQ(int_detQ > M_per_dim - 1) = M_per_dim - 1;

            try
                mapp_rxQ = de2bi(int_detQ, Nbps_per_dim, 'left-msb');
            catch ME
                 rethrow(ME);
            end

            if size(mapp_rxQ, 1) ~= Nsymb || size(mapp_rxQ, 2) ~= Nbps_per_dim
                 if Nsymb == 0 && isempty(mapp_rxQ)
                     mapp_rxQ = zeros(0, Nbps_per_dim);
                 elseif size(mapp_rxQ, 2) < Nbps_per_dim
                    mapp_rxQ = [zeros(Nsymb, Nbps_per_dim - size(mapp_rxQ, 2)), mapp_rxQ];
                 elseif size(mapp_rxQ, 1) ~= Nsymb
                 end
            end

            bit_rx2Q = zeros(Nsymb, Nbps_per_dim);
            bit_rx2Q(:, 1) = mapp_rxQ(:, 1);
            for ii = 2:Nbps_per_dim
                bit_rx2Q(:, ii) = xor(mapp_rxQ(:, ii - 1), mapp_rxQ(:, ii));
            end

            bit_rx_combined = [bit_rx2I, bit_rx2Q]; 
            bit_rx = reshape(bit_rx_combined', Nsymb * Nbps, []); 
    end

    bit_rx = bit_rx(:);

end