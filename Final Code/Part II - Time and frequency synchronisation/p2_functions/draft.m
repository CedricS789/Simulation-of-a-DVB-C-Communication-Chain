% a = pilot;
% K = 8;
% N = length(a);
% D = zeros(length(K),length(symb_rx));
% sum_n_hat=zeros(1,length(symb_rx));
% y = [symb_rx.' zeros(1,N)];
% 
% for k = 1:K
%      % Compute the differential correlation for each value of k
%     for n = 1:(length(symb_rx) - (N-1))
%         sum_D_k = 0;
%         for l = k:N-1
%             sum_D_k = sum_D_k + (conj(y(n+l))*a(l+1)) * conj((conj(y(n+l-k))*a(l-k+1)));
%         end
%             D(k,n) = (1/(N-k))*sum_D_k;
%     end
%     sum_n_hat = sum_n_hat + abs(D(k,:));
% end
% 
% [~, toa] = max(sum_n_hat);
% 
% sum_delta_f = 0;
% for k = 1:K
%     sum_delta_f = sum_delta_f + (angle(D(k,toa))/(2*pi*k*Tsymb));
% end
% cfo = -(1/K)*sum_delta_f;