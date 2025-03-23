function out = DownS(in, OSF)
    N = floor(length(in)/OSF);
    out = zeros(1,N);
    for i = 1:N
        out(i) = in(i*OSF);
    end
end