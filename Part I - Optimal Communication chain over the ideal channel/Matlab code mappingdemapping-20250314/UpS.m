% OSF at least 2 (Nyquist)
function out = UpS(in, OSF)
    N = length(in);
    out = zeros(1,OSF*N);
    for i = 1:N
        out(OSF*i) = in(i);
    end
end