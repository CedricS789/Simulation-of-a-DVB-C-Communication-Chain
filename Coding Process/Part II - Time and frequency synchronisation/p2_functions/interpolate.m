function y = interpolate(y, x)
    x0 = floor(x);
    x1 = ceil(x);
    if x0 < 1 || x1 > length(y)
        y = 0;  % or NaN if you want to indicate invalid access
    elseif abs(x1 - x0) < eps
        y = y(x0);
    else
        alpha = (x - x0)/(x1-x0);
        y = (1 - alpha)*y(x0) + alpha*y(x1);
    end
end