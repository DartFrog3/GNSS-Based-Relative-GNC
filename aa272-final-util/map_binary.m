function x_pm = map_binary(x)
    % binary mapping as in hw

    x_pm = ones(size(x), 'double');
    x_pm(x == 1) = -1;
end
