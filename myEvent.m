function [value, isterminal, direction] = myEvent(T, Y, thresh)
    value      = norm(Y) >= thresh;
    isterminal = 1;   % Stop the integration
    direction  = 0;
end