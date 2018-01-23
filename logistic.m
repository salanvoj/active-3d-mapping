function p = logistic(x)
%LOGISTIC Logistic function
%
% p = logistic(x) = 1 / (1 + exp(-x))
%
% Input:
% - x  Any size, -inf <= x <= inf.
%
% Output:
% - p  Same size as x, 1 / (1 + exp(-x)), 0 <= p <= 1.
%
% Note that x = logistic(logit(x)).
%

p = 1 ./ (1 + exp(-x));

end
