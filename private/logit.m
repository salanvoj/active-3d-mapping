function x = logit(p)
%LOGIT Logit function
%
% x = logit(p)
%
% Input:
% - p  Any size, 0 <= p <= 1.
%
% Output:
% - x  Same size as p, log(p / (1 - p)), -inf <= x <= inf.
%
% Note that x = logit(logistic(x)).
%

x = log(p ./ (1 - p));

end
