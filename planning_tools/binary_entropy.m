function h = binary_entropy(p)
%BINARY_ENTROPY Binary entropy function (entropy of Bernoulli process)

% assert(all(p(:) >= 0));
% assert(all(p(:) <= 1));

h = -(p .* log(p) + (1 - p) .* log(1 - p));
h(isnan(h)) = 0;

end
