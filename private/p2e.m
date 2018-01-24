function x = p2e(x)
%P2E
%
% x = p2e(x)

x = bsxfun(@rdivide, x(1:end-1, :), x(end, :));

end
