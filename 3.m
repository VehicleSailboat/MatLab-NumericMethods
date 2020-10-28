function [res N] = sit(f,x0, eps)
syms x;
k=(subs(f,x0));
xp = x0; n=0;
while abs(k - xp) > eps
xp = k;
k = (subs(f,k));
n=n+1;
end
res = sym2poly(k);
N=n;
end