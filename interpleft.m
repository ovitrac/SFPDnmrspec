function yi=interpleft(x,y,xi)
%INTERPLEFT interpolates by taking the nearest value on the left
% yi=interpleft(x,y,xi)

% MS-MATLAB-WEB 1.0 - 30/04/09 - Olivier Vitrac - rev.

n = length(x);
if numel(y)~=n, error('INTERPLEFFT x and y must vectors of same length'), end
if n==1
    yi=y;
else
    if any(diff(x)==0), [x,u] = unique(x); y = y(u); n=length(u); end
    yi=y(max(1,min(n,fix(interp1(x,1:n,xi,'linear','extrap')))));
end