function  f = dnorm(x,m,s,saturation,flip)
%DNORM distribution normale cumulée
%		ex. f = dnorm(x)
%		options : f = dnorm(x,m,s,saturation,flip)

% Woodox 3.00 - 24/04/01 - Olivier Vitrac

if nargin<2, m = 0; end
if nargin<3 || s==0, s = 1; end
if nargin<4, saturation = 0; end
if nargin<5, flip = 0; end
f = exp(-0.5*((x-m)./s).^2)./(sqrt(2*pi)*s);
if saturation
   [fmax,imax] = max(f);
   f = f/fmax;
   if flip, f(1:imax)=1;
   else f(imax:end) = 1; end
end