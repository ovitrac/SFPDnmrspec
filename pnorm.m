function  p = pnorm(x,m,s,flip)
%PNORM distribution normale cumulée
%		ex. p = pnorm(x)
%		options : p = pnorm(x,m,s,flip)

% Woodox 3.00 - 05/04/01 - Olivier Vitrac - rev. 24/04/01

if nargin<2, m = 0; end
if nargin<3 || s==0, s = 1; end
if nargin<4, flip = 0; end
p = (1+erf((x-m)./(sqrt(2).*s)))/2;
if flip, p = 1-p; end
