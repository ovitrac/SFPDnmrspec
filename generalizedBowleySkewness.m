function skewness = generalizedBowleySkewness(x,y,ythresh,p)
%GENERALIZEDBOWLEYSKEWNESS calculate the generalized quartile skewness coefficient around zero
%   q = generalizedBowleySkewness(x,y [ythresh,p])
%       x,y = mx1 vectors
%       ythresh use only values y>thresh (default = 0.05)
%       p: percentile value to consider (default = 0.1);
%       q: 1x2 skewness definitions
%           q(1)=<xn>  distance skewness (symmetry)
%           q(2)=<xn^3> skewness
%           where xn is the standardized x and <> is the average weighted by y
%
%   source: http://mathworld.wolfram.com/BowleySkewness.html (not robust and finally removed)
%           http://en.wikipedia.org/wiki/Skewness           
%
%   Simple examples: cos(x) gives -2e-4 and sin(x) gives -0.4
%       x = linspace(-2*pi,2*pi,1e4);
%       generalizedBowleySkewness(x,cos(x))
%       generalizedBowleySkewness(x,sin(x))
%
%   Generalized example
%   x = linspace(-2*pi,2*pi,1e5)'; y=(cos(x+pi/3)+cos(4*x+pi)).*exp(-0.5*abs(x));
%   generalizedBowleySkewness(x,y), clf, plot(x,[y cumtrapz(x,abs(y))])


% RMNSPEC v 0.1 - 30/01/2014 - INRA\Olivier Vitrac - rev. 01/02/2014

% Revision history
% 01/02/2014 standard implementation, fix examples

% Default
ythresh_default = 0.05;
p_default = 0.1;


% check arguments
if nargin<2, error('2 arguments are required'), end
if nargin<3, ythresh =[]; end
if nargin<4, p =[]; end
if isempty(ythresh), ythresh = ythresh_default; end
if isempty(p), p = p_default; end
x = x(:); y = y(:); m = length(x);
if length(y)~=m, error('x and y must be vectors of same size'); end

% do calculations
ymin = min(y); ymax = max(y);
y = (y - ymin)/(ymax-ymin); 
y(y<ythresh)=0;
ycum = cumtrapz(x,y);
L = interpleft(ycum/ycum(end),x,[p 1-p]);
central = (x>=L(1))&(x<=L(2));
x = x(central);
y = y(central);
y = y/trapz(x,y);
s = sqrt(trapz(x,y.*x.^2));
m = trapz(x,y.*x);
m3 = trapz(x,y.*x.^3);

skewness = [m/s m3/s^3];