function X = ratio2fraction(R)
%RATIO2FRACTION founds x1+x2+...+xn=1 based on ratios of [x1/x2;x2/x3;...;xn-1/xn]
% SYNTAXES
%       X = ratio(R)
% INPUTS
%       R: n x 1 array of ratios (rank in order x1/x2; x2/x3; xn-1/xn)
% OUTPUTS
%       X: n+1 x 1 array of fractions 
% 
% EXAMPLE
% R = [1.5 2 5];
% X = ratio2fraction(R) 
% 
% F =
% 
%     0.4839
%     0.3226
%     0.1613
%     0.0323
%
% See also: nmrlsqnonneg, createcombination
%
% RMNSPEC v 0.1 - 10/12/2013 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev. 
%
% history
%
% argcheck
if nargin < 1, error('1 argument is required'), end
if isvector(R), R = R(:); end
if any(R==0) || any(R==Inf) || any(isnan(R)), error('Array R (ratios) can not contain 0, Inf ou NaN value'), end
% R(R==Inf) = 0; R(isnan(R)) = 0;
% main
A = diag(ones(1,length(R)+1)) + diag(-R,1); % matrix of passage containing ratios
A(end,:) = 1;
B = zeros(length(R)+1,1); B(end) = 1;
X =  A\B;
