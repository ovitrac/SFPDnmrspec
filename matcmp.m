function out = matcmp(m1,m2,recursion)
%COMPARE numerical arrays together (with the rule NaN==NaN is true)
%   syntax: result = matcmp(m1,m2)

% MS 2.0 - 09/04/11 - INRA\Olivier Vitrac - rev.

% Revision history

% arg heck
if nargin<2 || nargin>3, error('2 arguments are required'); end
if nargin<3, recursion = true; end
if ~isnumeric(m1), error('m1 must be a numeric array'); end
if ~isnumeric(m2), error('m2 must be a numeric array'); end
n = numel(m1);
if n~=numel(m2), out = false; return; end
if recursion
    if ~matcmp(size(m1),size(m2),false), out = false; return; end
end


% compare NaN values
inan1 = isnan(m1);
inan2 = isnan(m2);
if ~all(inan1==inan2), out = false; return; end

% compare ~NaN values
out = all(m1(~inan1)==m2(~inan2));