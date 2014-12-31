function y = removesymmetryinstrings(x)
%REMOVESYMMETRYINSTRINGS remove chains that are symmetric in a list
%   y = removesymmetryinstrings(x)
%
%   Example
%       removesymmetryinstrings({'123' 'EEEEAAAEE' '95647'})
%       gives
%           '123'    'EEAAAEEEE'    '74659'

% MS 2.1 -  28/04/2014 - INRA\Olivier Vitrac -rev. 

% Revision history
% 28/04/2014 RC, flip is forced if half is crossed (it was the reverse in RC)

% arg check
if ~nargin, error('one argument is required'), end
if ~ischar(x) && ~iscellstr(x), error('x must be a string or a string array'), end
if iscellstr(x) % recursion if needed
    y = cell(size(x));
    for j=1:numel(x)
        y{j} = removesymmetryinstrings(x{j});
    end
    return
end

% do calculations
nx   = length(x);
imax = nx/2;
completed = false;
flip      = true;
i    = 0;
while ~completed && i<imax
    i = i + 1;
    if x(i)==x(end-i+1)
        completed = false;
    elseif x(i)<x(end-i+1)
        completed = true;
        flip   = false;
    else
        completed = true;
    end
end
if flip, y = fliplr(x); else y=x; end