function out = cellcmp(c1,c2)
%COMPARE cell arrays together (note: main classes have been implemented but not all)
%   syntax: result = cellcmp(c1,c2)

% MS 2.0 - 09/04/11 - INRA\Olivier Vitrac - rev.

% Revision history

% arg heck
if nargin~=2, error('2 arguments are required'); end
if ~iscell(c1), error('c1 must be a cell array'); end
if ~iscell(c2), error('c2 must be a cell array'); end
n = numel(c1);
if n~=numel(c2), out = false; return; end


% compare values
for i=1:numel(c1)
    if ~strcmp(class(c1{i}),class(c2{i})), out = false; return; end
    if ~matcmp(size(c1{i}),size(c2{i})), out = false; return; end
    if isnumeric(c1{i}) && ~matcmp(c1{i},c2{i}), out = false; return; end
    if iscellstr(c1{i}) && ~all(strcmp(c1{i},c2{i})), out = false; return; end
    if isstruct(c1{i}) && ~structcmp(c1{i},c2{i}), out = false; return; end
    if ischar(c1{i}) && ~strcmp(c1{i},c2{i}), out = false; return; end
end
out = true;