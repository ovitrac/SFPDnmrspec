function a = substructarray(s,f,dim)
% SUBSTRUCTARRAY do [struct.(field)] or {struct.(field)} or [struct.(field1) struct.(field2)...] or {struct.(field1) struct.(field2)...}
% Syntax: array = substructarray(struct,f [,dim])
%
%   See also STRUCTTAB2STRUCT, STRUCT2STRUCTTAB

% MS 2.1 - 24/12/12 - INRA\Olivier Vitrac - rev. 05/06/2014

% Revision history
% 05/03/2013 add [struct.(field1) struct.(field2)...]
% 05/06/2014 fix logical types, add see also

% arg check
if nargin<2, error('2 arguments are at least required'), end
if nargin<3, dim = 1; end
if ~isstruct(s), error('the first argument must be a structure'), end

% main()
if length(s)>1 && ischar(f)
    if ~isfield(s,f), error('''%s'' is not a valid field',f), end
    a = {s.(f)}; % default behavior (24/12/12)
elseif length(s)==1
    if ~iscell(f), f= {f}; end
    if ~iscellstr(f), error('f must be a cell array of strings'); end
    a = cellfun(@(g) s.(g),f,'UniformOutput',false);
else
    error('f must be a valid field name and not a list of field names')
end
if all(cellfun(@isnumeric,a))
    a = cat(dim,a{:});
elseif all(cellfun(@islogical,a))
    a = cat(dim,a{:});
else
    a = uncell(cat(dim,a(:)));
end