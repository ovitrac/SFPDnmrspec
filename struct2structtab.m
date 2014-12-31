function t=struct2structtab(s)
%STRUCT2STRUCTTAB cast a structure into a structure array
%      t = struct2structtab(s)
%
%   See also STRUCTTAB2STRUCT, SUBSTRUCTARRAY


% MS 2.1 - 20/01/12 - INRA\Olivier Vitrac rev.  05/06/2014

% Revision history
% 20/03/2013 fix mixed column-wise and row-wise arrays and cell arrays
% 14/10/2013 fix logical fields
% 15/10/2013 fix empty fields
% 05/06/2014 add see also

% argcheck
if nargin<1 || nargin>1, error('one argument is required'), end
if ~isstruct(s), error('the argument must be a structure'), end
if isempty(s) || numel(s)>1, error('structure arrays or empty structure are not authorized'), end
siz = structfun(@numel,s); sizmax = max(siz);
if any(diff(siz(siz>0))), error('all fields must be of a same size'), end
f = fieldnames(s);

% convert numeric fields and force column wise
% fn = f(structfun(@isnumeric,s));
% for i=1:length(fn)
%     s.(fn{i}) = num2cell(s.(fn{i})(:),2);
% end
for eachf = f'
    if isempty(s.(eachf{1}))
        s.(eachf{1}) = repmat({s.(eachf{1})},sizmax,1);
    elseif isnumeric(s.(eachf{1})) || islogical(s.(eachf{1}))
        s.(eachf{1}) = num2cell(s.(eachf{1})(:),2);
    else
        s.(eachf{1}) = s.(eachf{1})(:);
    end
end
v = struct2cell(s);

% recast
tmp = [f v]';
t = struct(tmp{:});