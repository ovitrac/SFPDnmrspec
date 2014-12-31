function s = structtab2struct(sarray,dim)
%STRUCTTAB2STRUCT converts a struct array into a structure with array/cell field contents (reverse operation of struct2stucttab)
%   syntax: s = structtab2struct(sarray [,dim])
%   example: structtab2struct(dir)
%
%   See also STRUCTTAB2STRUCT, SUBSTRUCTARRAY


% MS 2.1 - 29/11/13 - INRA\Olivier Vitrac - rev. 05/06/2014

% Revision history
%  05/06/2014 add see also

if nargin<1, error('one argument is required'), end
if ~isstruct(sarray), error('the first argument must be a structure'), end
if nargin<2, dim = 1; end

fn = fieldnames(sarray);
s = cell2struct(cellfun(@(f) substructarray(sarray,f,dim),fn,'UniformOutput',false),fn);