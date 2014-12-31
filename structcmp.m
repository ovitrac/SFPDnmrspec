function out = structcmp(s1,s2,fieldstoremove)
%COMPARE structures together (structure arrays are supported) whatever the order of fields
%   syntax: result = structcmp(s1,s2 [,fieldstoremove])
%           s1 and s2 structures to be compared
%           fieldstoremove: cell array of field names not to consider for  the comparison

% MS 2.0 - 09/04/11 - INRA\Olivier Vitrac - rev.

% Revision history

% arg heck
if nargin<2 || nargin>3, error('2 arguments are required'); end
if ~isstruct(s1), error('s1 must be a structure'); end
if ~isstruct(s2), error('s2 must be a structure'); end
if ~matcmp(size(s1),size(s2)), out = false; return; end
if nargin<3, fieldstoremove = {}; end
if ~iscell(fieldstoremove), fieldstoremove={fieldstoremove}; end

% remove fieldnames if asked
if ~isempty(fieldstoremove);
    for i=1:length(fieldstoremove);
        if isfield(s1,fieldstoremove{i}), s1 = rmfield(s1,fieldstoremove{i}); end
        if isfield(s2,fieldstoremove{i}), s2 = rmfield(s2,fieldstoremove{i}); end
    end
end

% compare fieldnames
f1 = fieldnames(s1);
f2 = fieldnames(s2);
if ~isempty(setdiff(f1,f2)), out = false; return; end
[~,i1] = sort(f1);
[~,i2] = sort(f2);

% compare sizes
if ~matcmp(size(s1),size(s2)), out = false; return; end

% compare values
for i = numel(s1)
    v1 = struct2cell(s1(i));
    v2 = struct2cell(s2(i));
    if ~cellcmp(v1(i1),v2(i2)); out = false; return; end
end
out = true;