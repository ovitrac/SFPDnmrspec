function [listex,iout,jout]=expandtextaslist(list,splitstr)
%expandtextaslist expands cell array of strings containing semi-colon separated lists e.g. {'element1;element2;element3' 'txt1;txt2;txt3'}
%    Syntax: listex=expandtextaslist(list [,splitstr])
%   Options: [listex,i,j]=expandtextaslist(...)
% INPUTS
%      list: string or mx1 or 1xm cell array (default separator = ';' or '|', multiple separators are discarded, extra spaces are trimmed)
%  splitstr: regular expression used to split list (default='\s*[;\|]+\s*')
% OUTPUTS
%    listex: nx1 list with expanded elements (n>=m, n==m when all elements are singleton)
%         i: nx1 index vector (with values between 1 and m) such that listex is member of list(i)
%         j: mx1 cell array such that j(k)=(1:lk)' where lk is the number of elements in the kth sublist (k=1..m )
%
% NB: trailing separator(s) string are discarded
%
%
%   See also: FMECAENGINE KEY2KEY ISMEMBERLIST

% Migration 2.0 (Fmecaengine v0.50) - 16/07/2011 - INRA\Olivier Vitrac rev. 22/05/2014

% Revision history
% 17/07/2011 add | as separator, remove trailing separator and spaces
% 22/05/2014 force strtrim, returns correctly '' and []


% default
splitstr_default = '\s*[;\|]+\s*';

% arg check
if nargin<1, error('one argument is required'); end
if nargin<2, splitstr=[]; end
if isempty(splitstr), splitstr = splitstr_default; end
if isempty(list), listex = list; return, end
if ischar(list), list = {list}; end
if ~iscellstr(list), error('list must be a cell array or a string'), end
m = numel(list);

% split
c = regexp(strtrim(regexprep(list,['(^' splitstr ')|(' splitstr '$)'],'')),splitstr,'split');
nc = cellfun('length',c);
nnc = sum(nc);
listex = cell(nnc,1);
idx  = zeros(nnc,1);
jdx = cell(m,1);
pos = 1;
for k=1:m
    u = pos+(0:nc(k)-1);
    listex(u) = strtrim(c{k}); % added strtrim 22/05/2014
    idx(u) = k;
    jdx{k} = (1:nc(k))';
    pos = pos + nc(k);
end

% outputs
if nargout>1, iout = idx; end
if nargout>2, jout = jdx; end