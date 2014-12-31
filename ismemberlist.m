function [tf,locout] = ismemberlist(A,B,splitstr,keepBorder,regularexpr)
%ISMEMBERLIST generalizes ismember to semi-colon separated lists (work only on cell arrays of strings or strings)
%   Syntax: [tf,loc] = ismemberlist(A,B [,splitstr,keepBorderforsingleton])
%          A: 1xnA or nAx1 cell array of strings
%          B: 1xnB or nBx1 cell array of strings
%   splitstr: regular expression used to split list (default='\s*[;\|]+\s*')
% keepBorder: flag to force to keep order when B contains a non-singleton list (default = false)
%regularexpr: flag to force to accept regular expression in B{1} (default = false)
%
%   See also: FMECAENGINE KEY2KEY EXPANDTEXTASLIST

% Migration 2.0 (Fmecaengine v0.42) - 17/07/2011 - INRA\Olivier Vitrac rev. 18/07/11

%Revision history
%18/07/11 replace intersect by a second ismember, add regularexpr

% default
isregularexp = '^\s*\\(.*)\\\s*$';
splitstr_default = '\s*[;\|]+\s*';
keepBorder_default = false;
regularexpr_default = false;

% arg check
if nargin<2, error('2 arguments are required'), end
if nargin<3, splitstr = []; end
if nargin<4, keepBorder = []; end
if nargin<5, regularexpr = []; end
if ischar(A), A = {A}; end
if ischar(B), B = {B}; end
if ~iscellstr(A), error('A must be a string or cell array of strings'); end
if ~iscellstr(B), error('B must be a string or cell array of strings'); end
if isempty(splitstr), splitstr = splitstr_default; end
if isempty(keepBorder), keepBorder = keepBorder_default; end
if isempty(regularexpr), regularexpr = regularexpr_default; end

% expand A
[Ae,iA] = expandtextaslist(A,splitstr);

% expand B if it is not a regular expression (i.e. delimited between '\')
if regularexpr && (numel(B)==1) && ~isempty(regexp(B{1},isregularexp,'once'))
    tmp = regexp(B{1},isregularexp,'tokens');
    Be = tmp{1};
else
    Be = expandtextaslist(B,splitstr);
    regularexpr = false;
end

% join columns
boolind = true; % default behavior (as ismember does)
if regularexpr % join with regexp
    % STEP1: protection of {{ }}
    tmp = regexprep(Be{1},{'\{{2}' '\}{2}'},{'$($' '$)$'});
    tmp = regexprep(tmp,{'#' '{' '}'},{'*' '(' ')'});
    j = ~cellfun('isempty',regexp(Ae,regexprep(tmp,{'\$\(\$' '\$\)\$'},{'{' '}'})))'; % ismember is replaced by regexp
    loc = j;
else % join with ismember (include repeated values), intersect needs to be avoided
    [j,loc] = ismember(Ae,Be); % j is indexed as Ae
    % keep Be order if needed
    if keepBorder && (numel(B)==1) && (length(Be)>1) % B contains a non-singleton list (order of Be needs to be kept)
        j = find(j);
        [~,loc2] = ismember(Ae(j),Be); % the second ismember is used to retrieve the position in Be (with redundancy)
        [~,k] = sort(loc2);
        loc = loc(j(k));
        j = j(k);
        boolind = false;
    end
end

% indices in A (default=boolean, integers are used when the order needs to be kept)
if boolind
    tf = false(size(A));
    tf(iA(j)) = true;
else
    tf = iA(j);                         % full indices to keep the order
end

if nargout>1, locout = loc; end



