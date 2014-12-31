function s=wraptext(s0,width,newline,tol,trim,forcecell,pattern,pattern_alternative)
%WRAPTEXT return a wrapped text with arbitrary newline sequence and accepting tolerance
%   syntax: s=wraptext(s0 [,width,newline,tol,trim])
%       width: required width (default=32)
%     newline: newline sequence (default='<br />')
%   forcecell: flag to force cell output ()
%         tol: +/-tolerance (default=5)
%        trim: flag (default=true), remove leading and trailing spaces

% Migration 2.0 - 08/05/2011 - INRA\Olivier Vitrac - rev. 08/10/13

% Revision history
% 28/01/12 add pattern, pattern_alternative (for internal purpose)
%           wraptext('specificvolumenative',10,'*',2,true,'[^aeiouAEIOU]{2,}') yields specific*volumenati*ve
% 08/10/2013 add trim in recursion (it was omitted), add forcecell

%default
pattern_default = '[\s-_\<\>\(\)\{\}\[\]\:\,\;\.]';
pattern_alternative_default = '[aeiouAEIOU]';

% arg check
if nargin<2, width = 32; end
if nargin<3, newline = '<br />'; end
if nargin<4, tol =5; end
if nargin<5, trim =true; end
if nargin<6, forcecell = false; end
if nargin<7, pattern=''; end
if nargin<8, pattern_alternative = ''; end
if isempty(pattern), pattern = pattern_default; end
if isempty(pattern_alternative), pattern_alternative = pattern_alternative_default; end

% wrap
if iscell(s0), s = cellfun(@(si) wraptext(si,width,newline,tol,trim,forcecell,pattern,pattern_alternative), s0,'UniformOutput',false); return, end
n = length(s0);
if n<=width, s=s0; return, end
if trim, s=strtrim(s0); else s=s0; end
p = regexp(s,pattern);
if isempty(p) && (length(s0)>width+tol) && ~isempty(pattern_alternative)
    p = regexp(s,pattern_alternative);
    shift = 0;
else
    shift = +1;
end
[dmin,imin] = min(abs(p-width));
if dmin<=tol, p=p(imin); else p=width; end
s = [s(1:p+shift-1) newline wraptext(s(p+shift:end),width,newline,tol,trim,false,pattern,pattern_alternative) ];
if forcecell, s = regexp(s,newline,'split'); end
