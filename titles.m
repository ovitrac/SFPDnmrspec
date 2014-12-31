function htout = titles(hs,txt,varargin)
%TITLES display titles superposed to figures instead of above (as with TITLE)
%   syntax: ht=titles(hs,txt,'property',value,...)
%    hs: array of axes handles (default = axes in the current figure)
%   txt: cell array of text (circular permutations applied if size of hs and txt are not compatible)
%        default value = {'a', 'b' ... 'z'}
%   keyword = 'transpose' forces hs=hs'
%   property/value
%      x: array of x normalized positions (default=0.1)
%      y: array of y normalized positions (default=0.9)
% prefix: prefix used for numbering (default='')
% suffix: suffix used for numbering (default='.')
%   any valid pair text property/value ('fontsize',12)
%    ht = mx1 structure array with fields:
%       ht(i).axes axes handles superposed
%       ht(i).text text handles

% MS 2.1 09/12/09 - INRA\Olivier Vitrac - rev. 26/12/09

% revision history
% 10/12/09 add searchkeywords
% 26/12/09 use argcheck (remove private function)

% default values
default = struct(...
    'x', 0.1,...
    'y', 0.9,...
    'prefix', '',...
    'suffix', '.' ...
    ); 

% arg check
if nargin<1, hs = []; end
if isempty(hs)
    hs = get(gcf,'children');
    hs = hs(strcmp(get(hs,'type'),'axes'));
    hs = hs(end:-1:1);
end
if ~all(ishandle(hs)), error('some handles are invalid'), end
if ~all(ismember(get(hs,'type'),'axes')), error('some handles do not refer to axes'), end
nax = numel(hs);
if nargin<2, txt=''; end
[param,varargin] = argcheck(varargin,default,'transpose');

% keywords check
if param.transpose, hs = hs'; end
param.x = param.x(mod((1:nax)-1,numel(param.x))+1);
param.y = param.y(mod((1:nax)-1,numel(param.y))+1);
if isempty(txt), txt = arrayfun(@(i)sprintf('%s%s%s',param.prefix,char(96+i),param.suffix),1:nax,'UniformOutput',false); end
if ischar(txt),  txt = cellstr(txt); end
txt = txt(mod((1:nax)-1,numel(txt))+1);

% plots
ht = repmat(struct('axes',[],'text',[]),nax,1);
for ia=1:nax
    ht(ia).axes = axes('position',get(hs(ia),'position'),'visible','off');
    ht(ia).text = text(param.x(ia),param.y(ia),txt{ia},varargin{:});
end

% output
if nargout, htout = ht; end

end