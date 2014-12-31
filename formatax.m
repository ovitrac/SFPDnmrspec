function formatax(hs,varargin)
%FORMATAX format axes for publication (to match GRACE design: http://en.wikipedia.org/wiki/Grace_(plotting_tool))
%   syntax: formatax
%           formatax(hs,'property1',value1,'property2',value2,'property3',value3,..)
%      hs: axes handles (default = all handles in the current figure)
%      property/value any valid axes properties

% MS 2.1 - 27/12/09  - INRA\Olivier Vitrac - rev. 19/12/12

% Revision history
%19/02/11 rethrow error message
%19/12/111 fix TickDir

arg_default = struct(...
'FontSize'  ,12,...
'Box'       ,'on',...
'FontName'  ,'Arial',...
'FontUnits' ,'points',...
'FontWeight','normal', ...
'LineWidth' ,1.5,...
'TickLength',[.01 .0250]*2,...
'TickDirMode','auto',...
'TickDir','in' ...
);

[arg,remain] = argcheck(varargin,arg_default);

% arg check
if nargin<1, hs = []; end
if isempty(hs)
    hs = get(gcf,'children');
    hs = hs(strcmp(get(hs,'type'),'axes'));
    hs = hs(end:-1:1);
end
if isempty(hs), hs = gca; end

set(hs,arg)
if ~isempty(remain)
    try
        set(hs,remain{:})
    catch %#ok<CTCH>
        error('invalid additional (axes) arguments\n\t%s',lasterr) %#ok<LERR>
    end
end