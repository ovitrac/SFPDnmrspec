function hnewax = addax(hax,funfcn,pos,varargin)
% ADDAX add a new ax (x or y) to an existing axe
%   syntax: hnewax = addax(hax,funfcn,pos,p1,p2...)
%       hax: vector of valid axe handles
%       funfcn: function handle, which defined the new axe scale
%       pos: 'x' or 'y' (default = 'x')
%       p1,p2...: aditional arguments for funfcn
%       hnewax: vector of new axe handles

% Example (see pub2, figure 10 for details):
%     hs = subplots([1 1],[1 1 1],.05,.05);
%     set(hs(:,2),'xlim',[3 2e3],'xticklabelmode','auto','xscale','log')
%     hns2 = addax(hs(:,2),@activation_energy,'X',300,1e-12,'kT');
%     set(hns2(2:end),'xticklabel',' ')
%     set(hns2,'xminortick','on')
%     set(get(hns2(1),'xlabel'),'string','free energy/k_B\cdotT')
%     set(hs(1:end-1,2),'xticklabel',' ')
%     set(gcf,'paperposition',[-2.3281    1.3552   25.6403   26.9670])


% MS-MATLAB 1.0 - 23/07/05 - Olivier Vitrac - rev.

% definitions
pos_default = 'x';

% arg check
if ~nargchk(3,2,nargin), error('2 arguments are at least required (max=3)'), end
if nargin<3, pos = []; end
if isempty(pos), pos = pos_default; end
pos = upper(pos);
if ~isa(funfcn,'function_handle'), error('invalid function handle'), end
if ~ischar(pos), error('invalid axe name: pos = ''X'' or ''Y'''), end
switch pos
    case 'X'
        npos = 'Y';
    case 'Y'
        npos = 'X';
    otherwise
        error('invalid axe name: pos = ''X'' or ''Y''')
end

% calculations
nax = length(hax);
hnewax = zeros(size(hax));
set(hax,'box','off')
for i=1:nax
    axpos = get(hax(i),'position');
    axlim = get(hax(i),[pos 'lim']);
    newaxlim = funfcn(axlim,varargin{:});
    hnewax(i) = axes(   'position',axpos,...
                        'XAxisLocation','top',...
                        'YAxisLocation','right',...
                        'Color','none',...
                        'box','off',...
                        'tickdir','out');
    set(hnewax(i),  [npos 'ticklabel'],' ',...
                    [npos 'tick'],[],...
                    [pos 'lim'],newaxlim,...
                    [pos 'tickmode'],'auto',...
                    [pos 'ticklabelmode'],'auto');
end

