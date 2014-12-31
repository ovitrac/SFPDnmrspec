function [haxes_out,hlines_out,ht_out] = refslope(haxorigin,haxdestination,slope,direction,lineoptions,textoptions,box,vars,defaultposition)
%REFSLOPE create reference slopes in a new axes based on a coordinates system of an existing axes
%   syntax: refslope([haxorigin,haxdestination/position,slope,direction,lineoptions,textoptions,box,vars,defaultposition])
%           [haxes,hlines] = refslope(...)
%   Inputs:
%           haxorigin: handles of initial axes (mx1 vector) (default = gca)
%                    >> semi-log axes generate an error (linear and log-log plots are accepted)
%           haxdestination: handles of subplots, where reference slopes must be plotted (line segments)
%           position: position of destination axes to be created in haxorigin
%                     1: upper right corner, 2: upper left, 3: bottom left, 4: bottom right (default)
%           slope: nx1 vector of slopes(default value = [.5 1 2])
%           direction: nx1 vector of directions in destination axes (circular permutation are applied)
%                     1: from upper right corner, 2: from upper left, 3: from bottom left(default), 4: from bottom right
%           lineoptions: a property object created with plotpub (ex. lineoptions=plotpub('linewidth',1))
%                   >> circular permutation are applied if its length does not match n
%           textoptions: a structure accepted by text (ex. textoptions = struct('fontsize',12);
%                   >> circular permutation are applied if its length does not match n
%           box: an integer code betwee 0 and 3 (circular permutation are applied)
%               0: no legend, no axes
%               1: legend only
%               2: axes only
%               3: legend and axes (default)
%           vars: variable names, default = {'x' 'y'}
%           defaultposition: template position to create refslope axes (assume position=3), default = [0.2 0.2 0.2 0.2]
%    Outputs:
%           haxes: as haxdestination (usefull when axes are created)
%       >>  haxes are taggeg as 'refslope'
%           hlines: handles of the reference lines (mxn structure array)
%
%    Note:   prior referece slopes are not deleted (this behevior is different from LEGENDPUB)
%           use delete(hr) to remove the reference slopes with handle h
%
%    Example 1: plots 10 reference slopes in a previously created axes
%       clf, ho=gca; plot([0 10],[0 20]), hd = axes('position',[.6 .2 .2 .2]);
%       refslope(ho,hd,.5:.5:5,3,[],struct('fontsize',7))
%    Example 2: log-log plot, 2 subplots, the destination axes are created
%       clf, ho=subplots(1,[1 1],0,.1); subplot(ho(1)),plot([1 100],[1 10]), subplot(ho(2)),plot([1 10],[1 100]), set(ho,'xscale','log','yscale','log')
%       [hax,hlines] = refslope(ho,[2 4],.5:.5:2)
%    Example 3: semi-log plot, 2 subplots, the destination axes are created
%       x = linspace(0,5,1000);
%       clf, ho=subplots(1,[1 1],0,.1);
%       subplot(ho(1)),plot(x,exp(x)), subplot(ho(2)),plot(exp(1+x),x)
%       set(ho(1),'xscale','linear','yscale','log'), set(ho(2),'xscale','log','yscale','linear')
%       [hax,hlines] = refslope(ho,[2 4],.5:.5:2)
%
%   See also: PLOTPUB LEGENDPUB

% MS 2.0 - 29/08/07 - Olivier Vitrac - rev.  07/03/12

% Revision
% 30/08/07: add position (create axes when required), box, vars, defaultposition
% 23/09/07: increase the linewidth of axes by 0.2
% 19/05/09: add ht_out
% 26/03/11: add semilog
% 07/03/12 force unit scaling for axes

%
haxdestination_default = 4; % default position
slope_default = [.5 1 2];
direction_default = 3; % default direction
textoptions_default = struct('fontsize',10);
sym = struct('center',[1 1;0 1;0 0;1 0],'coeff',[-1 -1;1 -1;1 1;-1 1 ]);
rmargin = [1.05 1.05]; % 5% margin
box_default = 3;
vars_default = {'x' 'y'};
defaultposition_default = [.2 .2 .2 .2];

% arg check
if nargin>9, error('syntax: hr=refslope(haxorigin,[haxdestination/position,slope,lineoptions,textoptions,box,vars,defaultposition])'), end %#ok<NCHK>
if nargin<1, haxorigin = gca; end
if any(~ishandle(haxorigin)), error(' some inial axes are invalid'), end
nax = length(haxorigin);
if nargin<2, haxdestination = []; end
if nargin<3, slope = []; end
if nargin<4, direction = []; end
if nargin<5, lineoptions = []; end
if nargin<6, textoptions = []; end
if nargin<7, box=[]; end
if nargin<8, vars = []; end
if nargin<9, defaultposition = []; end
haxdestination = haxdestination(:);
haxdestination = [haxdestination;NaN(max(0,nax-length(haxdestination)),1)];
haxdestination(isnan(haxdestination)) = haxdestination_default;
if isempty(slope), slope = slope_default; end
slope = sort(slope(:));
nslope = length(slope);
if isempty(direction), direction = direction_default; end
for j=1:numel(direction)
    if ~ismember(direction(j),1:4), error('unknown direction value, must be an integer between 1 and 4'), end
end
direction = direction(mod(0:nax-1,numel(direction))+1);
if isempty(lineoptions)
    col = gray(nslope+1); col = flipud(col(1:end-1,:));
    lineoptions = plotpub('linestyle','-','linewidth',.4,'color',col,'xscale','linear','yscale','linear','marker','none','box','off');
end
lineoptions = plotpub(lineoptions,'marker','none');
lineoptionsaxes = plotpub(lineoptions,'linewidth',lineoptions.linewidth+.2,'color','k');
if isempty(textoptions), textoptions = textoptions_default; end
textoptions = textoptions(mod(0:nslope-1,numel(textoptions))+1);
if isempty(box), box = box_default; end
for j=1:numel(box)
    if ~ismember(box(j),0:3), error('unknown box value, must be an integer between 0 and 3'), end
end
box = box(mod(0:nax-1,numel(box))+1);
if isempty(vars), vars = vars_default; end
if ~iscell(vars) || ~length(vars)==2, error('vars must be a cell array such as {''label x'' ''label y}'''), end
if isempty(defaultposition), defaultposition = defaultposition_default; end
if ~isnumeric(defaultposition) || length(defaultposition)~=4, error('defaultposition = [xbottomleftcorner ybottomleftcorner width height]'), end

% for each axes
currentaxes = gca; % store current axes
hlines = [];
for iax = 1:nax
    % check all axes
    if ~strcmp(get(haxorigin(iax),'type'),'axes'), error('AXES(%d) is not a valid initial axes',iax), end
    if ~ishandle(haxdestination(iax)) || ~strcmp(get(haxdestination(iax),'type'),'axes') % create axes if required
        if ismember(haxdestination(iax),1:4)
            opos = outerpos(haxorigin(iax));
            dpos = [ opos(1)+ sym.center(haxdestination(iax),1)*opos(3)*(1-defaultposition(3)) + sym.coeff(haxdestination(iax),1)*opos(3)*defaultposition(1), ...
                     opos(2)+ sym.center(haxdestination(iax),2)*opos(4)*(1-defaultposition(4)) + sym.coeff(haxdestination(iax),2)*opos(4)*defaultposition(2), ...
                     opos(3)*defaultposition(3),...
                     opos(4)*defaultposition(4) ];
            haxdestination(iax) = axes('position',dpos);
        else
            error('AXES(%d): invalid destination axes',iax)
        end
    end
    % plot
    set(haxdestination(iax),'xlim',[0 1],'ylim',[0 1]) % force unit scaling for axes (added 7/3/12)
    xscale = get(haxorigin(iax),'xscale');
    yscale = get(haxorigin(iax),'yscale');
    % Old behavior
    %if ~strcmp(xscale,yscale), error('AXES(%d): semilog plots are not accepted',iax), end
    %islog = strcmp(xscale,'log');
    islogx = strcmp(xscale,'log');
    islogy = strcmp(yscale,'log');
    ratio = findsize(haxdestination(iax))./findsize(haxorigin(iax)); % <1
    newscale = findscale(haxorigin(iax)).*ratio;
    xplot = repmat({[0;0]},1,nslope); yplot = xplot;
    subplot(haxdestination(iax)), hold on
    if ismember(direction(iax),[1 4]), alignment = 'right'; else alignment = 'left'; end
    if ismember(direction(iax),[2 4]), rotation = -1; else rotation=1; end
    ht = cell(nslope,1);
    for j = 1:nslope
        theta = atan(slope(j)*newscale(1)/newscale(2));
        xplot{j}(1) = sym.center(direction(iax),1);
        yplot{j}(1) = sym.center(direction(iax),2);
        xplot{j}(2) = sym.center(direction(iax),1) + sym.coeff(direction(iax),1)*cos(theta);
        yplot{j}(2) = sym.center(direction(iax),2) + sym.coeff(direction(iax),2)*sin(theta);
        xtxt        = [ sym.center(direction(iax),1) + sym.coeff(direction(iax),1)*cos(theta)*rmargin(1),...
            sym.center(direction(iax),2) + sym.coeff(direction(iax),2)*sin(theta)*rmargin(2)];
        % Old behavior
        %if ~islog
        %    if slope(j)==1, txt = sprintf('%s=%s',vars{2},vars{1}); else txt = sprintf('%s=%0.3g\\cdot%s',vars{2},slope(j),vars{1}); end
        %else
        %    if slope(j)==1, txt = sprintf('%s\\propto%s',vars{2},vars{1}); else txt = sprintf('%s\\propto%s^{%0.3g}',vars{2},vars{1},slope(j)); end
        %end
        if ~islogx && ~islogy
            if slope(j)==1, txt = sprintf('%s=%s',vars{2},vars{1}); else txt = sprintf('%s=%0.3g\\cdot%s',vars{2},slope(j),vars{1}); end
        elseif ~islogx && islogy
            if slope(j)==1, txt = sprintf('log(%s)=%s',vars{2},vars{1}); else txt = sprintf('log(%s)=%0.3g\\cdot%s',vars{2},slope(j),vars{1}); end            
        elseif islogx && ~islogy
            if slope(j)==1, txt = sprintf('%s=log(%s)',vars{2},vars{1}); else txt = sprintf('%s=%0.3g\\cdotlog(%s)',vars{2},slope(j),vars{1}); end
        else
            if slope(j)==1, txt = sprintf('%s\\propto%s',vars{2},vars{1}); else txt = sprintf('%s\\propto%s^{%0.3g}',vars{2},vars{1},slope(j)); end
        end
        if ismember(box(iax),[1 3]) % legend
            ht{j} = text(xtxt(1),xtxt(2),txt,...
                    'horizontalalignment',alignment,'verticalalignment','middle','rotation',rotation*180*theta/pi,textoptions(j));
        end
    end
    hlines = [hlines;plotpub(xplot,yplot,lineoptions)]; %#ok<AGROW>
    if ismember(box(iax),3) % legend
        text(sym.center(direction(iax),1)+sym.coeff(direction(iax),1)*rmargin(1),sym.center(direction(iax),2),vars{1},'horizontalalignment',alignment,'verticalalignment','middle',textoptions(j))
        text(sym.center(direction(iax),1),sym.center(direction(iax),2)+sym.coeff(direction(iax),2)*rmargin(2),vars{2},'rotation',rotation*90,'horizontalalignment',alignment,'verticalalignment','middle',textoptions(j))
    end
    if ismember(box(iax),[2 3]) % axes
        plotpub(sym.center(direction(iax),1)+[0 sym.coeff(direction(iax),1)],[0 0]+sym.center(direction(iax),2),lineoptionsaxes)
        plotpub([0 0]+sym.center(direction(iax),1),sym.center(direction(iax),2)+[0 sym.coeff(direction(iax),2)],lineoptionsaxes)
    end
    subplot(haxdestination(iax)), axis tight, ax = axis;
    set(haxdestination(iax),'visible','off','xlim',sort(sign(ax(1:2))),'ylim',sort(sign(ax(3:4)))) % force unit scaling for axes (added 7/3/12)
%     hr(:,iax) = tmp
end

% outputs
subplot(currentaxes) % restore current axes
set(haxdestination,'tag','refslope')
if nargout, haxes_out = haxdestination; end
if nargout>1, hlines_out = hlines; end
if nargout>2, ht_out = ht; end 


%%%%%%%%%%%%%%%%%%%%%%
% private functions
%%%%%%%%%%%%%%%%%%%%%%
function pos=outerpos(h)
% find the outerposition in normalized units of any axes with handle h
newunits = 'normalized';
oldunits = get(h,'units');
isnew    = strcmp(oldunits,newunits);
if ~isnew, set(h,'units',newunits); end
pos = get(h,'outerposition');
if ~isnew, set(h,'units',oldunits); end

function siz=findsize(h)
% find the size in normalized units of any axes with handle h
newunits = 'normalized';
oldunits = get(h,'units');
isnew    = strcmp(oldunits,newunits);
if ~isnew, set(h,'units',newunits); end
pos = get(h,'position');
if ~isnew, set(h,'units',oldunits); end
siz = pos(3:4);

function scale=findscale(h)
% find the scale of any axes with handle h
xax = get(h,'xlim');
yax = get(h,'ylim');
if any(isinf(xax)), error('xlim contains +/- Inf values'), end
if any(isinf(yax)), error('ylim contains +/- Inf values'), end
%Old Fashion
% if strcmp(get(h,'xscale'),'log')
%     xax = log10(xax);
%     yax = log10(yax);
% end
if strcmp(get(h,'xscale'),'log'), xax = log(xax); end
if strcmp(get(h,'yscale'),'log'), yax = log(yax); end
scale = [diff(xax) diff(yax)];