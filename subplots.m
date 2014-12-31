function hout = subplots(x,y,xsep,ysep,varargin)
% SUBPLOTS generates subplots with fixed width and height [whithin the current axis]
%   syntax: h = subplots(x,y [,xsep,ysep,'property','value','strict'])
%    inputs
%       x:  vector of widths (see example below)
%       y:  vector of heights (idem)
%       xsep: scalar or vector of separators along x (if length does not match, a circular permutation is applied)
%       ysep: scalar or vector of separators along y (if length does not match, a circular permutation is applied)
%       [...] stands for extra value pairs properties (see plot or set for details)
%       Note 1: x and y are expressed in arbitrary units (normalized herafter) whereas xsep and ysep are expressed in
%               normalized units
%       Note 2: For very accurate positioning when multiple subplots are combined, a 'strict' mode normalizing both x
%              and xsep (idem for y and ysep) can be applied with the keyword 'strict':
%              h = subplots(...,'strict')
%       Extra pair properties (non standard in Matlab)
%       'alive':  remove all undesirable handles: delete(h(~alive))
%                 The 'alive' values are a logical ny x nx vector or matrix or a vector of indices
%       'position': handle or 1x4 vector
%               empty value, use the current figure with position = [0.1300    0.1100    0.7750    0.8150]
%               figure handle, use this figure with position = [0.1300    0.1100    0.7750    0.8150]
%               axes handle, replace these axes by the new combination of subplots
%               [xbottomleft ybottomleft width height]: position vector in normalized units
%   
%    outputs
%       h: matrix of handles of size length(y) x length(x)
%       >> to plot in the subplot corresponding to the ith row and jth column
%          uses subplot(h(i,j))
%
%   Example: 27 subplots fully legended (please update the example to your need)
%       % --------------
%       % Definitions
%       % --------------   
%       ni = 3; % to be updated
%       nj = 4; % to be updated
%       varleg = {'v_1' 'v_2' 'v_3' 'v_4'}; % to be updated
%       predleg = 'D (m^2\cdots^{-1})';     % to be updated
%       % --------------
%       % axes creation
%       % --------------
%       h = subplots([.5 1 1 1 1],[1 .5 1 .5 1 .5],.02,[.01 .03],'fontsize',8); % generates the subplots
%       delete(h(2:2:end,1)) % remove the undesired subplots
%       % --------
%       % subplots
%       % --------
%       for i = 1:ni % scan all groups rowwise
%           for j=1:nj % scan columnwise 
%               subplot(h(2*(i-1)+1,j+1)) % main
%               plot(1:4,1:4,'ro','markersize',6,'color','k') % to be replaced
%               axis tight
%               subplot(h(2*i,j+1)) % horiz distributions
%               histfitlog(1+randn(100,1)*3) % to be replaced
%               set(h(2*i,j+1),'ydir','reverse')
%           end
%           subplot(h(2*(i-1)+1,1)) % vert distributions
%           histfitlog(1+randn(100,1)*3,[],1,0,1) % to be replaced
%           set(h(2*(i-1)+1,1),'xdir','reverse','yaxislocation','right')
%       end
%       % -----------------------------
%       % axis range, ylim for each row
%       % -----------------------------
%       for i = 1:ni
%           yrange = get(h(2*(i-1)+1,1:end),'ylim');
%           yrange = reshape([yrange{:}],2,nj+1)';
%           set(h(2*(i-1)+1,:),'ylim',[min(yrange(:,1)) max(yrange(:,2))])
%       end
%       % ---------------------------------
%       % axis range, xlim for each column
%       % ---------------------------------
%       for j = 1:nj
%           xrange = get(h(:,j+1),'xlim');
%           xrange = reshape([xrange{:}],2,2*ni)';
%           set(h(:,j+1),'xlim',[min(xrange(:,1)) max(xrange(:,2))])
%       end
%       % ----------------------
%       % axes label corrections
%       % ----------------------
%       set(h(1:2:end,1),'yticklabel',' ')     % vert distributions
%       set(h(1:end-1,2:end),'xticklabel',' ') % xtick labels 
%       set(h(:,3:end),'yticklabel',' ')       % ytick labels 
%       % labels
%       hxlabel = get(h(end,2:end),'xlabel');  % xlabels 
%       for j = 1:nj, set(hxlabel{j},'string',varleg{j}), end
%       hylabel = get(h(1:2:end,1),'ylabel');  % ylabels 
%       set(h(1:2:end,1),'yaxislocation','left')
%       set([hylabel{:}],'string',predleg)
%       htitle = get(h(1:2:end,1),'title');    % titles
%       for i = 1:ni, set(htitle{i},'string',[char(96+i) ')']), end
%       % ----------------------
%       %   paper format
%       % ----------------------
%       set(gcf,'paperorientation','landscape','paperposition',[0.5079    1.5212   28.6615   17.9416])
          

% MS-MATLAB 1.0 - 23/06/04 - INRA\Olivier Vitrac - rev. 28/01/11

% revision history
% 11/09/07 add keywords: 'alive' and 'position'
% 26/05/10 set x/y ticklabelmode to 'auto' for single axes (after alive)
% 28/01/11 major update (rationale): argcheck is used instead of internal controls, add strict 


% definitions
xsep_default = 0.01;
ysep_default = 0.01;
default = struct('alive',[],'position',[]);
keywordlist = 'strict';

% arg check
if nargin<2, error('subplots requires at least 2 arguments'), end
if nargin<3, xsep = []; end
if nargin<4, ysep = []; end
if isempty(xsep), xsep = xsep_default; end
if isempty(ysep), ysep = ysep_default; end
x = x(:); y = y (:);
nx = length(x);
ny = length(y);
if nx==0, error('SUBPLOTS: empty x'); end
if ny==0, error('SUBPLOTS: empty y'); end
[options,otheroptions] = argcheck(varargin,default,keywordlist);
alive = true(ny,nx);
nopos = true;

% alive
if ~isempty(options.alive)
    alivetmp = options.alive;
    if isnumeric(alivetmp) && all(alivetmp>0) && all(alivetmp<=nx*ny) && all(alivetmp==round(alivetmp))
        alive(:)=false;
        alive(alivetmp)=true;
    else
        if numel(alive)~=nx*ny, error('the size of parameter ''alive'' does not match the size of subplots %dx%d',ny,nx), end
        alive = (alivetmp>0);
    end
end
% position
if ~isempty(options.position)
    postmp = options.position;
    if length(postmp)==4
        pos = postmp;
    elseif ishandle(postmp)
        typ = lower(get(postmp,'type'));
        if strcmp(typ,'axes')
            pos = get(postmp,'position'); %delete(postmp)
        elseif strcmp(typ,'figure')
            pos = [0.1300    0.1100    0.7750    0.8150];
        else error('invalid type ''%s'' as position/handle',typ);
        end
    elseif isempty(postmp)
        pos = [0.1300    0.1100    0.7750    0.8150];
    else
        error('invalid position values or handle')
    end
    nopos = false;
end

if nopos
    % current
    big = gca;
    pos = get(big,'position'); %get(big,'position');
    delete(big)
end

% autosizing
nxsep = length(xsep);
nysep = length(ysep);
if nx>1, xsep = xsep(mod(0:nx-2,nxsep)+1); else xsep = 0; end
if ny>1, ysep = ysep(mod(0:ny-2,nysep)+1); else ysep = 0; end
if ~options.strict % default 
    x = x/sum(x); y = y/sum(y);
    x = (x - sum(xsep)*x)*pos(3);
    y = (y - sum(ysep)*y)*pos(4);
else
    xtot = sum(x)+sum(xsep);
    ytot = sum(y)+sum(ysep);
    x = pos(3)*x/xtot;
    y = pos(4)*y/ytot;
    xsep = pos(3)*xsep/xtot;
    ysep = pos(4)*ysep/ytot;
    
end

    

% subplots
h = zeros(ny,nx);
yc = pos(2)+pos(4)-y(1);
for i=1:ny
    %ig = ny-i+1;
    xc = pos(1);
    for j=1:nx
        h(i,j) = axes('position',[xc yc x(j) y(i)],'box','on','xlimmode','auto','ylimmode','auto','xgrid','off','ygrid','off');
        if i~=ny, set(h(i,j),'xticklabel',' '), end
        if j>1, set(h(i,j),'yticklabel',' '), end
        if j<nx, xc = xc + x(j) + xsep(j); end
        if ~isempty(otheroptions), set(h(i,j),otheroptions{:}), end
    end
    if i<ny, yc = yc - ysep(i) - y(i+1); end
end

% alive
if ~all(alive(:))
    delete(h(~alive))
    h = h(alive);
    if length(h)<2, set(h,'xticklabelmode','auto','yticklabelmode','auto'), end
end
% output
if nargout>0, hout = h; end