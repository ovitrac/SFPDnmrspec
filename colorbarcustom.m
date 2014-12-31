function [haxout,hpatch]=colorbarcustom(varargin)
%COLORBARCUSTOM add a customized colorbar
%   syntax: hax=colorbarcustom('property',value,keyword...)
%           [hax,hpatches] = colorbarcustom(...)
%      hax: handle of colorbar axes
%   hpatch: handles of patches 
%
%   Specific properties/values
%       'position': [x y W H] or existing axes handle
%       'colormap': nx3 RGB colormap (see colormap for details)
%          'ctick': colormap tick values
%    'orientation': 'horiz' (default) or 'vert'
%      'formatsci': structure to change the output of axes values as text (see formatsci to see the details of the fields)
%        'economy': to be combined with formatsci (default='economy')
%
%   Keywords: 'endbottom','endtop','endleft','endright','both' 
%              They add sharp ends instead of square ends.
%
%   Text properties (with sharp ends), list of implemented properties
%   'horizontalalignment', 'verticalalignment'
%   'fontangle','fontname','fontweight','rotation'
%
%
% %EXAMPLE 1: horizontal colorbar, log scale
%   figure
%   colorbarcustom('colormap',winter(12),'ctick',logspace(0,4,10),'xscale','log','fontsize',12,'position',[0.2 0.1 0.4 .03],'backgroundcolor','w')
% 
% %EXAMPLE 2: vertical colorbar, linear scale, lines between colors
%   [hc,hf]=colorbarcustom('orientation','vert','colormap',cbrewer('div','RdYlGn',15),'ctick',linspace(0,100,5),'fontsize',8,'position',[0.8 0.4 0.01 .4],'YAxisLocation','right');
%   formatax(hc,'fontsize',8), set(hf,'linestyle','-','edgecolor',[0.5 0.5 0.5])
%   ylabel('color scale','fontsize',8)
%
% %EXAMPLE 1a (as 1 with nice ends)
%   figure, colorbarcustom('colormap',winter(12),'ctick',logspace(0,4,10),'xscale','log','fontsize',12,'position',[0.2 0.1 0.4 .03],'both','xtick',[1 10 100 1000 1e4],'economy','none','rotation',45,'fontangle','italic')
%   
%
% %EXAMPLE 2a (as 2 with nice ends)
%   figure, colorbarcustom('orientation','vert','colormap',cbrewer('div','RdYlGn',10),'ctick',linspace(0,100,5),'fontsize',8,'position',[0.8 0.4 0.01 .4],'both','fontweight','bold','yaxislocation','right','rotation',45);

%MS 2.1 - 11/12/2012 - INRA\Olivier Vitrac - rev. 19/12/2012

% Revision history
% 13/12/2012 surf replaced by pacth to keep vectorial output when printing
% 12/12/2012 additional text properties
% 19/12/12 change default TickDir to 'out'
% 19/12/12 add backgroundcolor

% default
formatsci_default = struct(...
    'power10',@(x) floor(log10(x)),...
    'texpattern','%0.3g\cdot10^{%d}',...
    'pattern','%0.4g' ...
    );
defaultax = struct('fontsize',10,'tickdir','out');
defaultxt = struct('horizontalalignment','','verticalalignment','','fontangle','normal','fontname','Helvetica','fontweight','normal','rotation',0);
backgroundcolor = 'none';
keywords = {'endleft' 'endright' 'endup' 'enddown' 'both'};
default = struct('position',[],'colormap',jet(64),'ctick',[],'orientation','horiz','formatsci',formatsci_default,'economy','economy','backgroundcolor',backgroundcolor);

% arg check
[o,p]  = argcheck(varargin,default,keywords);
[ot,p] = argcheck(p,defaultxt);
p      = argcheck(p,defaultax,'','keep');
if isempty(o.position), o.position = gca; end
if length(o.position)==1 && ishandle(o.position) && strcmp(get(o.position,'Type'),'axes')
    hax = o.position;
    if ~isempty(p), set(hax,p{:}); end
else
    hax = axes('position',o.position,p);
end
if ~isempty(o.backgroundcolor) && ~strcmpi(o.backgroundcolor,'none')
    if strcmpi(o.backgroundcolor,'w'), o.backgroundcolor = [1 1 1]-sqrt(eps); end
    hax2 = axes('position',get(hax,'outerPosition'),'xticklabel',' ','ytickLabel',' ',...
        'color',o.backgroundcolor,'xColor',o.backgroundcolor,'yColor',o.backgroundcolor);
    uistack(hax,'top'), subplot(hax)
else
    hax2 = [];
end
n = size(o.colormap,1);
if ~ismatrix(o.colormap)~=2 && ~size(o.colormap,2)==3
    error('colormap should be a valid nx3 colormap')
end
ishoriz = strcmpi(o.orientation,'horiz');
if ishoriz,     tick = get(hax,'xtick'); scale = get(hax,'xscale'); ortholim = get(hax,'ylim'); orthotick = 'ytick';
else            tick = get(hax,'ytick'); scale = get(hax,'yscale'); ortholim = get(hax,'xlim'); orthotick = 'xtick';
end
if isempty(o.ctick), o.ctick = tick; end
p.(orthotick) = [-1 2];
if strcmp(scale,'linear'), val = linspace(o.ctick(1),o.ctick(end),n+1);
else                       val = logspace(log10(o.ctick(1)),log10(o.ctick(end)),n+1);
end

% main
[X,Y] = meshgrid(val,ortholim); Z = zeros(size(X));
if ~ishoriz; Z = X; X=Y; Y=Z; end
patches = surf2patch(X,Y,Z);
hold on
hs = zeros(n,1);
for i=1:n
    hs(i) = patch('Vertices',patches.vertices,'Faces',patches.faces(i,:),'FaceColor',o.colormap(i,:),'EdgeColor', 'none');
end

% change ends if required
endmodified = false;
if ishoriz
    if o.endleft || o.both
        delete(hs(1))
        patches.vertices(1,2) = mean(patches.vertices([1 4],2));
        hs(1) = patch('Vertices',patches.vertices,'Faces',patches.faces(1,[1 3 2]),'FaceColor',o.colormap(1,:),'EdgeColor', 'none');
        endmodified = true;
    end
    if o.endright || o.both
        delete(hs(end))
        patches.vertices(end,2) = mean(patches.vertices([end end-1],2));
        hs(end) = patch('Vertices',patches.vertices,'Faces',patches.faces(end,[1 3 4]),'FaceColor',o.colormap(end,:),'EdgeColor', 'none');
        endmodified = true;
    end
    if endmodified
        if strcmp(get(hax,'XAxisLocation'),'bottom')
            yref = 0;
            ot = argcheck(ot,struct('horizontalalignment','center','verticalalignment','top','fontsize',p.fontsize),'','keep');
        else
            yref = 1;
            ot = argcheck(ot,struct('horizontalalignment','left','verticalalignment','bottom','fontsize',p.fontsize),'','keep');
        end
        arrayfun(@(x) text(x,yref,formatsci(x,o.formatsci,o.economy),ot),get(hax,'xtick'));
        set(hax,'visible','off')
    end
else
    if o.enddown || o.both
        delete(hs(1))
        patches.vertices(1,1) = mean(patches.vertices([1 2],1));
        hs(1) = patch('Vertices',patches.vertices,'Faces',patches.faces(1,[1 3 2]),'FaceColor',o.colormap(1,:),'EdgeColor', 'none');
        endmodified = true;
    end
    if o.endup || o.both
        delete(hs(end))
        patches.vertices(end,1) = mean(patches.vertices([end end-1],1));
        hs(end) = patch('Vertices',patches.vertices,'Faces',patches.faces(end,[1 3 4]),'FaceColor',o.colormap(end,:),'EdgeColor', 'none');
        endmodified = true;
    end
    if endmodified
        if strcmp(get(hax,'YAxisLocation'),'left')
            xref = 0;
            ot = argcheck(ot,struct('horizontalalignment','right','verticalalignment','middle','fontsize',p.fontsize),'','keep');
        else
            xref = 1;
            ot = argcheck(ot,struct('horizontalalignment','left','verticalalignment','middle','fontsize',p.fontsize),'','keep');
        end
        arrayfun(@(y) text(xref,y,[formatsci(y,o.formatsci,o.economy) ' '],ot),get(hax,'ytick'));
        set(hax,'visible','off')
    end
end

% update axes
if ishoriz, set(hax,'yticklabel',' '); else set(hax,'xticklabel',' '); end
if ~endmodified, formatax(hax,p); end

% outputs
if nargout>0, haxout = [hax;hax2]; end
if nargout>1, hpatch = hs; end