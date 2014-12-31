function hleg=legendpubtab(hax,headers,txtparam,varargin)
%LEGENDPUBTAB display legend as a table
%   legendpubtab(hax,headers,txtparam,hleg1,leg1,hleg2,leg2,...)
%   hl = legendpubtab(...)
%            hax : valid axes handle
%        headers : cell array
%       txtparam : cell of property/value (for text objects)
%             hl : structure containing handles of different plotted objects
%
% MS 2.1 - 29/04/10 - Olivier Vitrac - rev. 26/03/11

% Revision history
% 09/07/10 fix variable column length
% 26/03/11 fix help


% Definitions
prop = struct( 'line',{{'linestyle','linewidth','color'}},...
               'marker',{{'marker','markersize','markeredgecolor','markerfacecolor'}},...
               'text',{{'fontsize','fontname','fontunits','fontweight','fontangle','linewidth','linestyle','backgroundcolor','edgecolor','color','interpreter'}} ...
              ); % list of accepted properties for legends


% arg check
if nargin<4, error('three outputs are at least required, syntax: legendpubtab(hax,headers,txtparam,hleg1,leg1,hleg2,leg2,...)'), end
if isempty(hax), hax = gca; end
if length(hax)>1 || ~ishandle(hax(1)), error('invalid axes handle'), end
if isempty(headers), title_relativeheight = 0; else title_relativeheight = 2; end
if ischar(headers), headers = {headers}; end
if isempty(txtparam), txtparam = {}; end
if ~iscell(txtparam), error('textparam must be a cell property/value'), end


% text
cax = gca;
subplot(hax), hold on
n = length(varargin); % number of columns
m = max(cellfun('length',varargin)); %length(varargin{1}); % number of rows
xleg = (0:n-1)/n; w = xleg(2)-xleg(1);
iline = cumsum([title_relativeheight;ones(m-1,1)]);
ylegall = (1-iline/iline(end)); hall = ylegall(1)-ylegall(2);
ylegall = ylegall + hall/2;
hleg = struct('cell',{cell(m,n)},'headers',zeros(1,n),'box',[]);

for j=1:n
    iline = cumsum([title_relativeheight;ones(length(varargin{j})-1,1)]);
    yleg = (1-iline/iline(end)); h = yleg(1)-yleg(2);
    yleg = yleg + h/2;
    txtheader = headers{min(j,length(headers))};
    if isstruct(varargin{j}) && isfield(varargin{j},'leg')
        allobjects = setdiff(fieldnames(varargin{j}),'leg');
        for i=1:min(m,length(varargin{j}))
            for object = allobjects(:)' % for all valid objects
                if any(varargin{j}(i).(object{1})) && ishandle(varargin{j}(i).(object{1})(1))
                    propvalues = get(varargin{j}(i).(object{1})(1),prop.(object{1})); % get current values
                    param = reshape([prop.(object{1})(:) propvalues(:)],length(propvalues),2)'; % pairwise property value
                    switch object{1}
                        case 'line'
                            hleg.cell{i,j}(end+1) = plot([1/4 3/4]*w+xleg(j),[0 0]+yleg(i),param{:});
                        case 'marker'
                            hleg.cell{i,j}(end+1) = plot(w/2+xleg(j),yleg(i),param{:});
                        case 'text'
                            txt = get(varargin{j}(i).text(1),'string');
                            hleg.cell{i,j}(end+1) = text(xleg(j),yleg(i),txt,'horizontalalignment','left','verticalalignment','middle',param{:});
                        otherwise
                            error('invalid legend object %s',object{1})
                    end
                end
            end % each object
        end
        xheader = xleg(j);
        hleg.headers(j) = text(xheader+w/2,1-hall/2,txtheader,'horizontalalignment','center','verticalalignment','middle',txtparam{:});
        xheader = xheader+w;
    elseif iscellstr(varargin{j})
        mi = min(m,length(varargin{j}));
        for i=1:mi
            hleg.cell{i,j} = text(xleg(j),yleg(i),varargin{j}(i),'horizontalalignment','left','verticalalignment','middle',txtparam{:});
        end
        hleg.headers(j)= text(xheader,1-hall/2,txtheader,'horizontalalignment','center','verticalalignment','middle',txtparam{:});
        xheader = xheader + w;
        hleg.box(end+1) = line([w;w]+xleg(j),[0;1],'color','k','linewidth',1.5);
    else
        error('object i=%d of %d is unknown',j,n)
    end
    hleg.box(end+1) = line([0;0],[0;1],'color','k','linewidth',1.5);
    hleg.box(end+1:end+3) = line([0;1],[[1;1]-hall [0;0] [1;1]],'color','k','linewidth',1.5);
end

set(hax,'xlim',[-.01 1],'ylim',[0 1.01],'visible','off')
subplot(cax)
