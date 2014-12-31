function scfd(data,varargin)
%  SCFD set current figure data (as collected with GCFD)
%   syntax: scfd(data)
%   options: scfd(data [,'noaxes','nolegend'])

% INRA\Olivier Vitrac 28/10/02 - rev. 14/12/11

% REVISION HISTORY
% 29/10/02 release candidate
% 09/11/09 add text object
% 20/02/11 add nolegend
% 31/08/11 fix nolegend
% 01/09/11 prevent empty strings from being printed (text(.5,.5,[]) generates an error)
% 14/12/11 add pacth

% object properties
propAXESlist = {'sXlabel','sYlabel','nXscale','nYscale','nXlim','nYlim'};
propLINElist = {'Color','LineStyle','LineWidth','Marker','MarkerSize','MarkerEdgeColor','MarkerFaceColor','Xdata','Ydata','Zdata'};
propIMAGElist = {'Xdata','Ydata','CData','AlphaData','CDataMapping','AlphaDataMapping','tag'};
propTEXTlist  = {'BackgroundColor','Color','EdgeColor','FontAngle','FontName','FontSize','FontUnits','FontWeight','HorizontalAlignment','Interpreter',...
    'LineStyle', 'LineWidth','Margin','Position','Rotation','String','Units','VerticalAlignment'};
propPATCHlist = {'CData','FaceVertexAlphaData','FaceVertexCData','EdgeAlpha','EdgeColor','FaceAlpha','FaceColor','Faces','LineStyle','LineWidth',...
    'Marker','MarkerEdgeColor','MarkerFaceColor','MarkerSize','Vertices','XData','YData','ZData'};
typeOBJECTlist = {'line','image','text','patch'};


% arg check
if nargin<1, error('... required 1 input'), end
if nargin>1, options = varargin; else options = {''}; end
if isempty(data), disp('... no data'), end
optionsl = lower(options);
noaxes = ismember('noaxes',optionsl);
nolegend = ismember('nolegend',optionsl);
fdata = fieldnames(data)';
validfields = intersect(fdata,typeOBJECTlist);
if isempty(validfields), error('SCFD: unknown object (no valid fields)'), end

% new figure
if ~noaxes
    figure
    w = 4;
    m = ceil(length(data)/w);
    n = min(length(data),w);
end

% scan all axes    
for i=1:length(data)
    
    % Axes
    if ~noaxes
        hs = subplot(m,n,i); set(hs,'box','on')
    end
    hold on

    % Objects
    for eachfield = validfields
        typeOBJECT = eachfield{1};
        if ~isempty(getfield(data,{i},typeOBJECT))    
            % Objects properties
            switch lower(typeOBJECT)
                case 'line'
                    propOBJECTlist = propLINElist;
                    nline = length(data(i).line);
                    hp = zeros(nline,1);
                    legp = cell(nline,1);
                    for j = 1:nline
                        if any(data(i).line(j).Xdata)
                            hp(j) = plot(data(i).line(j).Xdata,data(i).line(j).Ydata,'ko');
                        else
                            hp(j) = plot(data(i).line(j).Ydata,'ko');
                        end
                        legp{j} = sprintf('line %d',j);
                        for p = propOBJECTlist
                            value = getfield(data,{i},'line',{j},p{1});
                            if any(value), set(hp(j),p{1},value), end
                        end
                    end
                    hl = [];
                    if  ~nolegend
                        if nline>1 && nline<5, hl = legend(hp,legp,0);
                        elseif nline>4 && nline<10, hl = legend(hp,legp,-1); end
                        if any(hl), set(hl,'box','off'), end
                    end
                case 'image'
                    propOBJECTlist = propIMAGElist;
                    for j = 1:length(data(i).image)
                        hp = image(data(i).image(j).Xdata,data(i).image(j).Ydata,data(i).image(j).CData);
                        for p = propOBJECTlist
                            value = getfield(data,{i},'image',{j},p{1});
                            if any(value), set(hp,p{1},value), end
                        end
                    end
                case 'text'
                    propOBJECTlist = propTEXTlist;
                    for j = 1:length(data(i).text)
                        if ~isempty(data(i).text(j).String)
                            hp = text(data(i).text(j).Position(1),data(i).text(j).Position(2),data(i).text(j).Position(3),data(i).text(j).String);
                            for p = propOBJECTlist
                                value = getfield(data,{i},'text',{j},p{1});
                                if ~isempty(value), set(hp,p{1},value), end
                            end
                        end
                    end
                case 'patch'
                    propOBJECTlist = propPATCHlist;
                    for j = 1:length(data(i).patch)
                        tmp = [propOBJECTlist;struct2cell(data(i).patch(j))'];
                        tmp = tmp(:,~cellfun('isempty',tmp(2,:)));
                        patch(tmp{:});
                    end
            end % end switch
        end % endif ~isempty
    end % endfor eachfield
    
    % Axes properties
    if ~noaxes
        for p = propAXESlist
            value = getfield(data,{i},p{1}(2:end));
            if ~isempty(value) && p{1}(1)=='n'
                set(hs,p{1}(2:end),value);
            end
        end
        xlabel(data(i).Xlabel)
        ylabel(data(i).Ylabel)
        title(sprintf('AXES [\\bf%d\\rm]',i))
    end
end