function data = gcfd(figurehandle)
%  GCFD get current figure data
%   syntax: a=gcfd;
%           a=gcfd(gcf);
%
%   See also: scfd

% INRA\Olivier Vitrac 28/10/02 - rev. 03/05/14

% REVISION HISTORY
% 29/10/02 release candidate
% 09/11/09 add text
% 14/12/11 add figurehandle, add patch objects, help update
% 28/12/11 add Position, Units, fix texts on several lines
% 03/05/14 add handles to fields

% Default properties
propAXESlist  = {'sXlabel','sYlabel','nXscale','nYscale','nXlim','nYlim','nPosition','nUnits'};
propLINElist  = {'Color','LineStyle','LineWidth','Marker','MarkerSize','MarkerEdgeColor','MarkerFaceColor','Xdata','Ydata','Zdata'};
propIMAGElist = {'Xdata','Ydata','CData','AlphaData','CDataMapping','AlphaDataMapping','tag'};
propTEXTlist  = {'BackgroundColor','Color','EdgeColor','FontAngle','FontName','FontSize','FontUnits','FontWeight','HorizontalAlignment','Interpreter',...
    'LineStyle', 'LineWidth','Margin','Position','Rotation','String','Units','VerticalAlignment'};
propPATCHlist = {'CData','FaceVertexAlphaData','FaceVertexCData','EdgeAlpha','EdgeColor','FaceAlpha','FaceColor','Faces','LineStyle','LineWidth',...
    'Marker','MarkerEdgeColor','MarkerFaceColor','MarkerSize','Vertices','XData','YData','ZData'};
typeOBJECTlist = {'line','image','text','patch'};
excludedTAG = {'legend','Colorbar'};
data = [];
eol = char(10);

% inputs
if nargin<1, figurehandle = gcf; end
if ~ishandle(figurehandle), error('the input must be a handle'); end
if ~get(figurehandle,'type'), error('the object matching the supplied handle is not a figure'); end

% reposition and children
set(figurehandle,'units','normalized','position',[0 0 .5 .5])
hc = get(figurehandle,'children');
indisaxes = find(strcmp(get(hc,'type'),'axes'));
if any(indisaxes)
    n = length(indisaxes);
    dispf('...%d axes objects found in %d children',n,length(hc))
else
    disp('no axes objects'), return
end
i = 0;
% For each axes
for hi = flipud(hc(indisaxes))'
    tagname = get(hi,'tag');
    if ismember(tagname,excludedTAG)
        dispf('AXES [%d] = ''%s'' is excluded',i+1,tagname)
    else
        i = i+1; dispf('AXES [%d]',i)
        hp = flipud(get(hi,'children'));
        data = setfield(data,{i},'handles','axes','main',hi);
        % Axes properties
        for p = propAXESlist
            if any(get(hi,p{1}(2:end)))
                if p{1}(1)=='s'
                    haxprop = get(hi,p{1}(2:end));
                    data = setfield(data,{i},p{1}(2:end),get(haxprop,'string'));
                    data = setfield(data,{i},'handles','axes',p{1}(2:end),haxprop);
                elseif p{1}(1)=='n'
                    data = setfield(data,{i},p{1}(2:end),get(hi,p{1}(2:end)));
                end
            else
                data = setfield(data,{i},p{1}(2:end),[]);
            end
        end
        % Object properties
        for typeOBJECT = typeOBJECTlist
            switch lower(typeOBJECT{1})
            case 'line' , propOBJECTlist = propLINElist;
            case 'image', propOBJECTlist = propIMAGElist;
            case 'text' , propOBJECTlist = propTEXTlist;
            case 'patch', propOBJECTlist = propPATCHlist;
            end
            indisobject = find(strcmp(get(hp,'type'),typeOBJECT));
            m = length(indisobject);
            if any(indisobject)
                dispf('...%d ''%s'' objects found in %d children',m,typeOBJECT{1},length(hp));
                j = 0;
                for hl=hp(indisobject)'
                    j= j+1;
                    data = setfield(data,{i},'handles',typeOBJECT{1},{j},hl);
                    for p = propOBJECTlist
                        value = get(hl,p{1});
                         if ~isempty(value)
                            data = setfield(data,{i},typeOBJECT{1},{j},p{1},value);
                        else
                             data = setfield(data,{i},typeOBJECT{1},{j},p{1},[]);
                         end
                    end
                    % fix replace strings stored as array of strings by single strings including \n
                    if strcmp(typeOBJECT,'text') && size(data(i).text(j).String,1)>0
                        tmp = [cellstr(data(i).text(j).String)'; repmat({eol},1, size(data(i).text(j).String,1))];
                        data(i).text(j).String = [tmp{1:end-1}]; %#ok<AGROW>
                    end
                end
            end
        end
    end
end
