function tab2csv(tab,filename,varargin)
%TAB2CSV converts a 2D table (simple structure, 2D cell array, 2D numerical array) to a CSV file
%   syntax: tab2csv(tab,filename [, property1,value1,property2,value2,...])
%       tab: any 2D table (simple structure, 2D cell array, 2D numerical array)
%  filename: name of file (nb: the extension .csv is not appended)
%
%  Recognized property/value pairs (works best with openoffice)
%       'separator', ',' (default)
%   'textdelimiter', '"' (default)
%
%   See also: loadods, xlstblread
%   Additional information on CSV: http://creativyst.com/Doc/Articles/CSV/CSV01.htm

% MS 2.1 - 13/04/2011 INRA\Olivier Vitrac - rev. 20/07/11

% Revision History
% 20/06/11 recast structure array into a structure will cell fields
% 21/06/11 fix cell2struct used to recast structure arrays
% 20/07/11 fix rescan (insert abs)

% default
options_default = struct(...
               'separator', ',',...
               'textdelimiter', '"' ...
               );
if isunix, eol = '\n'; else eol = '\r\n'; end
           
% arg check
if nargin<2, error('2 arguments are required'); end
options = argcheck(varargin,options_default);

% Reprocess structures
if isstruct(tab)
    headers = fieldnames(tab)';
    if numel(tab)>1
        tab = tab(:); tab = cell2struct(cellfun(@(f) {tab.(f)},headers,'UniformOutput',false),headers,2);
    end
    nheaders = length(headers);
    tab = struct2cell(tab);
    tablength = cellfun(@(x) numel(x),tab);
    tmp = cell(max(tablength),nheaders);
    for i=1:nheaders
        if iscell(tab{i})
            tmp(1:tablength(i),i) = tab{i};
        elseif isnumeric(tab{i})
            tmp(1:tablength(i),i) = num2cell(tab{i}(:),2);
        else
            error('field ''%s'' is neither a cell nor a numeric array',headers{i})
        end
    end
    tab = [headers;tmp];
end

% Indentification of the type of each cell
typ = cellfun(@(x) isnumeric(x)*1 + ischar(x)*2 , tab ); % 1 if numeric, 2 if char
if any(~typ(:))
    error('only char and numeric types can be mixed.\n(e.g. cell or structures containing other cell or structures are not supported)')
end
nrow = size(typ,1);
rescan = [1; sum(abs(diff(typ,1)),2)~=0];

% Protect double quotes
itxt = (typ==2);
tab(itxt) = strrep(tab(itxt),'"','""');

% write row-wise
robot = {['%g' options.separator] [options.textdelimiter(1) '%s' options.textdelimiter(1) options.separator]};
fid = fopen(filename,'w');
for i = 1:nrow
    if rescan(i) % the format of the ith row is different from the format of the ith-1
        format = [robot{typ(i,:)}];
        format = [format(1:end-1) eol]; % remove trailing separator
    end
    fprintf(fid,format,tab{i,:});
end
if fclose(fid), error('unable to write ''%s''',filename); end