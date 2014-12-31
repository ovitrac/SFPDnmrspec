function dataout = loadods(infile,varargin)
% LOADODS: loads data from an open document spreadsheet (ods) file into a cell array or structure
% SYNTAX: data = loadods(filename [,options])
%         data = loadods(filename [,property1,value1,property2,value2,...])
%
% INPUTS:
%   filename:   String representing the path and name of the ods file
%   options:    structure with fields options.property = value
%     property/default value:
%         'sheetname',''    ... sheet to load (empty for the first sheet, 'all' for all sheets)
%             'blank',NaN   ... value for blank cells
%         'transpose',false ... true to transpose table/headers
%            'concat',true  ... true to concat column content
%           'headers',1     ... number of header lines
%           'struct',true   ... to transform the table as a structure with headers
%      'structarray',false  ... to transform the original structure of arrays into a structure array
%     'forceboolean',false  ... convert to a boolean any column with uniquely 0 and 1
%
% OUTPUTS:
%   data:       MxN cell array containing data from the spreadsheet
%               structure as data.header = column data
%               NB: not that header contains only a-z and A-Z characters
%
% TIP: NaN, Inf, -Inf, true, false are read litteraly
%
% SEE ALSO: XLSTBLREAD, XLSREAD

% MS 2.1 - 11/02/11 - INRA\Olivier Vitrac - rev. 17/01/14
% Initial code from (C) 2007 Alex Marten - alex.marten@gmail.com

% Revision history
%   09/02/11 major revision, release candidate
%   10/02/11 handle several sheets, NaN, Inf
%   11/02/11 fix the method to guess the name of the first worksheet
%   11/02/11 fix concat when all data are numeric
%   14/02/11 fix loading single sheet without headers
%   01/03/11 fix headers on multiple lines including empty cells
%   07/03/11 fix sheetname and row/column names starting with a digit
%   06/04/11 add structarray
%   13/04/11 fix headers on empty sheets
%   08/05/11 add MaxlengthField
%   16/09/11 force allsheets = true when the name of the first sheet cannot be determined
%   03/01/12 add forceboolean
%   24/01/12 fix header with item(3), based on a correct identification of office:body tag (starting from item 3)
%   17/01/14 fix test bounds when empty cells are mixed with NaN, etc.

% Set default options
validchars = '[^a-zA-Z0-9]'; % accepted characters for fields
MaxlengthField = 40;
options_default = struct(...
    'sheetname',''    ,... sheet to load (empty for the first sheet)
        'blank',NaN   ,... value for blank cells
    'transpose',false ,... true to transpose table/headers
       'concat',true  ,... true to concat column content
      'headers',1     ,... number of header lines
      'struct',true   ,...  to transform the table as a structure with headers
 'structarray',false,  .... to transform the original structure of arrays into a structure array
 'forceboolean',false ... convert to a boolean any column with uniquely 0 and 1
        );
options = argcheck(varargin,options_default);
if options.structarray, options.struct = true; end

% Check for the existence of the ods file
if ~ischar(infile), error('the first argument must be a string'); end
[infilestr,local]=lastdir(infile); if isempty(local), local = pwd; end
if isempty(regexpi(infilestr,'.ods$')), infile = [infile '.ods']; end
if ~exist(infile,'file'), error('the file ''%s'' does not exist in ''%s''',infilestr,local); end

% Create a temporary directory to unzip the ods file into
dispf('LOADODS:\tloading...')
fileinfo(infile);
dir_temp = fullfile(tempdir,infilestr);
if ~exist(dir_temp,'dir')
    if ~mkdir(tempdir,infilestr), error('Unable to create ''%s'' in ''',infilestr,tempdir); end
else
    dispf('\t reuse the existing directory ''%s''',dir_temp) 
end

% Unzip the contents of the ods file
unzip(infile,dir_temp);

% Load the XML file containing the spreadsheet data
try
    XMLfile = xmlread(fullfile(dir_temp,'content.xml'));
catch %#ok<CTCH>
    error('Unable to read the spreadsheet data');
end

% Parse down to the <office:spreadsheet> node
nodes = XMLfile.getChildNodes;
node = nodes.item(0);
nodes = node.getChildNodes;
item = 4; foundbody = false;
while (item>0) && ~foundbody
    item = item -1;
    if ~isempty(nodes.item(item))
        nodetype = nodes.item(item).toString;
        if strcmpi(nodetype,'[office:body: null]'), foundbody = true; end
    end
end
if ~foundbody, error('unable to find the tag office:body in file ''%s''. Is it a valid ODS file?',infile), end
if item<3, dispf('WARNING: non-standard ODS file, %d sections are missing but loadods is able to read it',3-item); end
node = nodes.item(item);
nodes = node.getChildNodes;
node = nodes.item(0);
nodes_sheets = node.getChildNodes;

% default sheetname (first sheet)
allsheets = false;
if isempty(options.sheetname)
    sheet = nodes_sheets.item(0);
    nodesinsheet = sheet.getChildNodes;
    if nodesinsheet.getLength==0 % work arround to remove table:calculation-settings (OV)
        sheet = nodes_sheets.item(1);
    end
    temp = get_attribute(sheet,'table:name'); % may return [] (OV: 16/09/11)
    if ~isempty(temp), options.sheetname = {temp}; end
elseif ischar(options.sheetname) && strcmp(options.sheetname,'all')
    allsheets = true;
end
if isempty(options.sheetname)
    dispf('\tWARNING: the name of the first sheet has not been determined, all sheets will be re read.\n\tSupply the name of the first sheet to remove this message.')
    allsheets = true;
end
if ~iscell(options.sheetname), options.sheetname = {options.sheetname}; end


% scan all sheets
numSheets = nodes_sheets.getLength;
ivalidsheet = 0;
dataout = [];
for isheet = 1:numSheets
    sheet = nodes_sheets.item(isheet-1);
    wksht = get_attribute(sheet,'table:name');
    nodesinsheet = sheet.getChildNodes;
    
    if ~isempty(wksht) && (nodesinsheet.getLength>0) && (allsheets || ismember(lower(wksht),lower(options.sheetname)))
        
        % convert sheetname into a valid field name
        ivalidsheet = ivalidsheet + 1;
        fs = regexprep(wksht,validchars,'');
        if ~isempty(regexp(fs,'^\d','once')), fs = sprintf('s%s',fs); end
        if isempty(fs), fs = sprintf('sheet%d',ivalidsheet); end

        % Get the number of columns
        nodes = sheet.getChildNodes;
        num_nodes = nodes.getLength;
        num_cols = 0;
        for count = 1:num_nodes
            node = nodes.item(count-1);
            if strcmp(char(node.getNodeName),'table:table-column')
                temp = get_attribute(node,'table:number-columns-repeated');
                if ~isempty(temp)
                    num_cols = num_cols+1; %num_cols = num_cols+str2num(char(num_cols));
                end
            elseif strcmp(char(node.getNodeName),'table:table-row')
                count = count-2; %#ok<FXSET>
                break
            end
        end
        
        % Get the number of rows
        num_rows = num_nodes-count-1;
        
        % Initialize memory for the data
        data = cell(num_rows,num_cols);
        
        % Extract the data for the sheet
        screen = ''; t0 = clock;
        for row_num = 1:num_rows
            screen = dispb(screen,'read row %d/%d in ''%s'' ...',row_num,num_rows,wksht);
            row = nodes.item(count+row_num);
            cols = row.getChildNodes;
            col_num = 0;
            num_items = cols.getLength-1;
            for item_num = 0:num_items
                col = cols.item(item_num);
                num_repeated = get_attribute(col,'table:number-columns-repeated');
                num_repeated = str2double(char(num_repeated));
                value_type = get_attribute(col,'office:value-type');
                if strcmp(value_type,'string')
                    temp = col.getChildNodes;
                    temp = temp.item(0);
                    temp = temp.getChildNodes;
                    temp = temp.item(0);
                    if any(strcmp(methods(temp),'getData'))
                        value = char(temp.getData);
                        switch lower(value)
                            case 'inf'
                                value = Inf;
                            case '-inf'
                                value = -Inf;
                            case 'nan'
                                value = NaN;
                            case 'true'
                                value = true;
                            case 'false'
                                value = false;
                        end
                    else
                        value = options.blank;
                    end
                elseif strcmp(value_type,'float')
                    value = str2double(get_attribute(col,'office:value'));
                else
                    value = options.blank;
                end
                if ~isempty(num_repeated)
                    for i = 1:num_repeated
                        col_num = col_num+1;
                        data{row_num,col_num} = value;
                    end
                else
                    col_num = col_num+1;
                    data{row_num,col_num} = value;
                end
            end
        end
        dispb(screen,'... read ''%s'' with %d rows in %0.4g s',wksht,num_rows,etime(clock,t0));
        
        % Fix the output
        % determine bounds
        % Before 17/01/2014
        % ok = reshape(cellfun(@(x) ischar(x) || (~isnan(x) && ~isempty(x)),data),size(data));
        % After 17/01/2014
        ok = reshape(cellfun(@(x) ischar(x) || (~all(isnan(x)) && ~isempty(x)),data),size(data));
        data = data(1:max(sum(ok,1)),1:max(sum(ok,2)));
        if options.transpose, data = data'; end
        if options.headers && ~isempty(data)
            headers = data(1:options.headers,:);
            data = data(options.headers+1:end,:);
        else
            headers = '';
        end
        dispf('LOADODS:\t''%s'' %dx%d worksheet with %d header lines',wksht,size(data),size(headers,1))
        n = size(data,2);
        % concatenate rows
        if options.concat && ~isempty(data)
            datatype = cellfun(@isnumeric,data(1,:));
            tmp = cell(1,n);
            for i =1:n
                if datatype(i) && all(cellfun(@isnumeric,data(:,i)))
                    tmp{i} = cat(1,data{:,i});
                    if options.forceboolean && all((tmp{i}==0) | (tmp{i}==1))
                        tmp{i} = logical(tmp{i});
                    end
                else
                    tmp{i} = data(:,i);
                end
            end
            if all(cellfun(@isnumeric,tmp)) && ~options.struct
                tmp = cat(2,tmp{:});
            end
            data = tmp;
        end
        % convert to struct
        if options.struct && ~isempty(data)
            if isempty(dataout), dataout = struct([]); end
            if isempty(headers), headers = repmat({''},1,n); end
            for i=1:n
               %if ischar(headers{:,i}), f = regexprep([headers{:,i}],validchars,''); else f = ''; end
                tmpheaders = [headers{cellfun(@(x) ischar(x),headers(:,i)),i}];
                if ischar(tmpheaders), f = regexprep(tmpheaders,validchars,''); else f = ''; end
                if ~isempty(regexp(f,'^\d','once'))
                    if options.transpose, f = sprintf('r%s',f); else f = sprintf('c%s',f); end
                end
                if isempty(f)
                    if options.transpose, f = sprintf('row%0.2d',i); else f = sprintf('col%0.2d',i); end
                end
                if length(f)>MaxlengthField, f = f(1:MaxlengthField); end
                if options.concat
                    dataout(1).(fs).(f) = data{i};
                else
                    dataout(1).(fs).(f) = data(:,i);
                end
            end
            if options.structarray
                currentfields = fieldnames(dataout(1).(fs));
                currentvalues = struct2cell(dataout(1).(fs));
                notcell = cellfun(@(x)~iscell(x),currentvalues);
                currentvalues(notcell) = cellfun(@(x) num2cell(x,2),currentvalues(notcell),'UniformOutput',false);
                tmp = [currentfields currentvalues]';
                dataout(1).(fs) = struct(tmp{:});
            end
        else
            if isempty(dataout), dataout = {}; end
            if ~isempty(data), dataout{end+1} = data; end %#ok<AGROW>
        end
        
    end % if ismember
end % next count

% remove worksheet information for single worksheetsheet
if ivalidsheet<1
    warning('no data in ''%s''',infile) %#ok<WNTAG>
    dataout = '';
elseif ivalidsheet==1
    if iscell(dataout), dataout = dataout{1};
    elseif isstruct(dataout), dataout = dataout.(fs);
    end
end

% Remove the temporary files
if ~rmdir(dir_temp,'s')
    warning('Temporary files in ''%s'' could not be removed',dir-temp); %#ok<WNTAG>
end


end % end function


%% PRIVATE FUNCTIONS

% Returns the value of the attribute_name in node
function attribute = get_attribute(node,attribute_name)
attribute = [];
if node.hasAttributes
    attributes = node.getAttributes;
    num_attributes = attributes.getLength;
    for count = 1:num_attributes
        item = attributes.item(count-1);
        if strcmp(char(item.getName),attribute_name)
            attribute = char(item.getValue);
        end
    end
end
end
