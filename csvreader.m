function outcell = csvreader(csvfile, varargin)
% return a dataset or cell of a delimited file using the robust opencsv java
% library.
%
% NOTE: you must run 'installCSVReader.m' first (or otherwise have the
% opencsv JAR on your MATLAB java class path)
%
% d = CSVREADER(csvfile) returns the contents of a csv file with 1 row of
% variable names as a dataset array
%
% d = CSVREADER(csvfile, 'Parameter', Value) passes parameter/value
% pairs to the CSVREADER. The following are valid parameters:
%
%   * 'delimiter'       - single character delimiter (default = ',')
%   * 'quotechar'       - single character to treat as quote character 
%                         (default = '"')
%   * 'strictQuotes'    - when true, ignore chars outside quotes (default =
%                         false)
%   * 'skiplines'       - number of lines at beginning of file to
%                         completely skip (NOT column heading lines).
%                         (default = 0)
%   * 'HeaderLines'     - number of lines after skip lines that contain column 
%                         heading names. Uses first line of header lines to
%                         name columns if output as a dataset.
%                         (default = 1)
%   * 'escapeChar'      - single character that escapes delimiter and quotes
%                         (default = '\')
%   * 'ignoreLeadingWhiteSpace' 
%                       - if true, parser ignores white space before field
%                         (default = true)
%   * 'dataset'         - if true, tries to convert from strings columns 
%                         that are numeric or logical and then returns a dataset.
%						  If false, the contents of the CSV are returned as a 
%						  cell array of all strings. (default = true)
%   
% Tips:
%   * All rows in the CSV must have the same number of columns
%   * If 'dataset' parameter is false, the CSV contents are returned as a
%     cell array.
%   * This file is slower than other built-in MATLAB functionality (e.g.
%     csvread or dlmread), but is more robust in its treatment of quoted
%     strings.  For example, the dataset constructor will fail if you have
%     quoted strings with a delimiter in the quotes.
%
% Examples
% --------
% % read in a normal CSV with 1 row of header lines as a dataset
% d = csvreader('test.csv')
%
% % read a tab delimited file instead
% d = csvreader('tabdelim.txt', 'delimiter', '\t'); 
%
% % return as a cell array (easier for messy CSVs without column headings)
% outcell = csvreader('messy.csv', 'dataset', false);
%
% See also: csvread
% 
% Copyright Kristofer D. Kusano, 9/13/2013
% Revised INRA\Olivier Vitrac - 24/10/2014
%
% This project uses the opencsv java library. The opencsv project
% website is http://opencsv.sourceforge.net.

%% Check toolbox installation
csvjavalib = 'opencsv-2.3.jar';
sourcecsvjavalib = 'http://www.mathworks.fr/matlabcentral/fileexchange/43551-csvreader';
if isempty(find_path_toolbox(csvjavalib,true))
    MSpath = find_path_toolbox('MS');
    if exist(fullfile(MSpath,csvjavalib),'file')
        warning('Force an installation of ''opencsv-2.3.jar'' from ''MS'' toolbox')
        installCSVReader(MSpath)
        dispf('Installation completed.\n')
    else
        error('please copy ''%s'' into ''%s'' from %s',csvjavalib,find_path_toolbox('MS'),sourcecsvjavalib)
    end
end

%% Parse Input
escape_ascii = char(0); % TODO: look up which ascii char is escape

p = inputParser; % built-in MATLAB input parser object

addRequired(p, 'csvfile', @ischar);

% name/val pairs (with defaults)
issinglechar = @(x) ~isempty(x) & ischar(x) & length(x) == 1 | x(1) == '\';
addParameter(p, 'delimiter', ',', issinglechar); % delimiter - doesn't need to be comma!
addParameter(p, 'quotechar', '"', issinglechar); % quote character
addParameter(p, 'strictQuotes', false, @islogical); % when true, ignore chars outside quotes
addParameter(p, 'skiplines', 0, @isnumeric); % number of lines to skip completely when reading file (NOT column headings)
addParameter(p, 'HeaderLines', 1, @isnumeric); % number of lines after skip lines that contain column headings. 
addParameter(p, 'escapeChar', escape_ascii, issinglechar); % char for escaping delimiter or quote 
addParameter(p, 'ignoreLeadingWhiteSpace', true, @islogical); % if true, parser ignores white space before field
addParameter(p, 'dataset', true, @islogical); % if true, tries to convert from strings columns that are numeric or logical

% run parser
parse(p, csvfile, varargin{:})
pr = p.Results; % the results

% fix potential escaped chars
pr.delimiter = sprintf(pr.delimiter);
pr.quotechar = sprintf(pr.quotechar);
if (strcmp(pr.escapeChar, '\'))
    pr.escapeChar = '\\';
end
pr.escapeChar = sprintf(pr.escapeChar);
%% Check JAR is on java class path
% cp = javaclasspath('-all');
% open_jar = csvjavalib; %'opencsv-2.3.jar'; % name of jar
% if (all(cellfun(@isempty, strfind(cp, open_jar))))
%     error('CSVReader:JARNotFound',...
%         '%s is not on the java class path. Have you run "installCSVReader" script yet?',...
%         open_jar)
% end
%% Open File
import au.com.bytecode.opencsv.CSVReader; % import java classes
import java.io.FileNotFoundException;

try
    reader = CSVReader(java.io.FileReader(csvfile),...
        pr.delimiter,...
        pr.quotechar,...
        pr.escapeChar,...
        int32(pr.skiplines),...
        pr.strictQuotes,...
        pr.ignoreLeadingWhiteSpace);
catch me
    error('csvreader:FileNotFound', 'The file %s does not exist', csvfile)
end

%% reall all lines in, convert to matlab cells
allLinesList = reader.readAll(); % read all lines
allLinesCell = cell(allLinesList.toArray()); % convert to array of Ojbects, to cell
try
    outcell = [allLinesCell{:}]'; % expand and rotate
catch me
    if (strcmp(me.identifier, 'MATLAB:catenate:dimensionMismatch'))
        error('csvreader:RowsNotSame', 'number of columns in each row are not the same')
    end
end

%% Try to guess data type of each column
if (~pr.dataset) % go no further if we don't want to find columns
    return
end

outcell_backup = outcell; % backup copy incase something goes wrong
try
    % first, look at the first row
    istart = pr.HeaderLines + 1;
    firstLine = outcell(istart, :);

    % next convert numeric 
    cols = size(firstLine, 2);
    cols_num = ~isnan(str2double(firstLine));
    for i = 1:cols
        if (cols_num(i))
            % conver string to double
            outcell(istart:end, i) = num2cell(str2double(outcell(istart:end, i)));
            
            % check for logical
            if (all([outcell{istart:end, i}] == 1 | [outcell{istart:end, i}] == 0))
                outcell(istart:end, i) = num2cell(logical(str2double(outcell(istart:end, i))));
            end
        end
    end
    
    % decide column headings
    if (pr.HeaderLines == 0)
        varnames = arrayfun(@(x) sprintf('Var%d', x), 1:cols);
    else
        varnames = outcell(1, :);
    end
    outcell = dataset({outcell(istart:end, :), varnames{:}});
catch me
    warning('csvreader:ParseFailed',...
        'Failed to guess the data types of the table: %s', me.message)
    outcell = outcell_backup; % ditch and return cell instead
end