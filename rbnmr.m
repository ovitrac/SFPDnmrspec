function A = rbnmr(D,silent,search)
% RBNMR     Reads processed Bruker NMR-data.
%
% SYNTAX    A = rbnmr;			% Reads 1r/2rr in the current working dir.
%			A = rbnmr;			% Finds 1r/2rr files (rec.) in pdata/1 
%			A = rbnmr(DS);		% Reads the files in the cell DS
%			A = rbnmr(file);	% Reads the data specified in the text file
%			A = rbnmr('<Datadirectory>/?0')
%								% Finds 1r/2rr-files with expno 
%								%	10, 20, ..., 90 in pdata/1
%			A = rbnmr('<Datadirectory>/?0/pdata/999')
%								% Finds 1r/2rr-files with expno 
%								%	10, 20, ..., 90 in pdata/999
%
% OUT       A: Struct with nmrdata:

% Nils Nyberg, SLU, 2001-05-02
% rewritten 2007-09-19, Matlab rel. R2007a
% partly rewritten 2008-09-18, rel. R2008a
%		New perl-based search function
%		New parameter checking and organization of functionsz
% 2008-11-13, Empty 'title'-file allowed
% 2009-08-12, New option 'search': Avoid 'search' on UNC pathnames
% 2009-09-29, Robust regexp ('Anything'-case) to handle strange parameter files
% 2009-09-29, ME-exception returned as 'title'.
% 2010-08-05, Robust regexp ('EmptyPar'-case) to handle parameter files with 
% whitespaces changed(?) into line breaks.
% 2011-05-18, 'Strange title-file:' returned as Title-field if textscan
% fails to read the title-file.
% 2011-06-01, Changed regexp to allow the new more relaxed data path
%
% Nils Nyberg, Copenhagen University, nn@sund.ku.dk

%% Init values
CurrentWorkingDir = pwd;


%% Check inputs
if (nargin < 3);
	search = 1;
	% Warning: DOS programs may not execute correctly when the current directory is a UNC pathname.
	if all(ismember([1 2], regexp(CurrentWorkingDir,'\\')))
		search = 0;
	end
end
if (nargin < 2); silent = 0; end
if (nargin < 1); D = []; end

if (iscell(D));
	input_class = 'cell';
elseif (isempty(D) && exist('1r','file')==2)
	input_class = 'datafiledirect'; File = '1r';
elseif (isempty(D) && exist('2rr','file')==2)
	input_class = 'datafiledirect'; File = '2rr';
elseif (~isempty(D) && ~isempty(regexp(D,'[12][ri]{1,2}$', 'once' )) && exist(D,'file')==2)
	input_class = 'datafile';
elseif (exist(D,'file')==2);
	input_class = 'textfile';
elseif (isempty(D))
	input_class = 'search'; D = pwd;
else
	input_class = 'search';
end

switch input_class
	case 'cell'
		% D is a cell with one full path to the processed data directory (or a file
		% in this directory) per cell.
		% Change to that directory and do a recursive call
		A = cell(length(D),1);
		for i=1:length(D);
			try
				A{i} = rbnmr(D{i});
			catch ME
				A{i}.Title = ['File/Data directory does not exist: ',D{i}];
			end
			if ~silent; fprintf('%3d/%-3d %-45s %s\n',...
					i,length(D),D{i},A{i}.Title);
			end
		end
		% Fix output, cell is only used when more than 
		% one data set is returned...
		if (length(A) == 1); A = A{1}; end
case 'textfile'
	% Read in the text file with one path per row.
	% Do recursive calls
	[fid, message] = fopen(D,'rt');
	if ~fid ; error(message); end
	DS = textscan(fid,'%s','CommentStyle','%'); DS = DS{1};
	fclose(fid);
	A = rbnmr(DS,silent);
case 'datafile'
	try cd(fileparts(D))
	catch ME
		A.Title = ME.message;
		cd(CurrentWorkingDir);
		return;
	end
	A = rbnmr([],silent); cd(CurrentWorkingDir);
case 'datafiledirect'
	try
		A = do_the_actual_import(File);
	catch ME
		A.Title = ME.message;
		cd(CurrentWorkingDir);
		return;
	end
case 'search'
	if search
		DS = search_data(D,silent);
	else
		DS = [];
	end
	if isempty(DS);
		cd(CurrentWorkingDir);
		error('RBNMR: Can not find default files (1r or 2rr) in %s or below.',...
			D);
	end

	% Ask if too many files
	reply = 'y';
	if ~silent && length(DS)>50
		reply = input(sprintf(...
			'%d NMR-files found. Do you want to continue? [y]',length(DS)),'s');
		if isempty(reply); reply = 'y'; end
	end
	if lower(reply)=='y'
		if ~silent;	disp('OK! Here we go...');end
		A = rbnmr(DS,silent);
	else
		fprintf('Yes, you are right. We''d better skip this... (%s)\n',reply);
		A=[];
	end
end

% Fix out put and clean up
cd(CurrentWorkingDir);
return








function A = do_the_actual_import(File)

%% Read first line of title-file
fid = fopen('title','r');
if (fid ~= -1);
        title = textscan(fid,'%s','Whitespace','\n','ReturnOnError',1);
        fclose(fid);
		if isempty(title{1})
			A.Title = '<Title file empty>';
		else
			try
				A.Title = title{1}{1};
			catch ME
				A.Title = sprintf('Strange title-file: %s',ME.message);
			end
		end
else
	A.Title = '<Title file empty>';
end
	

%% Date and file information
[path,name,ext] = fileparts(File);
if isempty(path); path = pwd; end;
A.Info.ImportDate = datestr(now);
A.Info.FileName = [name ext];
A.Info.FilePath = path;
A.Info.Title = A.Title;


%% Add relative path from 'path'
q = regexpi(...
    fullfile(A.Info.FilePath,A.Info.FileName),...
    'data[/\\].+[/\\]nmr[/\\](.+)[/\\](\d+)[/\\]pdata[/\\](\d+)[/\\](.+)','tokens');
s = '/';	% unix-style
try
	A.Info.RelativePath = [q{1}{1},s,q{1}{2},s,'pdata',s,q{1}{3},s,q{1}{4}];
catch ME
	% Do nothing. Relative path does not seem to make any sense...
end


%% Check parameter files and read parameters
if strcmp(File,'2rr') && exist('proc2s','file')==2;
    A.Proc2s = readnmrpar('proc2s');
end
if strcmp(File,'2rr') && exist('../../acqu2s','file')==2;
    A.Acqu2s = readnmrpar('../../acqu2s');
end
if exist('../../acqus','file')==2;
    A.Acqus = readnmrpar('../../acqus');
else
	error('RBNMR: Could not find ../../acqus')
end
if exist('procs','file')==2;
    A.Procs = readnmrpar('procs');
end

%% Add acq-date
% Converts time given in UTC (base 1970, seconds) as matlab serial time
% (base 0000, days)
TZ = str2double(regexp(A.Acqus.Stamp,'UT(.\d+)h','tokens','once'));
if isempty(TZ); TZ = 2; end;	% Assume UT+2h if not in stamp-field
A.Info.AcqSerialDate = A.Acqus.DATE/(60*60*24)+datenum([1970 01 01])+TZ/24;
A.Info.AcqDateTime = datestr(A.Info.AcqSerialDate);
A.Info.AcqDate = datestr(A.Info.AcqSerialDate,'yyyy-mm-dd');
% Convert serial date to text to keep format
A.Info.AcqSerialDate = sprintf('%.12f',A.Info.AcqSerialDate);

%% Add plotlabel from A.Acqus.Stamp-info
q = regexp(A.Acqus.Stamp,'data[/\\].+[/\\]nmr[/\\](.+)[/\\](\d+)[/\\]acqus','tokens');
if isempty(q)	% New, more relaxed, data path
	q = regexp(A.Acqus.Stamp,'#.+[/\\](.+)[/\\](\d+)[/\\]acqus','tokens');
end
if isempty(q)
	A.Info.PlotLabel = ['[',A.Info.FilePath,']'];
else
	A.Info.PlotLabel = ['[',q{1}{1},':',q{1}{2},']'];
end

%% Open and read file
if A.Procs.BYTORDP == 0
    endian = 'l';
else
    endian = 'b';
end

[FID, MESSAGE] = fopen(File,'r',endian);
if FID == -1
	disp(MESSAGE);
	error(['RBNMR: Error opening file (',File,').']);
end

A.Data = fread(FID,'int32');
fclose(FID);

%% Read imaginary data if the file 1i exists
if (exist('1i','file')==2)
    [FID, MESSAGE] = fopen('1i','r',endian);
    if FID == -1
        % Do nothing
    end
    A.IData = fread(FID,'int32');
    fclose(FID);
end    

%% Correct data for NC_proc-parameter
A.Data = A.Data/(2^-A.Procs.NC_proc);
if (isfield(A,'IData'))
    A.IData = A.IData/(2^-A.Procs.NC_proc);
end

A.Procs.NC_proc = 0;

%% Calculate x-axis
A.XAxis = linspace( A.Procs.OFFSET,...
                    A.Procs.OFFSET-A.Procs.SW_p./A.Procs.SF,...
                    A.Procs.SI)';

%% Additional axis and reordering if 2D-file
if isfield(A,'Proc2s')
    A.YAxis = linspace( A.Proc2s.OFFSET,...
                        A.Proc2s.OFFSET-A.Proc2s.SW_p./A.Proc2s.SF,...
                        A.Proc2s.SI)';

%% Reorder submatrixes (se XWinNMR-manual, chapter 17.5 (95.3))

		SI1 = A.Procs.SI; SI2 = A.Proc2s.SI;
		XDIM1 = A.Procs.XDIM; XDIM2 = A.Proc2s.XDIM;

		NoSM = SI1*SI2/(XDIM1*XDIM2);    % Total number of Submatrixes
		NoSM2 = SI2/XDIM2;		 			% No of SM along F1

		A.Data = reshape(...
				permute(...
					reshape(...
						permute(...
							reshape(A.Data,XDIM1,XDIM2,NoSM),...
						[2 1 3]),...
					XDIM2,SI1,NoSM2),...
				[2 1 3]),...
  			    SI1,SI2)';

%% Read the level file if it exists
% The old version (level) is a binary 
	if(exist('level','file')==2)
		[FID, MESSAGE] = fopen('level','r',endian);
		if FID == -1
			disp('READBNMR: Error opening level file');
			disp(MESSAGE);
		end

		L=fread(FID,'int32');
		fclose(FID);

		% The first two figures is the number of pos. and neg. levels
		A.Levels = L(3:end);
		% Adjust for NC-parameter
		A.Levels = A.Levels/(2^-A.Procs.NC_proc);
	end

%% Read the clevel file if it exists
% The new version (clevel) is a text file 
	if(exist('clevels','file')==2)
		L = readnmrpar('clevels');
		switch L.LEVSIGN
			case 0	% Positive only
				A.Levels = L.LEVELS(L.LEVELS > 0)';
			case 1	% Negative only
				A.Levels = L.LEVELS(L.LEVELS < 0)';
			case 2	% Both pos and neg.
				A.Levels = L.LEVELS(1:L.MAXLEV*2);
		end
	end

%% Check that A.Levels is not one (large) scalar. If so 'Contour' will crash.
	if (isfield(A,'Levels') && length(A.Levels) == 1)
		A.Levels = [A.Levels;A.Levels];
	end

end


function P = readnmrpar(FileName)
% RBNMRPAR      Reads BRUKER parameter files to a struct
%
% SYNTAX        P = readnmrpar(FileName);
%
% IN            FileName:	Name of parameterfile, e.g., acqus
%
% OUT           Structure array with parameter/value-pairs
%

% Read file
A = textread(FileName,'%s','whitespace','\n');

% Det. the kind of entry
TypeOfRow = cell(length(A),2);
    
R = {   ...
    '^##\$*(.+)=\ \(\d\.\.\d+\)(.+)', 'ParVecVal' ; ...
    '^##\$*(.+)=\ \(\d\.\.\d+\)$'   , 'ParVec'    ; ...
    '^##\$*(.+)=\ (.+)'             , 'ParVal'    ; ...
    '^([^\$#].*)'                   , 'Val'       ; ...
    '^\$\$(.*)'                     , 'Stamp'     ; ...
    '^##\$*(.+)='                   , 'EmptyPar'  ; ...
	'^(.+)'							, 'Anything'	...
    };

for i = 1:length(A)
    for j=1:size(R,1)
        [s,t]=regexp(A{i},R{j,1},'start','tokens');
        if (~isempty(s))
            TypeOfRow{i,1}=R{j,2};
            TypeOfRow{i,2}=t{1};
        break;
        end
    end
end

% Set up the struct
i=0;
while i < length(TypeOfRow)
    i=i+1;
    switch TypeOfRow{i,1}
        case 'ParVal'
            LastParameterName = TypeOfRow{i,2}{1};
            P.(LastParameterName)=TypeOfRow{i,2}{2};
        case {'ParVec','EmptyPar'}
            LastParameterName = TypeOfRow{i,2}{1};
            P.(LastParameterName)=[];
        case 'ParVecVal'
            LastParameterName = TypeOfRow{i,2}{1};
            P.(LastParameterName)=TypeOfRow{i,2}{2};
        case 'Stamp'
            if ~isfield(P,'Stamp') 
                P.Stamp=TypeOfRow{i,2}{1};
            else
                P.Stamp=[P.Stamp ' ## ' TypeOfRow{i,2}{1}];
            end
        case 'Val'
			if isempty(P.(LastParameterName))
				P.(LastParameterName) = TypeOfRow{i,2}{1};
			else
				P.(LastParameterName) = [P.(LastParameterName),' ',TypeOfRow{i,2}{1}];
			end
        case {'Empty','Anything'}
            % Do nothing
    end
end
    

% Convert strings to values
Fields = fieldnames(P);

for i=1:length(Fields);
    trystring = sprintf('P.%s = [%s];',Fields{i},P.(Fields{i}));
    try
        eval(trystring);
	catch %#ok<CTCH>
        % Let the string P.(Fields{i}) be unaltered
    end
end

function DS = search_data(pathRE,silent) %#ok<INUSD>
% Try to recursively find 1r or 2rr files
% The last part of the argument are used as search criterea
% in the expnumber part of the data path.

% The procno is restricted to procno=1

% The search is performed by a perl function (for portability).

if nargin < 2; silent = 0; end

path = pathRE;
while ~isempty(path) && not(exist(path,'dir')==7);
	path = fileparts(path);
end

	
if isempty(path); 
	path = '.'; RE = pathRE; 
else
	RE = pathRE(length(path)+2:end); %Empty if too short
end

if strcmp(path(end),filesep) % Double backslash belongs to RE
	path = path(1:end-1);
	if ~isempty(RE)
		RE = [filesep,RE];
	end
end

% Make the matching more 'dos'-like...
% ? => .?, * => .+
RE = regexprep(RE,'?','.');
RE = regexprep(RE,'*','.+');



try cd(path)
catch ME
	rethrow(ME);
end

% Fix for running GNU perl under windows
if ispc; my_filesep = '/'; else my_filesep = filesep; end

% Huristics for making the search pattern...
PerlStr1 = 'use File::Find;use vars qw/*name/;*name=*File::Find::name;';
if strcmp(path,'.') && isempty(regexp(RE,'pdata', 'once' ))
	PerlStr2 = sprintf(...
		'find( sub {print $name."\\n" if $name =~ m!%s%s%spdata%s1%s[12]r{1,2}$!},".");',...
		my_filesep,RE,my_filesep,my_filesep,my_filesep);
elseif ~isempty(regexp(RE,'pdata', 'once' ))
	PerlStr2 = sprintf(...
		'find( sub {print $name."\\n" if $name =~ m!%s%s%s[12]r{1,2}$!},".");',...
		my_filesep,RE,my_filesep);
else
	PerlStr2 = 'find( sub {print $name."\n" if m![12]r{1,2}$!},".");';
end
PerlFile = tempname;
fid = fopen(PerlFile,'wt');
fprintf(fid,'%s\n',PerlStr1);
fprintf(fid,'%s\n',PerlStr2);
fclose(fid);

% if ~silent
% 	disp('Search')
% 	disp('======')
% 	disp(sprintf('Search for the data files: path = %s, search expression = %s',path,RE));
% 	disp('')
% 	disp('Perl program:')
% 	type(PerlFile);
% 	disp('=============================================================')
% end

PerlStr = ['perl "' PerlFile '"'];
if ispc
	PerlCmd = fullfile(matlabroot, 'sys\perl\win32\bin\');
	PerlCmd = ['set PATH=',PerlCmd, ';%PATH%& ' PerlStr];
	[status, result] = dos(PerlCmd);
else
	[status, result] = unix(PerlStr);
end
if status; error('RBNMR: can not find/execute perl'); end
delete(PerlFile);


try DS = textscan(result,'%s');
catch ME
	DS = [];
	return;
end

DS = sortDS(DS{1});

function DS = sortDS(SortThis)
DS = SortThis;
my_filesep = '/';

for i=1:6
	re = sprintf('%s(\\d{%d})%s',my_filesep,i,my_filesep);
	SortThis = regexprep(SortThis,re,'/0$1/');
end

[tmp,I]=sort(SortThis); %#ok<ASGLU>
DS = DS(I);


 
