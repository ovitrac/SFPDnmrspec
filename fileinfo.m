function [rout,strout] = fileinfo(filename,ppath,dispon)
%FILEINFO returns a structure that contains the main information of a given file
%  Syntax: info = fileinfo(filename,[path])
%  Option: [info,str] = fileinfo(filename,[path])
%   filename can include a full path or not
%   path is used if filename does not include it (the current path is used otherwise)
%   info = structure
%   str = formated string

% MS-MATLAB 1.0 - 20/04/04 - INRA\Olivier Vitrac - rev. 12/01/11

% revision history
% 19/08/04: display if no ouput, added str
% 11/09/07 fix filename as cell array
% 12/01/11 fix versn in FILEPARTS for Matlab later than 7.11 (not supported any more)

% arg check
if nargin<2, ppath = []; end
if nargin<3, dispon = []; end
if ~nargout || isempty(dispon), dispon = true; end
oldmatlab = verLessThan('matlab','7.11');

if iscell(filename)
    m = length(filename);
    str = cell(m,1);
    [r,str{1}] = fileinfo(filename{1},ppath,dispon);
    for i=2:m, [r(i),str{i}] = fileinfo(filename{i},ppath,dispon); end
    dispon = false;
else
    % other checks
    if oldmatlab
        [pathstr,name,ext,versn] = fileparts(filename); %#ok<FPART>
    else
        [pathstr,name,ext] = fileparts(filename);
        versn = [];
    end
    if isempty(pathstr)
        if isempty(ppath)
            pathstr = cd;
        else
            pathstr = ppath;
        end
    end
    if ~exist(pathstr,'dir')
        error(sprintf('the directory ''%s'' does not exist',pathstr))
    end
    fullfilename = fullfile(pathstr,[name ext versn]);
    if ~exist(fullfilename,'file')
        error(sprintf('the file ''%s'' does not exist in ''%s''',[name ext versn],pathstr))
    end

    % info
    f = dir(fullfilename);

    r =struct(   'filename', [name ext versn], ...
        'name',     name, ...
        'ext',      ext, ...
        'ver',      versn, ...
        'date',     f.date, ...
        'bytes',    f.bytes, ...
        'path',     pathstr...
        );

    if f.bytes>1024*1024
        str = sprintf('\t%s%s\t\t%s\t\t%0.1f MBytes\t\t%s',name,ext,f.date,f.bytes/(1024*1024),pathstr);
    elseif f.bytes>1024
        str = sprintf('\t%s%s\t\t%s\t\t%0.1f kBytes\t\t%s',name,ext,f.date,f.bytes/1024,pathstr);
    else
        str = sprintf('\t%s%s\t\t%s\t\t%d Bytes\t\t%s',name,ext,f.date,f.bytes,pathstr);
    end
end

% output
if nargout>0, rout = r; end
if nargout>1, strout = str; end
if dispon, disp(str), end
