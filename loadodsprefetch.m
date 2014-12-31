function [data,isupdated]=loadodsprefetch(filename,varargin)
%LOADODSPREFETCH loadods surrogate to use/manage prefetch files when they exist
%      data = loadodsprefetch(...)
%         It uses the same syntax as loadods.
%         Additional properties/values are:
%           prefetchprefix: prefetch extension (default = 'PREFETCH_')
%             prefetchpath: path of the prefetch file (default=tempdir)            
%               noprefetch: flag to forcce the prefetch to be updated (default=false);
%       [data,isupdated]=loadodsprefetch(...)
%           isupdated returns true if the data has been updated

% MS 2.1 - 20/01/12 - INRA\Olivier Vitrac rev. 17/01/14

% Revision history
% 24/01/12 add a comparison based on requested sheetnames
% 26/01/12 use 'case' with argcheck to keep case in folder and file names.
% 26/01/12 fix non-empty loadodsoptions (as a list instead of a structure)
% 12/06/12 fix error message file missing
% 20/09/13 add isupdated
% 08/12/13 force prefetchupdate = true when new spreadsheets are requested
% 17/01/14 force columnwise loadodsoptions

% default
default = struct(...
    'prefetchprefix','PREFETCH_',...
    'prefetchpath',tempdir,...
    'noprefetch',false,...
    'sheetname',[] ...
    );

% arg check
if nargin<1, error('one argument is at least required'), end
[options,loadodsoptions] = argcheck(varargin,default,'','case');
if isempty(options.sheetname), options.sheetname=''; end
if ~isempty(options.sheetname), loadodsoptions(end+1:end+2) = {'sheetname';options.sheetname}; end % %propagate sheetname
loadodsoptions = loadodsoptions(:);
[~,prefetchfile] = fileparts(filename);
if ~exist(filename,'file'), error('the supplied file ''%s'' does not exist',filename); end
prefetchfile = fullfile(options.prefetchpath,[options.prefetchprefix prefetchfile '.mat']);
prefetchupdate = false;

% check that prefetch is up to date or match current needs
useprefetch = false; % default behavior
if options.noprefetch
    dispf('LOADODSPREFETCH: noprefetch option is used. The prefetch file is not used')
    prefetchupdate = true;
elseif ~exist(prefetchfile,'file')
    dispf('LOADODSPREFETCH: no prefetchfile detected')
    prefetchupdate = true;
else
    ref = dir(filename);
    pre = dir(prefetchfile);
    load(prefetchfile,'nfo');
    if (ref.datenum<pre.datenum) && (ref.datenum==nfo.datenum) && (ref.bytes==nfo.bytes) %#ok<NODEF> % prefetch up-to-date and same size
        load(prefetchfile,'sheetname');
        if ischar(sheetname) %#ok<NODEF>
            if strcmp(sheetname,'all')
                tmp=load(prefetchfile,'data');
                if isstruct(tmp.data), 
                    sheetname = fieldnames(tmp.data);
                    if ischar(options.sheetname) && strcmp(options.sheetname,'all'), options.sheetname = sheetname; end
                else
                    sheetname = {sheetname};
                end
            else
                sheetname = {sheetname};
            end
        end 
        if ischar(options.sheetname), options.sheetname = {options.sheetname}; end
        if isempty(setdiff(options.sheetname,sheetname))
            useprefetch = true;
        else
            dispf('LOADODSPREFETCH: prefetchfile is not used as additional worksheets are asked to be loaded.')
            prefetchupdate = true;
        end
    else
        dispf('LOADODSPREFETCH: prefetchfile is obsolete')
        prefetchupdate = true;
    end
end

% load data
if useprefetch
    dispf('LOADODSPREFETCH: use the prefetchfile below')
    fileinfo(prefetchfile)
    load(prefetchfile,'data') % load data
    if length(sheetname)>1 || strcmp(sheetname,'all') 
        data = rmfield(data,setdiff(sheetname,options.sheetname)); %#ok<NODEF> % remove unwanted worksheets
    end
    fdata = fieldnames(data);
    if length(fdata)==1 && isstruct(data.(fdata{1})), data = data.(fdata{1}); end
else
    loadodsoptions = [{filename};loadodsoptions];
    data = loadods(loadodsoptions{:});
end

% update prefetch file if needed
if prefetchupdate
    nfo = dir(filename); %#ok<NASGU>
    sheetname = options.sheetname; %#ok<NASGU>
    save(prefetchfile,'data','nfo','sheetname')
    dispf('LOADODSPREFETCH: the prefetch file below has been updated')
    fileinfo(prefetchfile)
end

if nargout>1, isupdated=prefetchupdate; end