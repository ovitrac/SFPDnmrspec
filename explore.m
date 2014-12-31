function list = explore(format,pathstr,profondeur,abbreviate)
% EXPLORE look for files matching a given or a list of formats
% 		Syntax: filelist = explore(format, [path],[depth],['abbreviate' or 'fullabbreviate'])
%		format =  string ou cells of string (see isformat for details, very powerfull)
%                 format are not case sensitive (required on WIN machines)
%                 Main metacharacters based on the old syntax of isformat
%                       ?  = any character
%                       #  = any digit
%                       #* = any number (list of digits)
%                       ~  = not a digit
%                       @  = not a word character ~([a-z_A-Z0-9])              
%                       *  = any list of characters
%                       ==> any combinations and occurences are possible
%       path = any valid path (local or remote) or current directory
%		maxdepth = maximum of examined levels (par d�faut = 10)
%       'abbreviate' generate a simplified structure of filenames
%       'fullabbreviate' generate only a cell string of filenames
%
%       Conversion between formats is possible with:
%           filelist = explore(filelist) which is equivalent to filelist = explore(filelist,'fullabbreviate')
%           filelist = explore(filelist,'abbreviate');

% TCPIP v. 1.0 - 26/02/01 - INRA\Olivier Vitrac rev. 04/01/11

% Revision history
%   25/02/04 abbreviate format
%   01/03/04 subpath
%   14/05/04 fix bug when an empty list occures with abbreviate format
%   25/01/05 fix case profondeur <= 1
%   10/02/05 fix fieldnames ~pathstr -> ~path
%   06/03/05 add fields bytes and date
%   04/12/05 update to new formats (see isformat), add 'full abbreviate'
%   20/01/08 add post conversions
%   24/01/08 improve compatibility with Matlab 7.5
%   04/01/11 fix versn in FILEPARTS for Matlab later than 7.11 (not supported any more)
%   04/01/11 add datenum

% default
profondeur_default = 19;
found = false;
oldmatlab = verLessThan('matlab','7.11');

% arg check
if nargin<2, pathstr = []; end
if isstruct(format), list = abbreviateform(format,pathstr,oldmatlab); return, end
if nargin<3, profondeur = []; end
if nargin<4, abbreviate = 'none'; end
if isempty(pathstr), pathstr = cd; end
if isempty(profondeur), profondeur = profondeur_default; end
sizeofroot = length(pathstr);
currentfilesep = filesep;

% exploration
if profondeur <=1
    rawlist = dir(fullfile(pathstr,format));
    j = 0;
    list = repmat(struct('path','','file','','subpath','','name','','ext','','ver','','date','','datenum',[],'bytes',[]),0,0);
    for i = 1:length(rawlist)
        if ~rawlist(i).isdir
            j = j + 1;
            if oldmatlab
                [p,n,e,v] = fileparts(rawlist(i).name); %#ok<*FPART,ASGLU>
                date_num = datenum(rawlist(i).date);
            else
                [p,n,e] = fileparts(rawlist(i).name); %#ok<ASGLU>
                v = NaN;
                date_num = rawlist(i).datenum;
            end
            list(j) = struct(   'path', pathstr, ...
                                'file', rawlist(i).name, ...
                                'subpath','', ...
                                'name',n, ...
                                'ext',e, ...
                                'ver',v,...
                                'date',rawlist(i).date,...
                                'datenum',date_num,...
                                'bytes',rawlist(i).bytes);
            found = true;
        end                     
    end
    list = abbreviateform(list,abbreviate,oldmatlab);
else
    origine = cd; cd(pathstr) %#ok<*MCCD>
    i_prof = 1 ;
    fich = cell(profondeur,1);
    fich{i_prof} = dir;
    i_fich = zeros(profondeur,1);
    i_fich(i_prof) = 1;
    i_list = 1;
    list = [];
    while (i_prof <= profondeur) && (i_prof > 0)
        if i_fich(i_prof) > length(fich{i_prof})
            cd('..')
            i_prof = i_prof - 1;
        else
            format_match = isformat(fich{i_prof}(i_fich(i_prof)).name,format);
            if fich{i_prof}(i_fich(i_prof)).isdir && fich{i_prof}(i_fich(i_prof)).name(1) ~= '.'
                cd(fich{i_prof}(i_fich(i_prof)).name)
                i_prof = i_prof + 1;
                i_fich(i_prof) = 0;
                fich{i_prof} = dir;
            elseif any(format_match) && fich{i_prof}(i_fich(i_prof)).bytes > 0
                found = true;
                list(i_list).format_match = format_match; %#ok<*AGROW>
                list(i_list).info.path = cd;
                list(i_list).info.file = fich{i_prof}(i_fich(i_prof)).name;
                if length(list(i_list).info.path)>sizeofroot+1
                    list(i_list).info.subpath = list(i_list).info.path(sizeofroot+2:end);
                else
                    list(i_list).info.subpath = currentfilesep;
                end
                list(i_list).info.date = fich{i_prof}(i_fich(i_prof)).date;
                if oldmatlab
                    list(i_list).info.datenum = datenum(list(i_list).info.date);
                else
                    list(i_list).info.datenum = fich{i_prof}(i_fich(i_prof)).datenum;
                end
                list(i_list).info.bytes = fich{i_prof}(i_fich(i_prof)).bytes;
                i_list = i_list + 1;
            end
        end
        if i_prof>0, i_fich(i_prof) = i_fich(i_prof) + 1; end
    end
    cd(origine);

    % abbreviate
    if any(abbreviate), list = abbreviateform(list,abbreviate,oldmatlab); end
end

if ~found, list = []; end

function listout = abbreviateform(list,abbreviate,oldmatlab)
% ABBREVIATE FORM
if isempty(abbreviate), abbreviate  = 'fullabbreviate'; end
switch lower(abbreviate)
    case 'abbreviate'
        if isempty(list) || ~isfield(list,'info')
            listout=list;
        else
            listout = struct2cell(list);
            listout = [listout{2,:,:}];
            for i=1:length(listout)
                if oldmatlab
                    [p,n,e,v] = fileparts(listout(i).file); %#ok<ASGLU>
                else
                    [p,n,e] = fileparts(listout(i).file); %#ok<ASGLU>
                    v = NaN;
                end
                listout(i).name = n;
                listout(i).ext = e(2:end);
                listout(i).ver = v;
            end
        end
    case 'fullabbreviate'
        nlist = length(list);
        listout = cell(nlist,1);
        if isfield(list,'info')
            for i=1:nlist, listout{i} = fullfile(list(i).info.path,list(i).info.file); end
        elseif isfield(list,'file')
            for i=1:nlist, listout{i} = fullfile(list(i).path,list(i).file); end
        else
            listout = list;
        end
    otherwise
        listout=list;
end
% function list = explore(format,pathstr,profondeur,abbreviate)
% % EXPLORE look for files matching a given or a list of formats
% % 		Syntax: filelist = explore(format, [path],[depth],['abbreviate' or 'fullabbreviate'])
% %		format =  string ou cells of string (see isformat for details, very powerfull)
% %                 format are not case sensitive (required on WIN machines)
% %                 Main metacharacters based on the old syntax of isformat
% %                       ?  = any character
% %                       #  = any digit
% %                       #* = any number (list of digits)
% %                       ~  = not a digit
% %                       @  = not a word character ~([a-z_A-Z0-9])              
% %                       *  = any list of characters
% %                       ==> any combinations and occurences are possible
% %       path = any valid path (local or remote) or current directory
% %		maxdepth = maximum of examined levels (par d�faut = 10)
% %       'abbreviate' generate a simplified structure of filenames
% %       'fullabbreviate' generate only a cell string of filenames
% %
% %       Conversion between formats is possible with:
% %           filelist = explore(filelist) which is equivalent to filelist = explore(filelist,'fullabbreviate')
% %           filelist = explore(filelist,'abbreviate');
% 
% % TCPIP v. 1.0 - 26/02/01 - INRA\Olivier Vitrac rev. 06/05/05
% 
% % Revision history
% %   25/02/04 abbreviate format
% %   01/03/04 subpath
% %   14/05/04 fix bug when an empty list occures with abbreviate format
% %   25/01/05 fix case profondeur <= 1
% %   10/02/05 fix fieldnames ~pathstr -> ~path
% %   06/03/05 add fields bytes and date
% %   04/12/05 update to new formats (see isformat), add 'full abbreviate'
% %   20/01/08 add post conversions     
% 
% % default
% profondeur_default = 19;
% extsep = '.';
% found = false;
% 
% % arg check
% if nargin<2, pathstr = []; end
% if isstruct(format), list = abbreviateform(format,pathstr); return, end
% if nargin<3, profondeur = []; end
% if nargin<4, abbreviate = 'none'; end
% if isempty(pathstr), pathstr = cd; end
% if isempty(profondeur), profondeur = profondeur_default; end
% sizeofroot = length(pathstr);
% currentfilesep = filesep;
% 
% % exploration
% if profondeur <=1
%     rawlist = dir(fullfile(pathstr,format));
%     j = 0;
%     for i = 1:length(rawlist)
%         if ~rawlist(i).isdir
%             j = j + 1;
%             [p,n,e,v] = fileparts(rawlist(i).name);
%             list(j) = struct(   'path', pathstr, ...
%                                 'file', rawlist(i).name, ...
%                                 'subpath','', ...
%                                 'name',n, ...
%                                 'ext',e, ...
%                                 'ver',v,...
%                                 'date',rawlist(i).date,...
%                                 'bytes',rawlist(i).bytes);
%             found = true;
%         end                     
%     end
%     list = abbreviateform(list,abbreviate);
% else
%     origine = cd; cd(pathstr)
%     i_prof = 1 ; pos = zeros(0,1);
%     fich = cell(profondeur,1);
%     fich{i_prof} = dir;
%     i_fich = zeros(profondeur,1);
%     i_fich(i_prof) = 1;
%     i_list = 1;
%     list = [];
%     while (i_prof <= profondeur) & (i_prof > 0)
%         if i_fich(i_prof) > length(fich{i_prof})
%             cd('..')
%             i_prof = i_prof - 1;
%         else
%             format_match = isformat(fich{i_prof}(i_fich(i_prof)).name,format);
%             if fich{i_prof}(i_fich(i_prof)).isdir & fich{i_prof}(i_fich(i_prof)).name(1) ~= '.'
%                 cd(fich{i_prof}(i_fich(i_prof)).name)
%                 i_prof = i_prof + 1;
%                 i_fich(i_prof) = 0;
%                 fich{i_prof} = dir;
%             elseif any(format_match) & fich{i_prof}(i_fich(i_prof)).bytes > 0
%                 found = true;
%                 list(i_list).format_match = format_match;
%                 list(i_list).info.path = cd;
%                 list(i_list).info.file = fich{i_prof}(i_fich(i_prof)).name;
%                 if length(list(i_list).info.path)>sizeofroot+1
%                     list(i_list).info.subpath = list(i_list).info.path(sizeofroot+2:end);
%                 else
%                     list(i_list).info.subpath = currentfilesep;
%                 end
%                 list(i_list).info.date = fich{i_prof}(i_fich(i_prof)).date;
%                 list(i_list).info.bytes = fich{i_prof}(i_fich(i_prof)).bytes;
%                 i_list = i_list + 1;
%             end
%         end
%         if i_prof>0, i_fich(i_prof) = i_fich(i_prof) + 1; end
%     end
%     cd(origine);
% 
%     % abbreviate
%     if any(abbreviate), list = abbreviateform(list,abbreviate); end
% end
% 
% if ~found, list = []; end
% 
% function listout = abbreviateform(list,abbreviate)
% % ABBREVIATE FORM
% if isempty(abbreviate), abbreviate  = 'fullabbreviate'; end
% switch lower(abbreviate)
%     case 'abbreviate'
%         if isempty(list) || ~isfield(list,'info')
%             listout=list;
%         else
%             listout = struct2cell(list);
%             listout = [listout{2,:,:}];
%             for i=1:length(listout)
%                 [p,n,e,v] = fileparts(listout(i).file);
%                 listout(i).name = n;
%                 listout(i).ext = e(2:end);
%                 listout(i).ver = v;
%             end
%         end
%     case 'fullabbreviate'
%         nlist = length(list);
%         listout = cell(nlist,1);
%         if isfield(list,'info')
%             for i=1:nlist, listout{i} = fullfile(list(i).info.path,list(i).info.file); end
%         elseif isfield(list,'file')
%             for i=1:nlist, listout{i} = fullfile(list(i).path,list(i).file); end
%         else
%             listout = list;
%         end
%     otherwise
%         listout=list;
% end