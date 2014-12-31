function checktoolboxinstall(TOOLBOX_LOCAL_NAME,localpath)
%  CHECKTOOLBOXINSTALL check if a given is correctly installed (if not, try to install it temporarly)
%   Syntax: checktoolboxinstall(TOOLBOX [,localpath])

% QSPR 1.0 - 13/02/04 - Olivier Vitrac - rev. 24/11/14

% Revsions
% 16/02/04 remove case sensitivity
% 24/11/14 some code refresh, add localpath

currentfolder = pwd;
if nargin<2, localpath = ''; end
if isempty(localpath), localpath = currentfolder; end

TOOLBOX_LOCAL_PATH = find_path_toolbox(TOOLBOX_LOCAL_NAME);
OK = 0;
% try
    if isempty(TOOLBOX_LOCAL_PATH)
        [fold,root] = lastdir(localpath);
        if exist(TOOLBOX_LOCAL_NAME,'dir') % try actual and go down
            if strcmpi(fold,TOOLBOX_LOCAL_NAME) % try actual position
                addpath(cd)
                OK = 2;
            else
                addpath([cd filesep TOOLBOX_LOCAL_NAME])
                cd(TOOLBOX_LOCAL_NAME) % go down
                OK = 3;
            end
        else % try up
            if exist(root,'dir')
                addpath(root)
                cd('..') % go up
                OK = 4;
            end
        end
    else
        OK = 1;
    end
% end

if ~OK
    error('The toolbox ''%s'' is not installed (add it with addpath command)',TOOLBOX_LOCAL_NAME)
elseif OK>1
    dispf('The toolbox ''%s'' has been temporarly added to the current path',TOOLBOX_LOCAL_NAME)
end
cd(currentfolder)
