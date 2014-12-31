function SFPDnmr(varargin)
%SFPDnmr launcher of application for deformulation of mixtures by 1H NMR in SFPD project
% Syntax
%       SFPDnmr()
%       SFPDnmr(apppath)
%       SFPDnmr(appconffile)
%       SFPDnmr(localpath,'property1',value1,'property2',value2,...)
%       SFPDnmr(localpath,'property1',value1,'property2',value2,...)
% Properties
%       'version': application version, (default=0.1)
%           'ext': common extensions for config files (default=SFPDnmr.xml)
%   'prefetchext': extension as above for prefetch files (default=SFPDnmr.mat)
%     'appconfig': default filename of appconfig (default=application) 
%     'nmrdbfile': default spectral database filename(default=nmrdb.mat) 
%  'deconvmethod': default deconvolution method filename (default=deconvolution_method)
%    'outputpath':  output folder (default=output)
%
% FOR USERS
% This application is delivered with some configuration files:
%   - appconfig: applcation configuration file loacated in 'apppath\private\' (default name 'application.SFPDnmr.xml')
%   - methodfile:  methodology/deformulation configuration file loacated in 'apppath\database\'  (default name 'deconvolution_method.SFPDnmr.xml')
%   -> Suplied files necessary to build nmrdbfile which is loacated in 'apppath\database\'(default name 'nmrdb.mat')
%   - massfile: file of mass matrix  loacated in 'apppath\database\' (default name 'nmropendbcorrtable.mat')
%   - loadingfile: file of loading base  loacated in 'apppath\database\' (default name 'nmopenbase.xml')
%       
% Output folder (to store result bases and figures) : default = fullfile(root,'output')
%
% 15/10/2014 - INRA\Olivier Vitrac, Mai Nguyen - rev.
% History


%% install the current toolbox if needed
currentapplication = mfilename('fullpath');
currentpathtoolbox = rootdir(currentapplication); currenttoolbox = lastdir(currentpathtoolbox);
checktoolboxinstall(currenttoolbox,currentpathtoolbox)

%% default arguments and check syntax
default = struct('version',0.1,...
                 'ext','.SFPDnmr.xml',... common extensions for config files
                 'prefetchext','.prefetch.mat',...as above for prefetch files
                 'appconfig','application',... default filename of appconfig
                 'nmrdbfile','nmrdb.mat',... default spectral database
                 'deconvmethod','deconvolution_method',... default deconvolution method
                 'outputpath','output'); % output folder
if nargin, firstparameter = varargin{1}; else firstparameter = ''; end
if ~ischar(firstparameter), error('syntax: SFPDnmr configfile'), end
o = argcheck(varargin(2:end),default);
appconfigprefetch = fullfile(tempdir,[o.appconfig o.prefetchext]);
[apppath,appconfig,ext] = fileparts(firstparameter);
if ~isempty(appconfig) % use existing configuration file
    if isempty(apppath), apppath = pwd; end
    if isempty(ext), ext = o.ext; end
    appconfig = fullfile(apppath,[appconfig ext]);
elseif  ~isempty(apppath) % use existing folder
    appconfig = fullfile(apppath,[o.appconfig o.ext]);
else % no path, no configuration file
    % try to reuse a previous prefetch file
    if exist(appconfigprefetch,'file')
        load(appconfigprefetch,'-mat','apppath')
        load(appconfigprefetch,'-mat','appconfig')
    else % use the current folder
        apppath = pwd;
        appconfig = fullfile(apppath,'private',[o.appconfig o.ext]);
    end
end
% check the validity of apppath
if ~exist(apppath,'dir'), error('the supplied application folder ''%s'' does not exist',apppath); end
% refresh appconfigprefetch
appconfigprefetch = fullfile(tempdir,fnamexml2mat(appconfig,o.prefetchext));

%% generate a default configuration if it does not exist and save the corresponding prefetch
% Several possibilities
%       fresh install
%       prefetch appconfig file in temp dir
%       appconfig in currentpathtoolbox/private/
header = {'SafeFoodPack Design Project (ANR-10-ALIA-009)'
          'Authors: Mai Nguyen, Cedric Lyathaud, Olivier Vitrac'
          'Contact: olivier.vitrac@agroparistech.fr'
          datestr(now)
          localname};
% save in xml format: appconfig
if ~exist(appconfig,'file')
    appdata = buildappdata(apppath,o);
    appdata = appenddefaultmethod(appdata);
    xml_save(appconfig,appdata,'on',header)
end
% save in preftech file .mat
if ~exist(appconfigprefetch,'file')
    appdata = buildappdata(apppath,o);
    appdata = appenddefaultmethod(appdata);
    save(appconfigprefetch,'-struct','appdata')
end

%% check whether the prefetch file needs to be refreshed
if olderthan(appconfigprefetch,appconfig)
    % load appconfig
    appdata = xml_load(appconfig);
    % save prefetch
    save(appconfigprefetch,'-struct','appdata')
end

%% load appconfig prefetch: appdata
appdata = load(appconfigprefetch);
% check whether apppath, appprivate are valid
if ~exist(appdata.apppath,'dir'), error('the supplied application folder ''%s'' does not exist',appdata.apppath); end
if ~exist(appdata.privatepath,'dir'), error('the supplied private folder ''%s'' does not exist',appdata.privatepath); end
% check and create if needed method and output folder
if ~exist(rootdir(appdata.methodfile),'dir'), mkdir(rootdir(appdata.methodfile)), end
if ~exist(appdata.outputpath,'dir'), mkdir(appdata.outputpath), end

%% create and save NMRDB for the first launch
if ~exist(appdata.nmrdbfile,'file')
    dispf('Create NMRDB file for the first application launching: %s',appdata.nmrdbfile)
    buildnmrdb(appdata.loadingfile,appdata.massfile,'outputpath',rootdir(appdata.nmrdbfile),'dbfilename',lastdir(appdata.nmrdbfile))
    dispf('Reference NMRDB file:'), fileinfo(appdata.nmrdbfile)
end

%% save defaultmethod in user xml file for the first launhching
if ~exist(appdata.methodfile,'file')
    dispf('Create deformulation method file for the first application launching: %s',appdata.methodfile)
    method = methoddlg2method(appdata.method);
    method.nfo = appdata.nfo;
    xml_save(appdata.methodfile,method,'on',header)
end

%% call interface NMR_interface with a constructor syntax
%   send appdata to be saved as 'userdata' in the interface
NMR_interface(appdata);
  

%% PRIVATE functions
% conversion xml filename in filename.mat
function fnamemat = fnamexml2mat(fnamexml,ext)
[~,fname] = fileparts(fnamexml); 
fnamemat = [fname ext];   

% returns true when file fa is older than file fb
function flag = olderthan(fa,fb)
dfa = dir(fa);
dfb = dir(fb);
flag = dfa.datenum < dfb.datenum;

% build default application configuration: appdata
function appdata = buildappdata(apppath,opt)
header = {'SafeFoodPack Design Project (ANR-10-ALIA-009)'
          'Authors: Mai Nguyen, Cedric Lyathaud, Olivier Vitrac'
          'Contact: olivier.vitrac@agroparistech.fr'};
private = struct('xmlheader', {header},'version',opt.version);
defaultappdata = @(apppath)  ...
      struct( ... general informations
             'version',opt.version,...
             'creationdate',datestr(now),...
             'machine',localname(),...
             'user',username(),...
             'apppath',apppath,...
             'privatepath',fullfile(apppath,'private'),...
             'extension',opt.ext,...
             'prefetchext',opt.prefetchext,...
             'loadingfile',fullfile(apppath,'database','nmropenbase.xml'),... loading file
             'massfile',fullfile(apppath,'database','nmropendbcorrtable.mat'),... mass matrix file
             'nmrdbfile',fullfile(apppath,'database',opt.nmrdbfile),...nmrdb file
             'appconfig',fullfile(apppath,'private',[opt.appconfig opt.ext]),... name of app configuration file
             'methodfile',fullfile(apppath,'method',[opt.deconvmethod opt.ext]),... name of methodology configuration file
             'outputpath',fullfile(apppath,opt.outputpath),... folder to store results
             'mixturespectrum','',... spectrum file name (fullpath) of last used mixture to be deformulated
             'projectname','',... name of project or fileparts of mixture name
             'nmolchoose',[10 1 52],... number of substances for quantitative deconvolution [likely min max]
             'rhothreshold',[0.75 .1 1],... threshold on correlation coeff to extract number of susbtances for quantitative deconvolution [likely min max]
             'identificationext','.identification.mat',... extension of database of identification
             'deconvolutionext','.deconvolution.mat',... extension of database of deconvolution
             'windowposition',[131.8000   30.5385  165.8000   37.0000],... window position of interface
             'private',private,...
             'nfo',struct('date',datestr(now),'version',opt.version,'host',localname(),'user',username()));
appdata = defaultappdata(apppath);

% constructor of default deformulation method
function appdata = appenddefaultmethod(appdata)
Identification = struct('worse',struct('value',10,'name','"*" factor applied on correlation coefficient','bounds',[1 50],'flag',true),...
                  'threshold',struct('value',0.6,'name','threshold on correlation coefficient for appliyng worse factor','bounds',[0.1 1],'flag',true),...
                  'idsubstance',struct('value','commonname','name','type of name of substances','flag',true));
Report1 = struct('figname',struct('value','report','name','figure name'),...
                        'paperposition',struct('value',[0.1920   -0.3226   20.6000   30.0000],'name','paper position'),...
                        'paperorientation',struct('value','portrait','name','paper orientation'),...
                        'marginlayout',struct('value',[.1 .1 .21],'name','margin of layout','bounds',[0 1]),...
                        'spacerheader',struct('value',0.01,'name','header spacer for plot area','bounds',[0 1]),...
                        'spacerfooter',struct('value',0.05,'name','footer spacer for plot area','bounds',[0 1]),...
                        'margincolorbar',struct('value',0.1,'name','margin of colorbar','bounds',[0 1]),...
                        'margingraph',struct('value',0.1,'name','margin of scenario tree graph','bounds',[0 1]),...
                        'headertext',struct('value','','name','header text'),...
                        'footertext',struct('value','','name','footer text'),...
                        'fontsizeheadertext',struct('value',12,'name','fontsize of header text','bounds',[1 100]),...
                        'fontsizefootertext',struct('value',8,'name','fontsize of footer text','bounds',[0 100]),...
                        'resolution',struct('value','300','name','print resolution','bounds',[0 1000]));
Report2 = Report1;
Report2.paperposition = struct('value',[3.0301    1.2380   23.6173   18.5080],'name','paper position');
Report2.paperorientation = struct('value','landscape','name','paper orientation');
Report2.spacerheater = struct('value',0.06,'name','heater spacer for plot area','bounds',[0 1]);
Report2.spacerfooter = struct('value',0.06,'name','footer spacer for plot area','bounds',[0 1]);
defaultmethod = struct( 'nclasstheorical',struct('value',18,'name','Number of theoretical classes','bounds',[1 52],'flag',true),...
                        'rhothresholdclass',struct('value',0.75,'name','threshold on correlation coefficient for ranking substances by class','flag',true,'bounds',[0 1]),...
                        'nmoltree',struct('value',5,'name','number of substances to be plotted in scenario tree','bounds',[1 52]),...
                        'nplotfit',struct('value',4,'name','number of spline approximation to be plotted','bounds',[1 52]),...
                        'concmax',struct('value',10000,'name','Maximal value of concentration for color scale','bounds',[1 1e6]),...
                        'colormapgraph',struct('value','cmaplog(50)','name','colormap graph'));
defaultmethod.Identification = struct('value',Identification);
defaultmethod.Report1 = struct('value',Report1); 
defaultmethod.Report2 = struct('value',Report2); 
appdata.method = defaultmethod;

% Create method structure to be saved in xml from a defaultmethod structure
function method = methoddlg2method(defaultmethoddlg)
for f = fieldnames(defaultmethoddlg)'
    if isstruct(defaultmethoddlg.(f{1}).value)
        submethod = methoddlg2method(defaultmethoddlg.(f{1}).value);
        method.(f{1}) = submethod;
    else
        method.(f{1}) = defaultmethoddlg.(f{1}).value; 
    end
end