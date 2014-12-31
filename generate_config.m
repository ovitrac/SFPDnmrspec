% script generate default configuration files for deformulation ued in opensource project
% FILE 1: Configuration file for application SFPDnmr
% FILE 2: Configuration file for interface and calculation parameters
% all parameters for deformulate by 2 functions : deformulateNMR and plotdeformulateNMR
% config structure will be used by StrucDlgPrepare and StructDlg (for display interface)
% 07/11/2014

%% building config structure for calculation
% structure
close all, clear all
outputpath = find_path_toolbox('SFPDnmrspec');
appconfigname = 'application';
methodfile = 'deconvolution_method';
extension = '.SFPDnmr.xml';
prefetchext = '.prefetch.mat';
nmrdbfile = 'nmrdb.mat'; 
Identification = struct('worse',struct('value',10,'name','"*" factor applied on correlation coefficient','bounds',[1 50],'flag',true),...
                  'threshold',struct('value',0.6,'name','threshold on correlation coefficient for appliyng worse factor','bounds',[0.1 1],'flag',true),...
                  'idsubstance',struct('value','commonname','name','type of name of substances','flag',true));
Report1 = struct('figname',struct('value','report','name','figure name'),...
                        'paperposition',struct('value',[0.5920    3.0424   19.8000   23.5927],'name','paper position'),...
                        'paperorientation',struct('value','portrait','name','paper orientation'),...
                        'marginlayout',struct('value',[.2 .1 .21],'name','margin of layout','bounds',[0 1]),...
                        'spacerheader',struct('value',0.01,'name','header spacer for plot area','bounds',[0 1]),...
                        'spacerfooter',struct('value',0.06,'name','footer spacer for plot area','bounds',[0 1]),...
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
Report2.spacerheater = struct('value',0.08,'name','heater spacer for plot area','bounds',[0 1]);
Report2.spacerfooter = struct('value',0.08,'name','footer spacer for plot area','bounds',[0 1]);
defaultmethod = struct( 'nclasstheorical',struct('value',18,'name','Number of theoretical classes','bounds',[1 52],'flag',true),...
                        'rhothresholdclass',struct('value',0.75,'name','threshold on correlation coefficient for ranking substances by class','flag',true,'bounds',[0 1]),...
                        'nmoltree',struct('value',5,'name','number of substances to be plotted in scenario tree','bounds',[1 52]),...
                        'nplotfit',struct('value',4,'name','number of spline approximation to be plotted','bounds',[1 52]),...
                        'concmax',struct('value',10000,'name','Maximal value of concentration for color scale','bounds',[1 1e6]),...
                        'colormapgraph',struct('value','cmaplog(50)','name','colormap graph'));
defaultmethod.Identification = struct('value',Identification);
defaultmethod.Report1 = struct('value',Report1); 
defaultmethod.Report2 = struct('value',Report2); 
s = StructDlgPrepare(defaultmethod);
method = StructDlg(s);

%% buiding config for SFPDnmr application
version = 0.1;
header = {'SafeFoodPack Design Project (ANR-10-ALIA-009)'
          'Authors: Mai Nguyen, Cedric Lyathaud, Olivier Vitrac'
          'Contact: olivier.vitrac@agroparistech.fr'};
private = struct('xmlheader', {header},'version',version);
appdata = struct('version',private.version,...
                 'creationdate',datestr(now),...
                 'machine',localname(),...
                 'user',username(),...
                 'apppath',outputpath,...
                 'privatepath',fullfile(outputpath,'private'),...
                 'extension',extension,...
                 'prefetchext',prefetchext,...
                 'loadingfile',fullfile(outputpath,'database','nmropenbase.xml'),... loading file
                 'massfile',fullfile(outputpath,'database','nmropendbcorrtable.mat'),... mass matrix fil
                 'nmrdbfile',fullfile(outputpath,'database',nmrdbfile),...
                 'appconfig',fullfile(outputpath,'private',[appconfigname extension]),... name of app configuration file
                 'methodfile',fullfile(outputpath,'method',[methodfile extension]),... name of methodology configuration file
                 'mixturespectrum','',... spectrum file name (fullpath) of last used mixture to be deformulated
                 'projectname','',... name of project or fileparts of mixture name
                 'outputpath',fullfile(outputpath,'output'),... folder to store results
                 'nmolchoose',[10 1 52],... number of substances for quantitative deconvolution
                 'rhothreshold',[0.75 .1 1],... threshold on correlation coeff to extract number of susbtances for quantitative deconvolution
                 'identificationext','.identification.mat',... extension of database of identification
                 'deconvolutionext','.deconvolution.mat',... extension of database of deconvolution
                 'windowposition',[131.8000   30.5385  165.8000   37.0000],... position of interface
                 'private',private,...
                 'method',defaultmethod...
                 );
             % windowposition 
%% save
% save app config
header = {'SafeFoodPack Design Project (ANR-10-ALIA-009)'
          'Authors: Mai Nguyen, Cedric Lyathaud, Olivier Vitrac'
          'Contact: olivier.vitrac@agroparistech.fr'
          datestr(now)
          localname};
xml_save(appdata.appconfig,appdata,'on',header)
[~,fname] = fileparts(appdata.appconfig);
save(fullfile(tempdir,[fname prefetchext]),'-struct','appdata')

% save method param
xml_save(appdata.methodfile,method,'on',header)
[~,fname] = fileparts(appdata.methodfile);
save(fullfile(appdata.apppath,'method',[fname prefetchext]),'-struct','method')
