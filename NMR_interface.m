function varargout = NMR_interface(varargin)
% NMR_INTERFACE MATLAB code for NMR_interface.fig
%      NMR_INTERFACE, by itself, creates a new NMR_INTERFACE or raises the existing
%      singleton*.
%
%      H = NMR_INTERFACE returns the handle to a new NMR_INTERFACE or the handle to
%      the existing singleton*.
%
%      NMR_INTERFACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NMR_INTERFACE.M with the given input arguments.
%
%      NMR_INTERFACE('Property','Value',...) creates a new NMR_INTERFACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NMR_interface_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NMR_interface_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NMR_interface

% Last Modified by GUIDE v2.5 26-Nov-2014 22:56:25


%% ===============================
%       APP CONFIGURATION
if ~nargin, SFPDnmr, return; end

% ===============================

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NMR_interface_OpeningFcn, ...
                   'gui_OutputFcn',  @NMR_interface_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before NMR_interface is made visible.
function NMR_interface_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NMR_interface (see VARARGIN)

% Choose default command line output for NMR_interface
handles.output = hObject;

% update old data and fill interface
if nargin>3
    appdata = varargin{1};
    if ~isempty(appdata) && isstruct(appdata) && isfield(appdata,'appconfig')
        set(handles.figure1,'UserData',appdata)
        set(handles.output_text,'String',appdata.outputpath,'UserData',appdata.outputpath) 
        % built structure for UserData of Callback user_config_text
        userconfig = struct('method',appdata.method,'methodfile',appdata.methodfile); % method in format Dlg stored in appdata 
        set(handles.user_config_text,'String',appdata.methodfile,'UserData',userconfig)
        set(handles.dbfile_text,'String',appdata.nmrdbfile,'UserData',appdata.nmrdbfile)
        set(handles.mix_text,'String',appdata.mixturespectrum,'UserData',appdata.mixturespectrum)
        set(handles.projname_text,'String',appdata.projectname,'UserData',appdata.projectname)
        set(handles.nmolchoosse_text,'String',appdata.nmolchoose(1),'UserData',appdata.nmolchoose(1))
        set(handles.rhothreshold_text,'String',appdata.rhothreshold(1),'UserData',appdata.rhothreshold(1))
    end
    subplot(handles.logo); imshow(imread(fullfile(appdata.privatepath,'logo.png')))
    subplot(handles.logo0); imshow(imread(fullfile(appdata.privatepath,'logo0.png')))
    set(handles.figure1,'Name',sprintf('SFPDnmrspec - v %0.2g',appdata.version))
    set(handles.appname_text,'String',sprintf('SFPDnmrspec - v %0.2g',appdata.version))
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NMR_interface wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NMR_interface_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% delete(hObject);


function projname_text_Callback(hObject, eventdata, handles)
% hObject    handle to projname_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of projname_text as text
%        str2double(get(hObject,'String')) returns contents of projname_text as a double
parser = '[^0-9a-zA-Z\+\-\_\s]';
pname = regexprep(get(hObject,'String'),parser,'');
set(hObject,'String',pname,'UserData',pname)

% --- Executes during object creation, after setting all properties.
function projname_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to projname_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function mix_text_Callback(hObject, eventdata, handles)
% hObject    handle to mix_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mix_text as text
%        str2double(get(hObject,'String')) returns contents of mix_text as a double
% dispf('mix_text_Callback: %s (%s)',get(hObject,'String'),get(hObject,'UserData'))
fname = substructarray(get(handles.figure1,'Userdata'),'mixturespectrum');
f = validfile(get(hObject,'String'), ... current text
              fname{1} ... default
              );
set(hObject,'String',f,'UserData',f)
% implement the dependence with projectname (the filename is used as project name)
p = makeprojectname(get(hObject,'String'));
set(handles.projname_text,'String',p,'UserData',p)

% --- Executes during object creation, after setting all properties.
function mix_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mix_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in mix_choose.
function mix_choose_Callback(hObject, eventdata, handles)
% hObject    handle to mix_choose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[mixfile,mixfilepath] = uigetfile('*.txt','Select the mixture NMR spectrum'); %get(handles.mix_text,'String')
set(handles.mix_text,'String',fullfile(mixfilepath,mixfile),'UserData',fullfile(mixfilepath,mixfile))
% implement the dependence with projectname (the filename is used as project name)
p = makeprojectname(mixfile);
set(handles.projname_text,'String',p,'UserData',p)


function output_text_Callback(hObject, eventdata, handles)
% hObject    handle to output_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output_text as text
%        str2double(get(hObject,'String')) returns contents of output_text as a double
% dispf('output_text_Callback: %s (%s)',get(hObject,'String'),get(hObject,'UserData'))
pname = substructarray(get(handles.figure1,'Userdata'),'outputpath');
p = validpath(get(hObject,'String'), ... current text
              pname{1} ... default
              );
set(hObject,'String',p,'UserData',p)

% --- Executes during object creation, after setting all properties.
function output_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in output_choose.
function output_choose_Callback(hObject, eventdata, handles)
% hObject    handle to output_choose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
outputpath = uigetdir('','Please choose ouputpath');
set(handles.output_text,'String',outputpath,'UserData',outputpath)


function user_config_text_Callback(hObject, eventdata, handles)
% hObject    handle to user_config_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of user_config_text as text
%        str2double(get(hObject,'String')) returns contents of user_config_text as a double
% dispf('user_config_text_Callback: %s (%d)',get(hObject,'String'),get(hObject,'UserData'))
% control valid name
fname = substructarray(get(handles.figure1,'Userdata'),'methodfile');
f = validfile(get(hObject,'String'), ... current text
              fname{1} ... default
              );
requiredfields = {'Identification' 'Report1' 'Report2' 'nclasstheorical' 'rhothresholdclass' 'nmoltree' 'nplotfit' 'concmax' 'colormapgraph'};          
method = safeload(f,'xml',requiredfields); % method in xml format (short)
% update userconfig 
userconfig = get(hObject,'UserData'); 
appdata = get(handles.figure1,'UserData');
userconfig.methodfile = f; userconfig.method = method2methoddlg(appdata.method,method);
set(hObject,'String',f,'UserData',userconfig)

% --- Executes during object creation, after setting all properties.
function user_config_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to user_config_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in user_config_choose.
function user_config_choose_Callback(hObject, eventdata, handles)
% hObject    handle to user_config_choose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[userconfigfile,userconfigpath] = uigetfile('*.xml','Select the deformulation method file',get(handles.user_config_text,'String'));
userconfig = get(handles.user_config_text,'UserData'); userconfig.methodfile = fullfile(userconfigpath,userconfigfile);
set(handles.user_config_text,'String',fullfile(userconfigpath,userconfigfile))
% control user_config_text
user_config_text_Callback(handles.config_user_text,eventdata,handles)

% --- Executes on button press in config_set.
function config_set_Callback(hObject, eventdata, handles)
% hObject    handle to config_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% refresh appdata
appdata = refreshappdata(get(handles.figure1,'UserData'),handles); %appdata = get(handles.figure1,'UserData');
requiredfields = {'Identification' 'Report1' 'Report2' 'nclasstheorical' 'rhothresholdclass' 'nmoltree' 'nplotfit' 'concmax' 'colormapgraph'};
method = safeload(get(handles.user_config_text,'string'),'xml',requiredfields);
% tranform into dlg structure
methoddlg = method2methoddlg(appdata.method,method);
% add methodfile so that user can specifie a nex name for saving
methoddlg.methodfile = struct(struct('value',get(handles.user_config_text,'String'),'name','method filename'));
r = StructDlg(StructDlgPrepare(methoddlg),'Deformulation parameters');
if ~isempty(r)    
    safesave(rmfield(r,'methodfile'),r.methodfile,'.xml','',[],handles.figure1,'Method file already exists. Do you want to overwrite it?','Deformulation method')
    set(handles.user_config_text,'string',r.methodfile)    
    % control user_config_text
    user_config_text_Callback(handles.user_config_text,eventdata,handles)
end


function dbfile_text_Callback(hObject, eventdata, handles)
% hObject    handle to dbfile_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dbfile_text as text
%        str2double(get(hObject,'String')) returns contents of dbfile_text as a double
% dispf('dbfile_text_Callback: %s (%d)',get(hObject,'String'),get(hObject,'UserData'))
fname = substructarray(get(handles.figure1,'Userdata'),'nmrdbfile');
f = validfile(get(hObject,'String'), ... current text
              fname{1} ... default
             );
set(hObject,'String',f,'UserData',f)

% --- Executes during object creation, after setting all properties.
function dbfile_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dbfile_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dbfile_choose.
function dbfile_choose_Callback(hObject, eventdata, handles)
% hObject    handle to dbfile_choose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[dbfile,dbfilepath] = uigetfile('*.mat','Select the NMR database',get(handles.dbfile_text,'String'));
set(handles.dbfile_text,'String',fullfile(dbfilepath,dbfile),'UserDaa',fullfile(dbfilepath,dbfile))
% control handles.dbfile_text 
dbfile_text_Callback(handles.dbfile_text,eventdata,handles)


% --- Executes on button press in dbfile_edit.
function dbfile_edit_Callback(hObject, eventdata, handles)
% hObject    handle to dbfile_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% refresh appdata
appdata = refreshappdata(get(handles.figure1,'UserData'),handles); %appdata = get(handles.figure1,'UserData');
p = struct('nmrdbfile',struct('value',appdata.nmrdbfile,'name','nmrdb filename'));
p.loadingfile = struct('value',appdata.loadingfile,'name','Choose a loading file','bounds',{{'*.xml','XML files'}},'flag','getsinglefile');
p.massfile = struct('value',appdata.massfile,'name','Choose a mass file','bounds',{{'*.mat'}},'flag','getsinglefile');
[s,figposition] = StructDlgPrepare(p);
r = StructDlg(s,'Choose files for create nmrdb',[],figposition);
if ~isempty(r) % button "ok": regenerate db
    if exist(appdata.nmrdbfile,'file')
        dispf('NMRDB already exists. File details:'),fileinfo(appdata.nmrdbfile) 
        answer = askoverwrite(appdata.nmrdbfile,'NMRDB file already exists. Do you want to overwrite it?','Spectral library',false); 
        if answer
            movefile(appdata.nmrdbfile,sprintf('%s~',appdata.nmrdbfile))
            dispf('Backup file: %s'), sprintf('%s~',appdata.nmrdbfile)
            buildnmrdb(r.loadingfile,r.massfile,'outputpath',rootdir(r.nmrdbfile),'dbfilename',lastdir(r.nmrdbfile))
        end 
    else buildnmrdb(r.loadingfile,r.massfile,'outputpath',rootdir(r.nmrdbfile),'dbfilename',lastdir(r.nmrdbfile))        
    end
    % control handles.dbfile_text 
    dbfile_text_Callback(handles.dbfile_text,eventdata,handles)
end


% --- Executes on button press in pair_iden_launch.
function pair_iden_launch_Callback(hObject, eventdata, handles)
% hObject    handle to pair_iden_launch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% refresh appdata
appdata = refreshappdata(get(handles.figure1,'UserData'),handles);
method = methoddlg2method(appdata.method);
% check prefetch .mat data of pair correlation coefficient: to be tored in apppath\output\projectname\
dbname = fullfile(appdata.outputpath,appdata.projectname,sprintf('%s%s',appdata.projectname,appdata.identificationext));
if exist(dbname,'file')  
    dispf('The identification file for mixture ''%s'' already exists. File details:',appdata.projectname),fileinfo(dbname) 
    answer = askoverwrite(dbname,'The identification file exists. Do you want to overwrite it?',sprintf('Project %s',appdata.projectname),false);   
    if answer 
       movefile(dbname,sprintf('%s~',dbname))
       dispf('Backup of older file: %s'), sprintf('%s~',dbname)
       dispf('Calculation in porgress...')
       R1 = nmrpairidentification(appdata.nmrdbfile,...
                            appdata.mixturespectrum,... 
                           'worse',method.Identification.worse,...
                           'idsubstance',method.Identification.idsubstance,... 
                           'threshold',method.Identification.threshold,... 
                           'outputpath',appdata.outputpath,...
                           'mixturename',appdata.projectname,...
                           'dbnameext',appdata.identificationext) ;
       R1.method = method;
       save(dbname,'-struct','R1')
       dispf('End of calculation')
    end
else
    dispf('Calculation in porgress...')
    R1 = nmrpairidentification(appdata.nmrdbfile,...
                           appdata.mixturespectrum,... 
                           'worse',method.Identification.worse,...
                           'idsubstance',method.Identification.idsubstance,... 
                           'threshold',method.Identification.threshold,... 
                           'outputpath',appdata.outputpath,...
                           'mixturename',appdata.projectname,...
                           'dbnameext',appdata.identificationext);
    R1.method = method;
    save(dbname,'-struct','R1')
    dispf('End of calculation.')
end


% --- Executes on button press in pair_iden_report.
function pair_iden_report_Callback(hObject, eventdata, handles)
% hObject    handle to pair_iden_report (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% refresh appdata
appdata = refreshappdata(get(handles.figure1,'UserData'),handles);
method = methoddlg2method(appdata.method);
% prepare header (project name, inputfile) and footer text (nmrdb, methodfile, output,)of reports 
headertext = sprintf('\\fontsize{14}Sample name: %s\n\\fontsize{8}Input spectrum: %s',appdata.projectname,appdata.mixturespectrum);
footertext = {sprintf('\\bfNMR database: \\rm%s',appdata.nmrdbfile)
              sprintf('\\bfMethod file: \\rm%s',appdata.methodfile)
              sprintf('\\bfIdentification data: \\rm%s%s',appdata.projectname,appdata.identificationext)
              sprintf('\\bfOutputpath \\rm%s',fullfile(appdata.outputpath,appdata.projectname))
             };
% check and load data of identification
dbname = fullfile(appdata.outputpath,appdata.projectname,sprintf('%s%s',appdata.projectname,appdata.identificationext));
if ~exist(dbname,'file')
    error('The identification file for mixture %s does not exists. Launch "PAIR IDENTIFICATION"',get(handles.projname_text,'string'))
else
    R1 = safeload(dbname,'mat',{'corrmeanweightsort'  'molsortbyrho' 'mixturename'}); 
end
uistack(handles.figure1, 'bottom')
out = plotdeformulateNMR(R1,'plotrho','margincolorbar',0.2,...
                 'paperposition',method.Report1.paperposition,...
                 'paperorientation',method.Report1.paperorientation,...
                 'marginlayout',method.Report1.marginlayout,...
                 'spacerheater',method.Report1.spacerheader,...
                 'spacerfooter',method.Report1.spacerfooter,...
                 'margincolorbar',method.Report1.margincolorbar,...
                 'margingraph',method.Report1.margingraph,...
                 'headertext',headertext,...
                 'footertext',footertext,...
                 'fontsizefootertext',method.Report1.fontsizefootertext,...
                 'fontsizeheadertext',method.Report1.fontsizeheadertext,...
                 'resolution',method.Report1.resolution,...
                 'outputpath',appdata.outputpath);
figure(out.hfig)
printfigcallback = @(src,evn) printfig(out.hfig,out.figname,out.printpath);
uicontrol('Style', 'pushbutton', 'String', 'Print','Position', [20 20 50 20],'Callback', printfigcallback);
set(out.hfig,'position',[257 9 1081 1107],'CloseRequestFcn',printfigcallback)
% disp('*')

function printfig(hfig,filename,path)
% disp(filename)
print_pdf(400,filename,path,'nocheck')
delete(hfig)

function nmolchoosse_text_Callback(hObject, eventdata, handles)
% hObject    handle to nmolchoosse_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nmolchoosse_text as text
%        str2double(get(hObject,'String')) returns contents of nmolchoosse_text as a double
% dispf('nmolchoosse_text_Callback: %s (%d)',get(hObject,'String'),get(hObject,'UserData'))
[num,fixed] = str2validnum(get(hObject,'String'), ... current text
                    substructarray(get(handles.figure1,'Userdata'),'nmolchoose'), ... [default min max]
                    get(hObject,'UserData') ... last valid entry
                    );
set(hObject,'UserData',num)
if fixed, set(hObject,'String',sprintf('%d',num)), end

% --- Executes during object creation, after setting all properties.
function nmolchoosse_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nmolchoosse_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function rhothreshold_text_Callback(hObject, eventdata, handles)
% hObject    handle to rhothreshold_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rhothreshold_text as text
%        str2double(get(hObject,'String')) returns contents of rhothreshold_text as a double
% dispf('rhothreshold_text_Callback: %s (%d)',get(hObject,'String'),get(hObject,'UserData'))
[num,fixed] = str2validnum(get(hObject,'String'), ... current text
                        substructarray(get(handles.figure1,'Userdata'),'rhothreshold'), ... [default min max]
                        get(hObject,'UserData') ... last valid entry
                          );
set(hObject,'String',sprintf('%0.2g',num),'UserData',num)

% --- Executes during object creation, after setting all properties.
function rhothreshold_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rhothreshold_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in deformulation_launch.
function deformulation_launch_Callback(hObject, eventdata, handles)
% hObject    handle to deformulation_launch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% refresh appdata
appdata = refreshappdata(get(handles.figure1,'USerData'),handles);
method = methoddlg2method(appdata.method);
% check identification data file
pairidentification = fullfile(appdata.outputpath,appdata.projectname,sprintf('%s%s',appdata.projectname,appdata.identificationext));
if ~exist(pairidentification,'file'), error('The identification file for mixture %s does not exist. Launch "PAIR IDENTIFICATION" (STEP 1)',get(handles.projname_text,'string')), end
% check prefetch .mat data of deconvolution: to be stored in apppath\output\projectname\
dbname = fullfile(appdata.outputpath,appdata.projectname,sprintf('%s%s',appdata.projectname,appdata.deconvolutionext));
if exist(dbname,'file')  
    dispf('The deconvolution file for mixture ''%s'' already exists. File details',get(handles.projname_text,'string')),fileinfo(dbname)
    answer = askoverwrite(dbname,'The deconvolution file already exists. Do you want to overwrite it?',...
                        sprintf('Project: %s',appdata.projectname),false);              
    if answer
        movefile(dbname,sprintf('%s~',dbname))
        dispf('Backup of older file: %s'), sprintf('%s~',dbname)  
        dispf('Calculation in progress ...')
        R2 = nmrdeconvolution(appdata.nmrdbfile,...
                     appdata.mixturespectrum,... 
                     pairidentification,...
                     'outputpath',appdata.outputpath,...
                     'mixturename',appdata.projectname,...
                     'nmolchoose',appdata.nmolchoose(1),...
                     'rhothreshold',appdata.rhothreshold(1),...
                     'dbnameext',appdata.deconvolutionext,...
                     'rhothresholdclass',method.rhothresholdclass,...
                     'nmoltree',method.nmoltree,...
                     'concmax',method.concmax,...
                     'nplotfit',method.nplotfit);
        R2.method = method;
        save(dbname,'-struct','R2')
        dispf('End of calculation')
    end      
else
    dispf('Calculation in porgress...')
    R2 = nmrdeconvolution(appdata.nmrdbfile,...
                     appdata.mixturespectrum,... 
                     pairidentification,...
                     'outputpath',appdata.outputpath,...
                     'mixturename',appdata.projectname,...
                     'nmolchoose',appdata.nmolchoose(1),...
                     'rhothreshold',appdata.rhothreshold(1),...
                     'dbnameext',appdata.deconvolutionext,...
                     'rhothresholdclass',method.rhothresholdclass,...
                     'nmoltree',method.nmoltree,...
                     'concmax',method.concmax,...
                     'nplotfit',method.nplotfit);
    R2.method = method;
    save(dbname,'-struct','R2')
    dispf('End of calculation')
end


% --- Executes on button press in deformulation_report.
function deformulation_report_Callback(hObject, eventdata, handles)
% hObject    handle to deformulation_report (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% refresh appdata
appdata = refreshappdata(get(handles.figure1,'USerData'),handles);
set(handles.figure1,'UserData',appdata)
method = methoddlg2method(appdata.method);
% check prefetch .mat data of deconvolution: to be stored in apppath\output\projectname\
dbname = fullfile(appdata.outputpath,appdata.projectname,sprintf('%s%s',appdata.projectname,appdata.deconvolutionext));
if ~exist(dbname,'file'), error('The deconvolution file for mixture %s does not exists. Launch "DEFORMULATION"',get(handles.projname_text,'string'))
else
    R2 = safeload(dbname,'mat',{'dbcorrmix' 'dbout' 'dataplot' 'rhotestsort' 'rang' 'mixturename' 'nplotfit'}); 
end
% plot
uistack(handles.figure1, 'bottom')
plotdeformulateNMR(R2,'plotconc',...                                       
                 'paperposition',method.Report2.paperposition,...
                 'paperorientation',method.Report2.paperorientation,...
                 'marginlayout',method.Report2.marginlayout,...
                 'spacerheater',method.Report2.spacerheader,...
                 'spacerfooter',method.Report2.spacerfooter,...
                 'margincolorbar',method.Report2.margincolorbar,...
                 'margingraph',method.Report2.margingraph,...
                 'headertext',method.Report2.headertext,...
                 'footertext',method.Report2.footertext,...
                 'fontsizefootertext',method.Report2.fontsizefootertext,...
                 'fontsizeheadertext',method.Report2.fontsizeheadertext,...
                 'resolution',method.Report2.resolution,...
                 'outputpath',appdata.outputpath);    


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
appdata = refreshappdata(get(handles.figure1,'UserData'),handles);
r = questdlg('Are you sure you want to close ?','Confirm close application','save and close','close without saving','cancel', 'save and close');
switch r
    case 'save and close'
        safesave(appdata,lastdir(appdata.appconfig),{'.xml' appdata.prefetchext},{appdata.privatepath tempdir},true,handles.figure1)       
        delete(hObject);
    case 'close without saving'
        delete(hObject);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% private UICONTROL functions
% [validnum,fixed] = str2validnum(text,[default min max],previous) returns always a valid number 
% fixed is true if the text has been replaced by previous (NaN) or default value (empty)
function [n,fixedout] = str2validnum(s,bounds,previous)
if isempty(s)
    n = bounds(1); fixed = true;
else
    fixed = false;
    n = str2double(s);
    if n<bounds(2)|| n>bounds(3), n = NaN; end
    if isnan(n), n = previous; fixed = true; end
end
if nargout>1, fixedout = fixed; end

% [f,fixed] = validfile(s,default) returns a valid file
function [f,fixedout] = validfile(s,default)
fixed = false;
if isempty(s), f = default; fixed = true; else f=s; end
if ~isempty(f) && ~exist(f,'file'), f = default; fixed = true; end
if nargout>1, fixedout = fixed; end

% [p,fixed] = validpath(s,default) returns a valid path
function [p,fixedout] = validpath(s,default)
fixed = false;
if isempty(s), p = default; fixed = true; else p=s; end
if ~isempty(p) && ~exist(p,'dir'), p = default; fixed = true; end
if nargout>1, fixedout = fixed; end

% p = makeprojectname(f) makes a valid project name from a full filename
function [p,fixedout] = makeprojectname(f)
if isempty(f), p = '';
else [~,p] = fileparts(f);
end
if nargout>1, fixedout = true; end

% Create method structure in small format (like in xml) from a method dialog structure
function method = methoddlg2method(methoddlg)
for f = fieldnames(methoddlg)'
    if isstruct(methoddlg.(f{1}).value)
        submethod = methoddlg2method(methoddlg.(f{1}).value);
        method.(f{1}) = submethod;
    else
        method.(f{1}) = methoddlg.(f{1}).value; 
    end
end

% prepare a structure to be used in StrucDlgPrepare (based on defaultmethod) with new values in method structure 
function methoddlg = method2methoddlg(defaultmethoddlg,method)
if nargin<2, method = struct([]); end
methoddlg = defaultmethoddlg;
if isempty(method), return, end
for f = intersect(fieldnames(defaultmethoddlg),fieldnames(method))'
    if isstruct(defaultmethoddlg.(f{1}).value)
        submethod = method2methoddlg(methoddlg.(f{1}).value,method.(f{1}));
        methoddlg.(f{1}).value = submethod;
    else
        methoddlg.(f{1}).value  = method.(f{1}); 
    end
end

% to refresh appdata to be stored in 'UserData' of handles.figure1
function appdata = refreshappdata(appdata,handles)
correspondingfield = struct(...
    'appfield',{'projectname' 'mixturespectrum' 'outputpath' 'methodfile' 'nmrdbfile' 'nmolchoose' 'rhothreshold' 'method' 'windowposition'},...
    'interfacefield',{'projname_text' 'mix_text' 'output_text' 'user_config_text' 'dbfile_text' 'nmolchoosse_text' 'rhothreshold_text' 'user_config_text' 'figure1'});
for i = 1:length(correspondingfield)
    if strcmp(correspondingfield(i).appfield,'windowposition')
        appdata.(correspondingfield(i).appfield) = get(handles.(correspondingfield(i).interfacefield),'position');
    elseif strcmp(correspondingfield(i).appfield,'methodfile') || strcmp(correspondingfield(i).appfield,'method')
        dtmp = substructarray(get(handles.(correspondingfield(i).interfacefield),'UserData'),(correspondingfield(i).appfield));
        appdata.(correspondingfield(i).appfield) = dtmp{1};
    elseif isnumeric(appdata.(correspondingfield(i).appfield))
        appdata.(correspondingfield(i).appfield)(1) = get(handles.(correspondingfield(i).interfacefield),'UserData');
    else
        appdata.(correspondingfield(i).appfield) = get(handles.(correspondingfield(i).interfacefield),'UserData');
    end
end
set(handles.figure1,'UserData',appdata)

% To load data from a file name (fullpath) which contains required fields
function data = safeload(filename,format,requiredfield)
data = struct([]);
% check if supplied file exists
if ~exist(filename,'file'), error('The supplied file does not exist: ''%s''',filename)
elseif strcmpi(format,'xml')
    dtmp = xml_load(filename);
elseif strcmpi(format,'mat')
    dtmp = load(filename);
end
if ~isstruct(dtmp), return, end
% check fields
ifield = 0; ok = true;
while ok && ifield < length(requiredfield)
    ifield = ifield + 1;
    ok = isfield(dtmp,requiredfield{ifield});
end
if isfield(dtmp,'nfo'), dtmp = rmfield(dtmp,'nfo'); end
if ok, data = dtmp; end

% to ask for overwitting a supplied file
function answer = askoverwrite(filename,question,title,overwrite)
if isempty(question), question = 'Do you want to ovilename,question,titleerwrite the existing file?'; end
if isempty(title), title = filename; end
if overwrite, answer = true; return, end
if ~exist(filename,'file'), answer = true; 
else
    r = questdlg(question,title,'cancel','overwrite','cancel');
    if strcmp(r,'overwrite'), answer = true; else answer = false; end
end

% to save/ overwrite a supplied file according to situation defined in data
% safesave(data,filename,ext,local,format,overwrite,private,question,title)
%   data: structure of data to be saved
%   filename: with/without path, extension
%   ext, local: empty strings, strings or cells
%   overwrite: flag (true by default forces overwrite)
%   private: private section of appdata or handle to application (handles.figure1)
%   question and title: see fuction askoverwrite
function safesave(data,filename,ext,local,overwrite,private,question,title)
if nargin < 3, ext = ''; end
if nargin < 4, local = ''; end
if nargin < 5, overwrite = true; end 
if nargin < 6, private = struct([]); end
if nargin < 7, question = ''; end
if nargin < 8, title = ''; end
if ~iscell(filename), filename = {filename}; end, filename = filename(:); 
if ~iscell(ext), ext = {ext}; end, ext = ext(:);
if ~iscell(local), local = {local}; end, local = local(:);
nfiles = max(max(length(filename),length(ext)),length(local));
filename = [filename;repmat(filename(end),nfiles-length(filename),1)];
ext = [ext;repmat(ext(end),nfiles-length(ext),1)];
local = [local;repmat(local(end),nfiles-length(local),1)];

% extract private if required
if ishandle(private), appdata = get(private,'UserData'); private = appdata.private; end
if isempty(private), vsn = []; vsnstr = ''; else vsn = private.version; vsnstr = sprintf('$SFPDnmrspec version: %0.3g$ ',vsn); end
% default header
nfo = struct('date',datestr(now),'version',vsn,'host',localname(),'user',username());
if isempty(private), header = {}; else header = private.xmlheader; end
header = [header; {sprintf('$Date: %s$ %s- $host: %s$ $user: %s',nfo.date,vsnstr,nfo.host,nfo.user)} ];
% build full filename (path,filename,extension), backup file and save
for i = 1:nfiles
    % extract current format (higher precedence for ext)
    if isempty(ext{i}), [~,~,currentext] = fileparts(filename{i}); else currentext = ext{i}; end
    format = uncell(regexp(currentext,'\.([^\.]*)$','tokens'));
    if isempty(format), format = {'mat'}; end    
    % build full filename
    fullname = fullfile(local{i},[removeext(filename{i}) ext{i}]);
    % get the consent of the user
    if askoverwrite(fullname,question,title,overwrite)
        data.nfo = nfo;
        dispf('Backup file: "%s~"',fullname)
        movefile(fullname,sprintf('%s~',fullname)) % backup file
        switch format{1} % do the save per format basis
            case 'xml', xml_save(fullname,data,'on',header)
            case 'mat', save(fullname,'-struct','data')     
        end
    end
end

% extract filename without extension 
function fnwithoutext = removeext(filename)
[p,fn] = fileparts(filename);
fnwithoutext = fullfile(p,fn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
