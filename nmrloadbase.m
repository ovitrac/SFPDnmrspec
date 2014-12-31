function [dbpur,dbxpur,dbmask,mask,dbmultiplet,dbmix,dbxmix] = nmrloadbase(varargin)
%NMRLOADBASE loads necessary nmr databases (mask, spectra NMR base...)  
% Syntax: [dbxpur,dbmask] = nmrloadbase(varargin)
% List of input pair parameters/values:
%       dbpur: NMR spectra base of pur substances  (default = [], see nmrloadascii)
%       dbxpur: concatened NMR spectra base of pur substances (default = [], see nmrloadascii)
%       dbmask: mask database of pur substances (default = [])
%       dbmix: NMR spectra base of mixtures  (default = [], see nmrloadascii)
%       dbxmix: concatened NMR spectra base of mixtures (default = [], see nmrloadascii)
%       nmrpurlocal: local directory (inside of nmrspec toolbox) where ODS and TXT files of NMR spectra of pur substances are located (default = '')
%       nmrpurfile: ods file containing substances information (default = '', see nmrloadascii)
%       nmrpursheetname: worksheet of ODS file to be read (default = '', see nmrloadascii)
%       nmrpurprimarykey: name of the column listing codes to be used as fieldnames in Matlab (default = '', see nmrloadascii)
%       nmrmixlocal: local directory (inside of nmrspec toolbox) where ODS and TXT files of NMR spectra of mixtures are located (default = '')
%       nmrmixfile: ods file containing mixtures information (default = '', see nmrloadascii)
%       nmrmixsheetname: worksheet of ODS file to be read (default = '', see nmrloadascii)
%       nmrmixprimarykey: name of the column listing codes to be used as fieldnames in Matlab (default = '', see nmrloadascii)
%       masklocal: local directory (inside of nmrspec toolbox) where ODS file of mask is located (default = '')
%       maskfile: ods file containing mask data (default = '')
%
% EXAMPLE
% [dbpur,dbxpur,dbmask,mask,dbmultiplet,dbmix,dbxmix] = nmrloadbase; % load default nmr databases
%
% see also : NMRLOADASCII, NMRLOADMASK
% RMNSPEC v 0.1 - 18/03/13 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev.01/11/13 
% history
% 01/11/13 refresh options of functions nmrloadascci and nmrloadmask
% default
default = struct(...
                'dbpur', [],...
                'dbxpur',[],...
                'dbmask',[],...
                'mask',[],...
                'dbmultiplet',[],...
                'dbmix',[],...
                'dbxmix',[],...'dbfit',[],...
                'nmrpurlocal','',...
                'nmrpurfile','',...
                'nmrpursheetname','',...
                'nmrpursheetnamenfo','',...
                'nmrpurprimarykey','',...
                'nmrmixlocal','',...
                'nmrmixfile','',...
                'nmrmixsheetname','',...
                'nmrmixsheetnamenfo','',...
                'nmrmixprimarykey','',...
                'masklocal','',...
                'maskfile',''...
                );
% keyword = 'calib';
% argcheck
o = argcheck(varargin,default); %,keyword);

% main
if isempty(o.dbpur) && evalin('base','exist(''dbpur'',''var'')'), o.dbpur = evalin('base','dbpur'); end
if isempty(o.dbxpur) && evalin('base','exist(''dbxpur'',''var'')'), o.dbxpur = evalin('base','dbxpur'); end
if isempty(o.dbmask) && evalin('base','exist(''dbmask'',''var'')'), o.dbmask = evalin('base','dbmask'); end
if isempty(o.dbmix) && evalin('base','exist(''dbmix'',''var'')'), o.dbmix = evalin('base','dbmix'); end
if isempty(o.dbxmix) && evalin('base','exist(''dbxmix'',''var'')'), o.dbxmix = evalin('base','dbxmix'); end

% without user's specific definition
if nargout<3
    if isempty(o.nmrpurfile) && isempty(o.nmrpurlocal)
        [dbpur,dbxpur] = nmrloadascii('path',fullfile(find_path_toolbox('rmnspec'),'data_pur'),'pur'); % load default dbpur as defined in nmrloadascii
%         if o.calib
%             [dbpur,dbxpur] = nmrloadascii('path',fullfile(find_path_toolbox('rmnspec'),'data_pur'),'pur','calib');
%         else 
%         end
    end
end
if nargout<6
    if isempty(o.nmrpurfile) && isempty(o.nmrpurlocal)
        [dbpur,dbxpur] = nmrloadascii('path',fullfile(find_path_toolbox('rmnspec'),'data_pur'),'pur');      
    end
    if isempty(o.maskfile) && isempty(o.masklocal)
        [dbmask,mask,dbmultiplet] = nmrloadmask('path',fullfile(find_path_toolbox('rmnspec'),'data_pur')); % load default dbmask as defined in nmrloadmask
    end
end
if nargout>5
    if isempty(o.nmrpurfile) && isempty(o.nmrpurlocal)
        [dbpur,dbxpur] = nmrloadascii('path',fullfile(find_path_toolbox('rmnspec'),'data_pur'),'pur');
    end
    if isempty(o.maskfile) && isempty(o.masklocal)
        [dbmask,mask,dbmultiplet] = nmrloadmask('path',fullfile(find_path_toolbox('rmnspec'),'data_pur')); % load default dbmask as defined in nmrloadmask
    end
    if isempty(o.nmrmixfile) && isempty(o.nmrmixlocal)
        [dbmix,dbxmix] = nmrloadascii('path',fullfile(find_path_toolbox('rmnspec'),'data_mixture'),'odsfile','extract.ods','sheetname','extract','sheetnamenfo','nfoextract','primarykey','reference'); % load default dbmix
    end
end

% user defines ods files or their local
if ~ischar(o.nmrpurfile) || ~ischar(o.maskfile) || ~ischar(o.nmrmixfile), error('ODS file names must be a string'), end
if nargout<3    
    if ~isempty(o.nmrpurfile) || ~isempty(o.nmrpurlocal)
        [dbpur,dbxpur] = nmrloadascii('path',fullfile(find_path_toolbox('rmnspec'),o.nmrpurlocal),'odsfile',o.nmrpurfile,'sheetname',o.nmrpursheetname,...
                                      'sheetnamenfo',o.nmrpursheetnamenfo,'primarykey',o.nmrpurprimarykey,'pur'); 
%         if o.calib
%             [dbpur,dbxpur] = nmrloadascii('path',fullfile(find_path_toolbox('rmnspec'),o.nmrpurlocal),'odsfile',o.nmrpurfile,'sheetname',o.nmrpursheetname,...
%                                           'sheetnamenfo',o.nmrpursheetnamenfo,'primarykey',o.nmrpurprimarykey,'pur','calib','dbfit',o.dbfit);
%         else  
%         end
    end
end
if nargout<6
    if ~isempty(o.nmrpurfile) || ~isempty(o.nmrpurlocal)
        [dbpur,dbxpur] = nmrloadascii('path',fullfile(find_path_toolbox('rmnspec'),o.nmrpurlocal),'odsfile',o.nmrpurfile,'sheetname',o.nmrpursheetname,...
                                           'sheetnamenfo',o.nmrpursheetnamenfo,'primarykey',o.nmrpurprimarykey,'pur'); 
    end
    if ~isempty(o.maskfile) || ~isempty(o.masklocal)
        [dbmask,mask,dbmultiplet] = nmrloadmask('path',fullfile(find_path_toolbox('rmnspec'),o.masklocal),'odsfile',o.maskfile);
    end
end
if nargout>5
    if ~isempty(o.nmrpurfile) || ~isempty(o.nmrpurlocal)
        [dbpur,dbxpur] = nmrloadascii('path',fullfile(find_path_toolbox('rmnspec'),o.nmrpurlocal),'odsfile',o.nmrpurfile,'sheetname',o.nmrpursheetname,...
                                           'sheetnamenfo',o.nmrpursheetnamenfo,'primarykey',o.nmrpurprimarykey,'pur'); 
    end
    if ~isempty(o.maskfile) || ~isempty(o.masklocal)
        [dbmask,mask,dbmultiplet] = nmrloadmask('path',fullfile(find_path_toolbox('rmnspec'),o.masklocal),'odsfile',o.maskfile);
    end
    if  ~isempty(o.nmrmixfile) || ~isempty(o.nmrmixlocal)
        [dbmix,dbxmix] = nmrloadascii('path',fullfile(find_path_toolbox('rmnspec'),o.nmrmixlocal),'odsfile',o.nmrmixfile,'sheetname',o.nmrmixsheetname,'sheetnamenfo',o.nmrmixsheetnamenfo,'primarykey',o.nmrmixprimarykey);  
    end
end
