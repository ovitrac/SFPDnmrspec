function print_pdf(resolution,fichier,chemin,options,varargin)
%PRINT_PDF print interactively the current figure as PDF document
%   print_PDF([resolution,filename,fullpath,options])
%   print_PDF generates a standard pdf based on figure properties
%   print_PDF(resolution,fichier,options)
%   print_PDF(resolution,fichier,chemin,options)
%
%   default options = '-loose'
%   options 'nocheck' forces printing without any dialog
%   options 'append' is equivalent to nocheck but with an append feature (via a first PostScript 2.1 printing)
%   for 3D without PLOTOBJ, options = '-opengl' is recommened (see PRINT for details)
%
% Advance use of 'append' (PDF conversion from PostScript is performed for efficiency only when the last page is printed)
%       print_pdf(resolution,fichier,chemin,options,'pagenum',currentpage,'pagestart',1,'pagestop',finalpage)
%           currentpage = page index to print (when currentpage = pagestart value, the existing postscript file is removed)
%             finalpage = page to initiate PS-to-PDF printing
%   Example of use (print all figures in a single PDF)
%       listoffigs = sort(findobj('type','figure'))'; numfig = length(listoffigs);
%       for i=listoffigs
%           figure(i), print_pdf(600,'myPDFfile.pdf',pwd,'append','pagenum',i,'pagestart',1,'pagestop',numfig)
%       end

% Woodox 1.0/MS 2.1 - 13/08/07 - Olivier Vitrac - rev. 03/03/12

% revision
% 25/07/07 add options and fix chemin option
% 13/08/07 fix path ambiguity when both fichier and chemin contain a path
% 11/12/07 add fileinfo, preview and dialog
% 24/06/09 add nocheck
% 12/12/09 add evince &  for unix systems (in replacement of winopen)
% 20/07/11 fix chemin with 'nocheck'
% 21/07/11 new fix for 'nocheck'
% 22/07/11 fix extension with 'nocheck'
% 03/03/12 add append
% 12/03/12 help improvements, add example using append


% definitions
resolution_default = 600; % PDF
ext_default = '.pdf';
currentfig = gcf;
options_default = '-loose';
pmode = get(currentfig,'paperpositionmode');
psdefault = struct('pagenum',NaN,'pagestart',1,'pagestop',Inf);

% arg check
if nargin<1, resolution = []; end
if nargin<2, fichier = ''; end
if nargin<3, chemin = ''; end
if nargin<4, options = ''; end
psoptions = argcheck(varargin,psdefault);
if isempty(resolution), resolution = resolution_default; end
if isempty(fichier), [~,fichier] = fileparts(get(currentfig,'filename')); end
if isempty(fichier), [~,fichier] = fileparts(get(currentfig,'name')); end
if isempty(fichier), fichier = sprintf('Figure_%d',currentfig); end
if isempty(options), options = options_default; end
if any(chemin) 
    if chemin(1)=='-', options = chemin; chemin = ''; end
    [pathstr,name,ext] = fileparts(fichier);
    if ~isempty(pathstr) && exist(fullfile(chemin,pathstr),'dir')
        chemin = fullfile(chemin,pathstr);
    end
else
    [chemin,name,ext] = fileparts(fichier);
end
if isempty(chemin), chemin = pwd; end
if ~exist(chemin,'dir'), error('the path ''%s'' does not exist',chemin), end
if ~strcmp(ext,'.pdf'), ext = ext_default; end
fichier = fullfile(chemin,[name ext]);

% printing with no check
if strcmpi(options,'nocheck')
    print(['-r' int2str(resolution)],'-dpdf','',fichier), return
elseif strcmpi(options,'append')
    tmppsfile = regexprep(fichier,['\' ext_default '$'],'.ps','ignorecase');
    okps2pdf = true;
    if ~isnan(psoptions.pagenum) && (psoptions.pagenum>=psoptions.pagestart)
        if isinf(psoptions.pagestop) || isnan(psoptions.pagestop)
             psmsg = sprintf('page %d (first page=page %d)',psoptions.pagenum,psoptions.pagestart);
        else psmsg = sprintf('page %d/%d (first page=page %d)',psoptions.pagenum,psoptions.pagestop,psoptions.pagestart);
             okps2pdf = psoptions.pagenum>=psoptions.pagestop;
        end
    else     psmsg = '(specify ''pagenum'', ''pagestart'', ''pagestop'' for optimization)';
    end
    if exist(tmppsfile,'file') && ~isnan(psoptions.pagenum) && (psoptions.pagenum<=psoptions.pagestart), delete(tmppsfile), end
    % PS print
    dispf('PRINTPDF %s with options=''append'': Postscript level 2 printing...',psmsg)
    print('-append',['-r' int2str(resolution)],'-dpsc2','',tmppsfile)
    % PS to PDF conversion (if required)
    if okps2pdf
        dispf('PRINTPDF %s with options=''append'': PS->PDF...',psmsg)
        if exist(fichier,'file'), delete(fichier), end
        ps2pdf('psfile', tmppsfile, 'pdffile', fichier, 'gspapersize', 'a4');
        if ~isnan(psoptions.pagestop) && ~isinf(psoptions.pagestop)  && (psoptions.pagenum>=psoptions.pagestop)
            delete(tmppsfile)
            fileinfo(fichier)
        end
    else
        dispf('PRINTPDF %s with options=''append'': conversion PS->PDF postponed...',psmsg)
        dispf('\t ps2pdf(''psfile'',''%s'',''pdffile'',''%s'',''gspapersize'',''a4'')\n',tmppsfile,fichier);
    end
    return
end

% printing with controls
printon = true;
ok = ~exist(fichier,'file');
while ~ok
    answer = questdlg({sprintf('the file ''%s'' already exist',[name ext]),sprintf('in ''%s''',chemin)},'Overwrite an existing PDF file','overwrite','new file','cancel','overwrite');
    if strcmp(answer,'new file')
        [fichier,chemin] = uiputfile('*.pdf','Choose a new filename for your pdf');
        [chemin,name,ext] = fileparts(fullfile(chemin,fichier));
        if ~strcmp(ext,'.pdf'), ext = ext_default; end
        fichier = fullfile(chemin,[name ext]);
        ok = ~exist(fichier,'file');
    else
        printon = ~strcmp(answer,'cancel');
        ok = true;
    end
end
if printon && ~strcmp(pmode,'auto')
    answer=questdlg({ sprintf('Current ''paperpositionmode'' is ''%s''',pmode)
            ' '
            'If you set it as ''auto'', you can change the paper layout'
            'by resizing the figure and you can check the result by refreshing'
            'the preview in the printpreview panel.'
            'With ''auto'' printing starts only after you close the the preview panel'
            'You do not need to change the options in the preview panel.'
            ' '
            'Select an option' },'PAPERPOSITIONMODE: auto or not ?','auto','no change','cancel','no change');
    if strcmp(answer,'auto')
        set(currentfig,'paperpositionmode','auto');
    elseif strcmp(answer,'cancel')
        printon = false;
    else
        dispf('no change, paperpositionmode: %s',pmode)
    end
end
 
if printon
    if get(currentfig,'paperpositionmode')
        disp('Printing starts only after you close the the preview panel...')
        uiwait(printpreview)
    end
    start = clock;
    dispf('printing ''%s''....',fichier)
    print(['-r' int2str(resolution)],'-dpdf',options, fichier)
    dispf('... end in %0.3g s',etime(clock,start))
    fileinfo(fichier)
    if isunix
       system(['evince ' fichier ' &']);
    else
        winopen(fichier)
    end
else
    disp('PRINT_PDF canceled')
end