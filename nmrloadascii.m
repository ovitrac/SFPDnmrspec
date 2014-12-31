function [db,dbx] = nmrloadascii(varargin)
%NMRLOAD load a collection of NMR spectra set in an ODS file exported as ASCII files (2D tables ppm itensity coded as '%f %f')
%     Syntax: db = nmrloadascii('parameter1',value1,'parameter2',value2,...)
%
%     List of input pair parameters/values:
%              'path': full path where ODS and TXT files are located (default = fullfile(find_path_toolbox('rmnspec'),'datapur'))
%           'odsfile': ODS file listing available spectra (default='mask.ods')
%         'sheetname': worksheet containing the table listing code of molecules and their properties (default='substance')
%      'sheetnamenfo': worksheet containing help for db (default='nfo_substance'}
%        'primarykey': name of the column listing codes to be used as fieldnames in Matlab (default = 'commonname')
%      'spectrumfile': name of the column listing filenames storing spectra in TXT files (default = 'reference')
%      'prefetchfile': file to be used for fast loading (default=odsfile_sheetname.mat)
%            'ppmmin': truncate ppm values lower than ppmmin (default = -1)
%            'ppmmax': truncate ppm vaules greater than ppmmax (default = 13)
%        'resolution': number of points considered in the RMN spectrum (default = 2^15)
%            'method': interpolation method (default = 'cubic'), see INTERP1
%         'extrapval': extrapolation value (default = 'extrap'), see INTERP1
%      'ppmstandards': nx2 array to remove segments in spectra corresponding to standards (TMS, DOH, D-chlrophorm), where n is the number of standards
%                  default = [-0.1 0.2; 1.53 1.575; 7.2 7.3];
% 'ppmbufferbaseline': buffer in ppm to identify the begining and end of peaks (default=0.01 ppm, see mfilt in nmrbaseline)
% 'ppmbufferstandard': buffer in ppm to be used for removing standards (default=0.001ppm)
%         'ppmbuffer': buffer in ppm to identify the begining and end of peaks (default=0.0076 ppm (corresponding mfilt=20), see mfilt in monotone2peaks)
%          'maskfile': ods file containing mask information (default = mask.ods)
%  'maskprefetchfile': mask file to be used for fast loading (default = dbmask.mat)
%
%     Output: db=structure with fields:
%     db.(primarykey).spectrumfile: name of the file storing the spectrum (without extension .TXT)
%     db.(primarykey).ppm chemical shifts
%     db.(primarykey).I corresponding intensity (raw values)
%     db.(primarykey).In intensity normalized so that energy is 1 and average baseline value is 0
%     db.(primaryley).normvalue: value used to normalize intensity (In = I/normvalue)
%     db.(primarykey).In0 as In but centered on its average
%     db.(primarykey).field1, db.(primarykey).field2,...
%     where field1, field2,... are columns content in odsfile
%     db.(primarykey).filename: filename of raw spectrum
%     db.(primarykey).mfiltbaseline: filtering width to indentify baseline
%     db.(primarykey).mfiltpeaks: filtering width to indentify peaks
%     db.(primarykey).If: filtered I with mfilt
%     db.(primarykey).p: posotion of siginficant peaks 
%     db.(primarykey).a: amplitude of peaks identified
%     db.(primarykey).rr: relative rank of the peak
%     db.(primarykey).pstart: position of the beginning of peaks
%     db.(primarykey).pstop: position of the end of peaks
%     db.(primarykey).preject: rejection percentage of peaks
%     db.(primarykey).areject: corresponding amplitude to preject
%     db.(primarykey).Id: Intensity for discretization 
%     db.(primarykey).gates: nx4 array from dbmask data for MULTIPB ([ppmmin ppmmax buffer weight]) 
%     db.(primarykey).weight: nx1 array of weight attributed by user in dbmask
%     db.(primarykey).weightn: nx1 array of normalized weight attributed by user in dbmask (weight/sum(weight)
%     db.(primarykey).H: nx1 cell of proton number, attributed by user in dbmask
%     db.(primarykey).Hn: nx1 cell of normalized proton number, attributed by user in dbmask (H/sum(H)
%     db.(primarykey).W: weighting function (for 'pur' database)
%     db.(primarykey).Iw: weighted intensity (for 'pur' database)
%     db.(primarykey).calibration: structure, containing coeff of a polynomial (P and S of polyfit) area = f(number of H)
%     db.(primarykey).InH: normalized intensity by total number of H of molecule (In/(pente*nH total)
%
%   Advanced use:
%       nmrloadascii(...,'noprefetch') forces the prefetch file to be refreshed
%       [db,dbx] = nmrloadascii(...) or [~,dbx] = nmrloadascii(...)
%       dbx: database for mathematical treatment with fields
%       dbx.filename: full path of raw spectra data files
%            dbx.ppm: mx1 array of ppm values (chemical displacement)
%           dbx.step: ppm step
%              dbx.m: number of ppm values
%              dbx.n: number of spectra
%              dbx.I: mxn array of normalized intensity values (as In in db)
%             dbx.I0: mxn array of centered and normalized intensity values (as In0 in db)
%  dbx.mfiltbaseline: filtering width to indentify baseline
%     dbx.mfiltpeaks: filtering width to indentify peaks
%             dbx.If: mxn filtered I with mfilt
%        dbx.preject: nx1 rejection percentage of peaks
%        dbx.areject: nx1 corresponding amplitude to preject
%             dbx.Id: mxn intensity for discretization  
%             dbx.Iw: mxn weighted intensity (for 'pur' database)
%              dbx.H: total number of H (for 'pur' database)
%
% See also: nmrcorr, nmrsubdb


% RMNSPEC v 0.1 - 31/08/12 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev. 01/11/12

% History
% 02/09/12 fix dbx.I on centered values
% 10/09/12 fix db.(...).ppm, ppmstandards
% 10/09/12 fix dbmol (ODS file) so that it is sorted as data files are
% 10/09/12 add field filename, to db and dbx
% 11/09/12 fix ppmstandards
% 11/09/12 add nmrbaseline to fit baseline of nrm spectra (add o.ppmbuffer). Add dbx.I = db.(...).I and dbx.I0 = db.(...).In0
% 30/09/12 first implementation of sgmt2pb, TODOLIST to be implemented
% 01/10/12 add o.ppmbufferstandard (used in sgmt2pb) for removing standards
%          add o.ppmbufferbaseline, fix prejection calculated by monotone2peaks according to mfilt (used in nmrbaseline for identifying baseline)        
% 02/10/12 add o.ppmbufferpeaks (different than o.ppmbufferbaseline, used for indentifying peaks)
%          add fields:p, a, rr, pstart, pstop, preject, areject, If, mfilt, Id in DB and DBX (as output fields of monotone2peaks) 
% 18/03/13 change permanently ppmmax = 13
%          add nmrloadmask, multipb for weighting nmr spectra (o.maskfile and o.maskprefetchfile)
%          ass keyword 'pur' for activate mask database (so that loading mixture database doesn't need mask)
%          add field mol.W (weight), mol.Iw (intensity weighted), dbx.Iw, dbx.W  
% 27/03/13 merge ods file 'substance' with mask.ods containing worksheet 'susbtance' -> change ods file name

% 04/04/13 add help
%          fix the mode to remove standards (variable according to mask file, parasite zone)
% 25/10/13 add fields proton and weight in dbpur (dbpur.H and dbpur.weight)
% 30/10/13 add fields normalized proton and normalized weight (dbpur.Hn,dbpur.weightn)
% 01/11/13 add field help for description of fields in db 
%          add db.(mol).calibration (P and S of polyfit), db.(mol).InH if dbfit is defined--> NON)
% 22/01/14 add db.(mol).normvalue

% default
default = struct(...
                'path',fullfile(find_path_toolbox('rmnspec'),'data_pur'),...
                'odsfile','mask.ods',... % substance.ods
                'sheetname','substance',...
                'sheetnamenfo','nfosubstance',...
                'primarykey','commonname',...
                'spectrumfile','reference',...
                'prefetchfile','',...
                'ppmmin',-1,...
                'ppmmax',13,... previous value: 13.9995
                'resolution',2^15,...
                'method','cubic',...
                'extrapval','extrap',...
                'ppmstandards',[-0.1 0.2; 1.53 1.575; 7.22 7.3; ],...
                'ppmbufferstandard',0.001,...
                'ppmbufferbaseline',0.01,...
                'ppmbufferpeaks',0.0076,...
                'maskfile','mask.ods',...
                'maskprefetchfile','dbmask.mat'); %'dbfit',[]);
                
keyword = {'noprefetch' 'pur'}; % 'pur' = load substance database, load mask alos 'calib'

% Argcheck
o = argcheck(varargin,default,keyword);
if isempty(o.prefetchfile), o.prefetchfile = sprintf('%s_%s.mat',o.odsfile,o.sheetname); end
% if o.calib && isempty(o.dbfit)
%     dbfit = nmrloaddbfit('path',fullfile(find_path_toolbox('rmnspec'),'data_pur'));
%     fitnfo = fileinfo(fullfile(find_path_toolbox('rmnspec'),'data_pur','dbfit.mat'));
%     fprintf('NMRLOADASCII:\t DBFIT is not defined, use following base:\n%s\ndate and bytes:%s%s',strcat(fitnfo.path,fitnfo.filename),fitnfo.date,fitnfo.bytes)
% end
ppm = linspace(o.ppmmin,o.ppmmax,o.resolution)';
mfiltbaseline = max(1,ceil(o.ppmbufferbaseline*o.resolution/(o.ppmmax - o.ppmmin)));
mfiltpeaks = max(1,ceil(o.ppmbufferpeaks*o.resolution/(o.ppmmax - o.ppmmin))); % mfilt for identifying peaks -> dbx.If / dbx.I0f
% Check if prefetch must be used
if exist(fullfile(o.path,o.prefetchfile),'file') && ~o.noprefetch    
    dispf('NMRLOADASCII:\t resuse the following prefetch file')
    fileinfo(fullfile(o.path,o.prefetchfile))
    load(fullfile(o.path,o.prefetchfile))

else % regenerate data
    dbmol = loadodsprefetch(fullfile(o.path,o.odsfile),'sheetname',{o.sheetname o.sheetnamenfo}); % load ods file of substances in database
    file = explore('*.txt',fullfile(o.path),[],'abbreviate'); % search all files .txt in local dir
    nfile = length(file);
    screen = '';
    % load mask if 'pur'
    if o.pur, dbmask = nmrloadmask('path',o.path,'odsfile',o.maskfile,'prefetchfile',o.maskprefetchfile); end    
    for i=1:nfile
        mol = bykeywords(dbmol.(o.sheetname),o.spectrumfile,file(i).name); % all info linked with "mol reference"
        % check/valid entry in table
        for f=fieldnames(mol)', if iscell(mol.(f{1})), mol.(f{1}) = mol.(f{1}){1}; end, end
        molname = deblank(mol.(o.primarykey)); % search name of molecule
        filename = fullfile(file(i).path,file(i).file);
        % load raw txt file
        screen = dispb(screen,'[%d/%d]\tloading spectrum of %s...',i,nfile,molname);
        fid = fopen(filename);
        raw = textscan(fid,'%f %f');
        fclose(fid);
        % store data in a structure
            % weight to remove standards
        if o.pur
            standards = dbmask.(molname).gates;
            standards = standards(standards(:,4)==0,:); % parasite peak with gates(:,4) = weight = 0
            ppmstandards = [standards(:,1) standards(:,2)]; % get ppmmin and ppmmax 
            mask = sgmt2pb(ppm,ppmstandards,'buffer',o.ppmbufferstandard); % put a specific buffer here
        else
            mask = sgmt2pb(ppm,o.ppmstandards,'buffer',o.ppmbufferstandard); % put a specific buffer here
        end
        mol.filename = filename;
        mol.ppm = ppm;
        mol.I  = interp1(raw{1},raw{2},ppm,o.method,o.extrapval); % raw intensity
        baseline = monotone2peaks(mol.I,'mfilt',[1 mfiltbaseline]);% research preject for fitting a baseline
        [isbase,Ib] = nmrbaseline(mol.I,'mfilt',mfiltbaseline,'prejection',baseline(2).preject); % fit a baseline
        mol.normvalue = trapz(ppm,mol.I.^2);
        mol.In = mask.*Ib/mol.normvalue;
        mol.mfiltbaseline = mfiltbaseline;
        peaks = monotone2peaks(mol.In,'mfilt',[1 mfiltpeaks]); % research significant peaks for discretization
        mol.mfiltpeaks = peaks(2).mfilt; % mfilt
        mol.If = peaks(2).If; % filtered I with mfilt
        mol.p = peaks(2).p; % posotion of peaks
        mol.a = peaks(2).a; % amplitude of peaks
        mol.rr = peaks(2).rr; % relative rank of the peak
        mol.pstart = peaks(2).pstart; %position of the beginning of peaks
        mol.pstop = peaks(2).pstop; % position of the end of peaks
        mol.preject = peaks(2).preject; % rejection percentage of peaks
        mol.areject = peaks(2).areject; % corresponding amplitude to preject
        % Intensity for discretization (Id)
        mol.Id = zeros(o.resolution,1);
        mol.Id(peaks(2).p) = mol.In(peaks(2).p);
        % remove standards
        % old method: for j=1:size(o.ppmstandards,1); mol.In((ppm>=o.ppmstandards(j,1)) & (ppm<=o.ppmstandards(j,2))) = 0; end        
        mol.In0 = mol.In-mean(mol.In);
        mol.isbase = isbase;
        % add weight for nmr peaks
        if o.pur
            mol.gates = dbmask.(molname).gates;
            mol.weight = dbmask.(molname).weight; % array of weight corresponding to band
            mol.weightn = dbmask.(molname).weight/sum(dbmask.(molname).weight); % array of normailzed weight
            tmpH = dbmask.(molname).proton; tmpH(isnan(tmpH)) = 0;            
            mol.H = tmpH ; % number of proton
            mol.Hn =  tmpH/sum(tmpH); % normalized number of proton
            mol.multiplicity = dbmask.(molname).multiplicity;
            mol.W = multipb(ppm,mol.gates); %function of weight
            mol.Iw = mol.In.*mol.W; % intensity weighted
%             if o.calib
%                 mol.calibration = dbfit.(molname).calibration; %calibration [p,S]=polyfit(area,H)
%                 mol.InH = mol.In/(dbfit.(molname).calibration.p(1)*sum(tmpH)); %In/(pente*nHtotal)
%             end
        end
        db.(molname) = mol;
    end
    % add help
    db.help = cell2struct(dbmol.(o.sheetnamenfo).description,dbmol.(o.sheetnamenfo).fieldname);
    
    % save
    dispf('NMRLOADASCII:\t new/updated prefetch file')
    save(fullfile(o.path,o.prefetchfile),'db','dbmol')
    fileinfo(fullfile(o.path,o.prefetchfile))
end

% Additional output: dbx
if nargout>0
    fdb = fieldnames(rmfield(db,'help'));
    dbx.filename = cellfun(@(m) db.(m).filename,fdb,'UniformOutput',false); % for traceability
    dbx.ppm = ppm;
    dbx.step = ppm(2)-ppm(1);
    dbx.m = o.resolution; % number of points (32768)
    dbx.n = length(fdb); % number of substances in database
    tmp = cellfun(@(m) db.(m).In,fdb,'UniformOutput',false);
    dbx.I = cat(2,tmp{:});
    tmp = cellfun(@(m) db.(m).In0,fdb,'UniformOutput',false);
    dbx.I0 = cat(2,tmp{:});
    dbx.mfiltbaseline = mfiltbaseline;
    dbx.mfiltpeak = mfiltpeaks;
    tmp = cellfun(@(m) db.(m).If,fdb,'UniformOutput',false);
    dbx.If = cat(2,tmp{:});
    tmp = cellfun(@(m) db.(m).preject,fdb,'UniformOutput',false);
    dbx.preject = cat(2,tmp{:});
    tmp = cellfun(@(m) db.(m).areject,fdb,'UniformOutput',false);
    dbx.areject = cat(2,tmp{:});
    tmp = cellfun(@(m) db.(m).Id,fdb,'UniformOutput',false); 
    dbx.Id = cat(2,tmp{:});
    if o.pur
        tmp = cellfun(@(m) db.(m).Iw,fdb,'UniformOutput',false); 
        dbx.Iw = cat(2,tmp{:});
        tmp = cellfun(@(m) db.(m).W,fdb,'UniformOutput',false); 
        dbx.W = cat(2,tmp{:});
        dbx.H = cell2mat(cellfun(@(m) sum(db.(m).H),fdb,'UniformOutput',false)); 
%         if any(strcmp('InH',fieldnames(db.(fdb{1}))))
%             cellfun(@(m) db.(m).InH,fdb,'UniformOutput',false); 
%             dbx.InH = cat(2,tmp{:});
%         end
    end
    % propagates fields supplied in ODS file
    primarykeysindb = cellfun(@(m) db.(m).(o.primarykey),fdb,'UniformOutput',false);
    dbx = orderfields(catstruct(bykeywords(dbmol.(o.sheetname),o.primarykey,primarykeysindb),dbx));
end