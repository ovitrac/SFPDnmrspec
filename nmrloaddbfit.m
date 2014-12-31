function dbfit = nmrloaddbfit(varargin)
%NMRLOADDBFIT buidls fitting database of NMR spectra (spectral db created with nmrloadascii, nmrloadbase...)
%     Syntax: dbfit = nmrloaddbfit('parameter1',value1,'parameter2',value2,...)
%
%     List of input pair parameters/values:
%              'path': full path where db is saved
%            'dbname': name of db (default = dbfit.mat)
%       keywords
%          noprefetch: force to refresh database
%
% outputs : 
% dbfit.(mol).multipletpresent : multiplet types present in susbtance spectrum
% struct dbfit.(mol).count : numberof multiplet in each multiplet type
% struct dbfit.(mol).all : all data from all present multiplet 
% struct dbfit.(mol).(multipletpresent)(i,j).(fitfields) (i = number of a considered
% multiplet type in substance spectrum, j = 1 or 2= 2 fitting models)
%         'ppm': chemical shitf
%         'Iraw': raw intensity
%         'Irawn': normalized raw intensity by total number of proton of subs
%         'Ifit': intensity after fitting (see monotonepeakfit)
%         'Ifitn': normalized intensity after fitting by total number of proton of subs
%         'area': surface area of peak
%         'weight': weigth attributed to peak (see monotonepeakfit)
%         'relativeweight': relative weight attributed (see monotonepeakfit)
%         'width': width of peak, in ppm unit (see monotonepeakfit)
%         'position': position (center) of peak, in ppm unit (see monotonepeakfit)
%
% see also: monotonepeak, monotonepeakfit
% RMNSPEC v 0.1 - 12/04/13 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev.05/11/14
% history
% 15/4/13 add mfiltratio as option (default = [])
% 21/5/13 change data saved in output for nsignificantspeaks and significantproba
% 01/11/13 add fields dbfit.(all,s,t,d...).area (area of motif) 
%          add dbfit.(all,s,t,d...).Irawn (normalized Iraw by total number of H)
%          add dbfit.(all,s,t,d...).Ifitn (normalized Ifit by total number of H)
%          add db.calibration containing P and S (coeff of polynomial when area vs nH)
% 05/11/14 add fields dbitf.(mol).(all,s,t,d...).model : model of fitted peaks, in order to reconstitute peaks after if needed
%          add fields dbitf.(mol).(all,s,t,d...).npeaktofit: number of peak to be fitted
% default
default = struct('path',find_path_toolbox('rmnspec'),...
                 'dbname','dbfit.mat',...
                 'filenfo','mask.ods',...
                 'sheetnamenfo','nfodbfit');
keyword = 'noprefetch';

% Argcheck
o = argcheck(varargin,default,keyword);
% main
% check if dbfit exits
if exist(fullfile(o.path,o.dbname),'file') && ~o.noprefetch 
    dispf('DBFIT: reuse existing fitting database')
    fileinfo(fullfile(o.path,o.dbname))
    load(fullfile(o.path,o.dbname))
else
    [dbpur,dbxpur,dbmask,~,dbmultiplet] = nmrloadbase;
    substances = dbxpur.commonname; nsubstances = length(substances);
    dbfit = struct([]);
    % storage by multiplet type
    multiplettype = fieldnames(dbmultiplet); % all multiplet types as defined in mask.ods file (sheetname = 'multiplicity')
    fitfield = {'weight' 'relativeweight' 'width' 'position' 'r2max'};
    for isub = 1:nsubstances
        count = zeros(length(multiplettype),1);
        ROI = dbpur.(substances{isub}).gates; % extract all RoI (nx4 array: ppmmin, ppmmax, buffer and weight)
        ROI = ROI(ROI(:,4)>0,:); % remove parasite peak with ROI(:,4) = weight = 0
        proton = dbpur.(substances{isub}).H(ROI(:,4)>0);
        multiplicity = dbpur.(substances{isub}).multiplicity(ROI(:,4)>0); % multiplicity = multiplicity; 
        peaktofit = dbmask.(substances{isub}).peaksignificant(ROI(:,4)>0); % number of peaks to be fitted;
        mfiltratio = dbmask.(substances{isub}).mfiltratio(ROI(:,4)>0); % mfiltratio = mfiltratio;
        baseline = dbmask.(substances{isub}).baseline(ROI(:,4)>0); % baseline = baseline;
        for imultiplet = 1:length(multiplettype)
            present = cellfun(@(x) strcmp(x,multiplettype{imultiplet}), multiplicity); % check if a tested multiplet type is present
            ROIpresent = ROI(present,:); % ROI for a tested multiplet type            
            if ~isempty([multiplicity{present}])
                count(imultiplet,1) = length(multiplicity(present));
                multipletpresent = multiplicity{present};
                dbfit(1).(substances{isub}).(multipletpresent) = repmat(struct('ngaussian',[],'npeaktofit',[],'model',[],'ppm',[],'Iraw',[],'Ifit',[],...
                                                                 'weight',[],'relativeweight',[],'r2max',[],...
                                                                 'width',[],'position',[],'window',[]),length(multiplicity(present)),2);
                for j = 1:length(multiplicity(present))
                    valid = ((dbxpur.ppm>=ROIpresent(j,1))&(dbxpur.ppm<=ROIpresent(j,2))); % considered ROI
                    x = dbxpur.ppm(valid);
                    y = dbxpur.I(valid,isub); 
                    mfiltratiopresent = mfiltratio(present);
                    baselinepresent = baseline(present);
                    if  ~isnan(mfiltratiopresent(j,1))                        
                        p = monotonepeak('x',x,'y',y,'mfilt',mfiltratiopresent(j,1)*length(y),'array','sort','descend');
                    else p = monotonepeak('x',x,'y',y,'array','sort','descend');
                    end        
                    if baselinepresent(j,1)==0
                       [gaussianpeak,model] = monotonepeakfit(p,'x',x,'y',y,'significant',.9,'sort');
                    elseif baselinepresent(j,1)==1, [gaussianpeak,model] = monotonepeakfit(p,'x',x,'y',y,'significant',.9,'baseline','sort','endforced');
                    end

                    for k = 1:2 % 2 fitting models 
                        npeaktofit = peaktofit(present);
                        dbfit(1).(substances{isub}).(multipletpresent)(j,k).ngaussian = size(gaussianpeak,1);
                        dbfit(1).(substances{isub}).(multipletpresent)(j,k).npeaktofit = npeaktofit;
                        dbfit(1).(substances{isub}).(multipletpresent)(j,k).model = model;
                        dbfit(1).(substances{isub}).(multipletpresent)(j,k).ppm = x;
                        dbfit(1).(substances{isub}).(multipletpresent)(j,k).Iraw = y;
                        dbfit(1).(substances{isub}).(multipletpresent)(j,k).Irawn = y/sum(proton);
                        dbfit(1).(substances{isub}).(multipletpresent)(j,k).window = gaussianpeak(1,k).window;
                        if  ~isnan(npeaktofit(j,1)),                        
                            dbfit(1).(substances{isub}).(multipletpresent)(j,k).Ifit = model(x,npeaktofit(j,1),k);
                            dbfit(1).(substances{isub}).(multipletpresent)(j,k).Ifitn = model(x,npeaktofit(j,1),k)/sum(proton);
                            dbfit(1).(substances{isub}).(multipletpresent)(j,k).area = trapz(fliplr(x),model(x,npeaktofit(j,1),k));
                            dbfit(1).(substances{isub}).(multipletpresent)(j,k).window.nsignificantpeaks = npeaktofit(j,1);
                            dbfit(1).(substances{isub}).(multipletpresent)(j,k).window.significantproba = NaN;
                        else
                            npeaksignificant = gaussianpeak(1,k).window.nsignificantpeaks;
                            dbfit(1).(substances{isub}).(multipletpresent)(j,k).Ifit = model(x,npeaksignificant,k);
                            dbfit(1).(substances{isub}).(multipletpresent)(j,k).Ifitn = model(x,npeaksignificant,k)/sum(proton);
                            dbfit(1).(substances{isub}).(multipletpresent)(j,k).area = trapz(fliplr(x),model(x,npeaksignificant,k));
                        end
                        for field = fitfield;
                            dbfit(1).(substances{isub}).(multipletpresent)(j,k).(field{1}) = [gaussianpeak(:,k).(field{1})];      
                        end % next field           
                    end % next fitting model    
                end % next multiplet 
            end 
        end % next multiplet type   
        dbfit(1).(substances{isub}).multipletpresent = unique(multiplicity);
        dbfit(1).(substances{isub}).count = count(count>0);
        l=1;
        mul = unique(multiplicity);
        for imul = 1:length(mul)
            for jmul=1:size(dbfit(1).(substances{isub}).(mul{imul}),1)     
                dbfit(1).(substances{isub}).all(l,:) = dbfit(1).(substances{isub}).(mul{imul})(jmul,:);
                l=l+1;
            end
        end
        window = [dbfit.(substances{isub}).all(:,1).window];
        pos = [window.center]; [~,isort] = sort(pos,'descend');
        dbfit.(substances{isub}).all = dbfit.(substances{isub}).all(isort,:);
    
    % calibration by area=f(nH)
    area = [dbfit.(substances{isub}).all(:,1).area]';
    [p,S] = polyfit(proton,area,1);
    dbfit.(substances{isub}).calibration.p = p;
    dbfit.(substances{isub}).calibration.S = S;
    end % next substance
    
    % add help
    nfo = loadodsprefetch(fullfile(o.path,o.filenfo),'sheetname',o.sheetnamenfo);
    dbfit.help = cell2struct(nfo.description,nfo.fieldname);
    % save
    save(fullfile(o.path,o.dbname),'dbfit')
    fileinfo(fullfile(o.path,o.dbname))
end