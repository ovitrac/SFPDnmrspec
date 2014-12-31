function [dydxout,dydxicout] = diffp(x,y,ordre,h,n,D,xD,alpha,ilist,fusion,wy,proba)
%DIFFP évalue les dérivées nièmes à partir d'une approximation polynomiale (y = K*p) locale régularisée
%	syntaxe : dydx = diffp(x,y)
%	options : [dydx,dydxic] = diffp(x,y,ordre,h,n,D,xD,alpha,i,fusion,wy,prob)
%	Inputs
%	======
%		x       = entrées, dim(x) = [m 1]
%		y       = sorties, dim(y) = [m ny] avec ny>= 1
%		ordre   = ordre de la différenciation (par défaut = 0:1)
%				==> si ordre est de type cell : ordre{2} contient les ordres requis pour le calcul des intervalles de confiance
%		h       : x +/- h = intervalle de x pris en compte pour l'approximation polynomiale (par défaut = 10 * mean(diff(x)))
%		n       = degré du polynome utilisé pour l'approximation (par défaut = min(5,m/2))
%		D       = définit les contraintes sur les différentielles de y kiemes en x
%  	 	          exemple de prototype [3 0; 4 0; 5 0] pour définir des dérivées 3iemes, 4iemes et 5iemes nulles
%   	          dans le cas de derivées kiemes toutes nulles (la 2nde colonne peut être omise)
%			  	==> par défaut : D = [], ie: absence de contraintes
%		xD      = définit les points pour lesquels les contraintes sont valides (par défaut = en tous les points)
%		alpha   = valeur du paramètre de régularisation (par défaut = 1)
%		i       = index des points [xi yi] pour lequels une approximation est recherchée (par défaut = 1:m)
%		fusion	= indique si les colonnes de y doivent être fusionnées lors de l'assemblage final (par défaut = 1)
%		wy		= poids respectif pour la fusion (par défaut = repmat(1,1,my))
%		proba	= seuil de probabilité pour le calcul des intervalles de confiance dydxic
%	Outputs
%	=======
%	dydx	= dérivées nièmes calculées
%	dydxic	= intervalles de confiance

	
% Thermique 1.1 - 13/06/01 - Olivier Vitrac - rev. 15/10/2013

% Revision history
% 28/06/03 weighting and symmetry
% 28/04/04 adaptable window expansion
% 28/04/04 bug: post xiamp correction
% 12/07/2005 additional modifications
% 23/11/2006 fix bug with several inputs (columns)
% 02/06/2010 fix multiple Y
% 24/07/2010 optimized code for 2009x rules (2010 rules not implemented to keep some backward-compatibility)
% 15/10/2013 implements variable h

% exemple : évaluation de la dérivée 1ère de sin(t) entre 0 et 10pi bruitée (à 14f)
% t = 0:0.1:10*pi;
% y = sin(t)+0.05*sin(14*t);
%
% sans régularisation
% dydt = diffp(t,y,1);
% figure(1), plot(t,y,'k-',t,cos(t),'b-',t,dydt,'r:')
%
% avec régularisation : dérivée troisième nulle
% [dydt,dydti] = diffp(t,y,1,[],[],3);
% figure(2), plot(t,y,'k-',t,cos(t),'b-',t,dydt,'r-',t,dydt-dydti,'r:',t,dydt+dydti,'r:')
%
% récaptitulatif :
% en noir : primitive
% en bleu : dérivée exacte
% en rouge ! dérivée calculée
% en rouge pointillé : intervalle de confiance de la dérivée calculée

%	=====================================================================================================
%		                 Principe de la méthode utilisée (procédure récursive)
%	=====================================================================================================
% Rappel sur l'approximation polynomiale : y = K*p
%		y = valeurs à prédire, dim(y) = [m 1]
%		p = polynome de degré n (Lagrange, Hermite)
%		K = noyau d'interpolation (généralisation de la matrice de la matrice de Vandermonde) relatif aux données x, dim(K) = [m n] avec Lagrange
%		identification par une méthode des moindres carrés (LS) : arg min (||y-p*X||2)^2
%		REM : Le choix d'une norme L2 (euclidienne) permet d'utiliser des méthodes d'inversion linéaire (i.e. espace de Hilbert)
% Définition d'une approximation locale en xi de yi = f(x) à partir d'une méthode LS pondérée (WLS)
% 		soit Wi = noyau de pondération, dim(Wi) = [m m] 
%		LS utilisant une norme pondérée Wi : arg min ((||y-pi*X||Wi)^2) = arg min ((||Wi(y)-pi*(Wi*X)||2)^2)
% Méthode de régularisation (Tychonov)
%		problème régularisé : arg min ( (||Wi(y)-pi(Wi*X)||2)^2 + alpha^2 * (||Li*pi-di||2)^2 )
%		Li = opérateur de différenciation de p relatif à xi, dim(Li) = [nL n] où nL est le nombre de contraintes
%		di = valeur(s) cible(s) de la différencielle Li*pi
%		alpha = scalaire positif contrôlant la régularisation
%	Assemblage du problème complet
%    
%             ||                                       ||2
%             || [   Wi*K   ]           [   Wi*y   ]   ||
%   arg min   || [          ] * pi   -  [          ]   ||  
%             || [ alpha*Li ]           [ alpha*di ]   ||
%             ||                                       ||2
% 				  =========			 	==========
%                     Ai					       Bi
%	soit arg min (||Ai*p-Bi||2)^2 qui est résolut par pi = pinv(Ai)*Bi pour chaque point xi
%
% Remarque : le choix judicieux de Wi permet généralement de limiter raisonnablement le nombre de points
% voisins retenu et réduit ainsi la dimension de K, ainsi que la complexité de la procédure d'inversion.
%
% Remarque2 : le pseudoinverse (pinv ou pinvs) de A est calculé sur la base d'une décomposition
% en vecteurs singuliers de K, seuls les vecteurs associés aux valeurs singulières les plus élevées
% sont utilisés pour le calcul de la solution finale (eq. filtre passe-bas).


% Constantes
ordre_defaut	= 0:1;
alpha_defaut	= 1;
degre_defaut	= 5;
nhdefaut		= 10;
fusion_defaut	= 0;
proba_defaut	= 0.95;
it_max          = 100;   % max iteration for expansion
hexpand         = 1.4;   % expansion factor (max expansion = it_max*hexpand)

% arg check
h_defaut = nhdefaut * mean(diff(x));
if nargin<4, h = h_defaut; end

% Symmetric boundary conditions (added 9/06/03)
x = x(:);
[my,ny]	= size(y);
if my==1 && ny>1,y= y'; [my,ny] = size(y);end %#ok<NASGU>
if nargin<3, hsym = nhdefaut * max(diff(x));
elseif isempty(h), hsym = nhdefaut * max(diff(x));
else hsym = h;
end
if numel(hsym)>1, hsym=max(hsym(:)); end
leftBC  = 2:find((x-x(1))<hsym, 1, 'last' );      mLBC = length(leftBC);
rightBC = find((x(end)-x)<hsym, 1, 'first'):my-1; mRBC = length(rightBC);
ioriginal = mLBC+(1:length(x));
x = [   2*x(1)-flipud(x(leftBC))
        x
        2*x(my)-flipud(x(rightBC)) ];
y = [   2*repmat(y(1,:),mLBC,1)-flipud(y(leftBC,:))
        y
        2*repmat(y(my,:),mRBC,1) - flipud(y(rightBC,:))
    ];

    
% Contrôle des entrées : diffp(x,y,ordre,h,n,D,alpha,i)
							x			= x(:);
							m			= length(x);
							[my,ny]	= size(y);
if my==1 && ny>1, 			y 			= y'; [my,ny] = size(y);					end
							h_defaut = nhdefaut * mean(diff(x));
							%defineKon = 0; %unused
if m~=my, error('x and y must have the same number of rows.'),						end
if nargin<3,			ordre		= ordre_defaut;									end
if nargin<4,			h			= h_defaut;										end
if nargin<5,			n			= degre_defaut;									end
if nargin<6,			D			= [];											end
if nargin<7,			xD			= []; else xD = xD(:);							end
if nargin<8,			alpha		= alpha_defaut;									end
if nargin<9,			ilist		= [];                                           end
if nargin<10,			fusion		= fusion_defaut;								end
if nargin<11, 			wy 			= ones(1,ny);                                   end
if nargin<12,			proba		= proba_defaut;									end
if isempty(ordre),		ordre		= ordre_defaut;									end
if isempty(h), 			h 			= h_defaut;										end
if isempty(n), 			n			= degre_defaut;									end
if isempty(alpha),		alpha		= alpha_defaut;									end
if isempty(ilist),		ilist		= ioriginal; else ilist = ioriginal(ilist);     end
if isempty(fusion),		fusion		= fusion_defaut;								end
if isempty(proba),		proba		= proba_defaut;									end
if length(n)>1,			error('n must be a scalar'),								end
						n 			= min(n,round(m/2));
						wy		= sparse(wy/sum(wy));
						compute_IC	= nargout>1;
if n<0 || n~=fix(n), error('n must be a positive integer'),							end

% variable h (15/10/2013)
if  length(h)>1
    nh = numel(h);
    if nargout>1, error('only the first output can be returned with vectorized h values'), end
    dydxout = struct( 'diffp',{repmat({[]},size(h))},...
                      'hmin',[],'diffphmin',[],...
                      'hmax',[],'diffphmax',[]);
    for i=1:nh
        try
            dydxout.diffp{i} = diffp(x,y,ordre,h(i),n,D,xD,alpha,ilist,fusion,wy,proba);
        catch
            dispf('DIFFP:: ERROR for h(%d)=%0.4g',i,h(i))
        end
    end
    [~,ind] = sort(h,'ascend');
    imin = find(~cellfun(@isempty,dydxout.diffp(ind)),1,'first');
    imax = find(~cellfun(@isempty,dydxout.diffp(ind)),1,'last');
    dydxout.diffphmin = dydxout.diffp{ind(imin)}; dydxout.hmin = h(imin);
    dydxout.diffphmax = dydxout.diffp{ind(imax)}; dydxout.hmax = h(imax);
    return
end

% Normalisation des contraintes	
if ~iscell(D), D = {D}; end
if ~iscell(xD), xD = {xD}; end
nDc = length(D);
if nDc~=length(xD), error('incompatible constrains'), end
nL	= zeros(nDc,1);
nxD	= zeros(nDc,1);
for i=1:nDc
	[nLi,k]	= size(D{i});
	nL(i)	= nLi;
	nxD(i)	= length(xD{i});
	if k==1,
		D{i} = [D{i} zeros(nLi,1)];
	elseif k>2
		error('D must be a nLx1 or nLx2 matrix')
	end
end
defineLon 	= any(nL);

% Définition des ordres pour le calcul des intervalles de confiance
if iscell(ordre)
	if length(ordre)>=2
		ordreic = ordre{2};
		ordre = ordre{1};
	else
		ordre = ordre{1};
		ordreic = ordre;
	end
else
	ordreic = ordre;
end

% Initialisation
%nilist	= length(ilist); %unused
nordre	= length(ordre);
nordreic= length(ordreic);
if fusion
	dydxout	= zeros(length(ilist),nordre);
	if compute_IC, dydxicout = zeros(length(ilist),nordreic); end
else
	dydxout	= zeros(length(ilist),nordre*ny);
	if compute_IC, dydxicout = zeros(length(ilist),ny*nordreic); end
end

% Fenetrage
isol = 0;
for i = ilist
    isol = isol + 1;
	% Assemblage des poids W[i]
	hi = h;
    it = 0;
    mi = 0;
    while (mi<1.5*n) && (it<it_max)
        if it==1, warning('DIFFP: sampling interval expanding for x=%0.3g',x(i)), end %#ok<WNTAG>
        if it>0, hi = hi * hexpand; end
        it = it + 1;
        choisir_voisins = 1;
        while choisir_voisins
            ui		= (x-x(i))/hi;
            wi		= max(	[1-abs(ui).^3,zeros(m,1)], [], 2 ).^3;					% [max(1-|u|^3,0)]^3
            indi = find(wi>0);														% points voisins
            ng		= i-indi(1);													% nb voisins à gauche de i
            nd		= indi(end)-i;												    % nb voisins à droite de i
            if (abs(nd-ng)>1) && (hi==h), hi = 2*hi*( 1 - min(ng,nd)/(ng+nd) );	        % élargissement
            else choisir_voisins = 0; end
        end
        mi = length(indi);															% nombre de points voisins retenus
    end
    if it==it_max
        error('DIFFP: bad weighting (likely reason: too high discrepancy in sampling rate)')
    end
	if fusion
		Wi = spdiags(repmat(wi(indi),ny,1),0,mi*ny,mi*ny);						% matrice creuse diagonale diag(Wi) = wi(indi)
	else
		Wi = spdiags(wi(indi),0,mi,mi);											% matrice creuse diagonale diag(Wi) = wi(indi)
	end
    ximid = x(i);
    xiamp = abs(x(indi(end))-x(indi(1)));
    xi   = sparse(x(indi));
	xired = sparse((xi-ximid)/xiamp);                                           % normalized abscissae
	yi = sparse(y(indi,:));
	% Assemblage de K (noyau d'interpolation)
	% les exposants sont classés dans l'ordre décroissant (conformément au choix retenu dans Matlab)
	Ki = diffker(xired,n,0);														% 0 = primitive
	% Assemblage des contraintes sur les dérivées kièmes : L[i]
	nxD_find = 0;
	if defineLon
		xD_find		= cell(nDc,1);
		nxD_find	= zeros(nDc,1);
		for ic = 1:nDc
			if ~nxD																% en l'absence de points spécifiques = en tout point
				xD_find{ic} = (x(i)-ximid)/xiamp;
				nxD_find(ic) = 1;
			else																% pour des points spécifiques
				xD_find{ic} = (intersect(xD{ic},xi)-ximid)/xiamp;
				nxD_find(ic) = length(xD_find{ic});
			end
		end
	end
	Li = sparse([]);
	di = sparse([]);
	icDfind = find(nxD_find);
	icDfind = icDfind(:)';
	if any(icDfind)
		for ic = icDfind
			Li = [ Li ; diffker(xD_find{ic} ,n,D{ic}(:,1)) ]; %#ok<AGROW>
			di = [ di ; sparse(repmat(D{ic}(:,2),nxD_find(ic),1)) ]; %#ok<AGROW>
		end
	end
	% Assemblage du problème complet
	if fusion
		wyKBi		= reshape(repmat(wy,mi,1),mi*ny,1);
		if nxD_find
			wyLBi	= reshape(repmat(wy,nxD_find,1),nxD_find*ny,1);
		else
			wyLBi	= sparse([]);
		end
		wyKAi		= repmat(wyKBi,1,n+1);
		wyLAi		= repmat(wyLBi,1,n+1);
        % before 23/11/06
		%Ai			= [	(Wi*repmat(Ki,ny,1))		.*wyKAi	;	repmat(alpha*Li,ny,1) .*wyLAi	];
		%Bi			= [	(Wi*reshape(yi,mi*ny,1))	.*wyKBi	;	repmat(alpha*di,ny,1) .*wyLBi	];
        % after 23/11/06
        %Ai			= [	(Wi*repmat(Ki,ny,1))		.*wyKAi	;	repmat(alpha*Li,ny,1) .* repmat(wyLAi,ny,1)	];
		%Bi			= [	(Wi*reshape(yi,mi*ny,1))	.*wyKBi	;	repmat(alpha*di,ny,1) .* repmat(wyLBi,ny,1)	];
        % after 02/06/10
        Ai			= [	(Wi*repmat(Ki,ny,1))		.*wyKAi	;	repmat(alpha*Li,ny,1) .* repmat(wyLAi,nL,1)	];
		Bi			= [	(Wi*reshape(yi,mi*ny,1))	.*wyKBi	;	repmat(alpha*di,ny,1) .* repmat(wyLBi,nL,1)	];
	else
		neq			= size(Ki,1) + size(Li,1);
		Ai			= [ Wi*Ki			; alpha*Li ];
		Bi			= zeros(neq,ny);
		for iy = 1:ny
			Bi(:,iy)	= [ Wi*yi(:,iy) ; alpha*di ];
		end
	end
	% Identification du polynome ou des polynomes (si ny>1 et ~fusion)
	% Résolution du système linéaire avec pinvs (cf. svds pour + d'info)
	% ==> les coefficients sont classés columwise
	pi			= pinvs(Ai)*Bi;
	% Calcul des dérivées requises en x[i]
	ki = diffker(0,n,ordre);	% noyau d'interpolation local
	dydx = ki*pi;
    % before 23/11/06
    %exponent = repmat(ordre(:),1,ny);
    % end before 23/11/06
	if fusion
		dydx = reshape(dydx,size(dydx,1)/nordre,nordre);
        % before 23/11/06
        %exponent = reshape(exponent,size(exponent,1)/nordre,nordre);
        % after 23/11/06
        exponent = ordre(:)';
    else
        % after 23//11/06
        exponent = repmat(ordre(:),1,ny);
        % end after 23//11/06
		dydx = reshape(dydx,size(dydx,1)/nordre,ny*nordre);
        exponent = reshape(exponent,size(exponent,1)/nordre,ny*nordre);
	end
    dydx = dydx ./ (xiamp .^ exponent);
    
	% Intervalle de confiance (algo adapté de polyconf)
	if compute_IC
		ki		= diffker(0,n,ordreic);	% noyau d'interpolation local
		yie		= repmat(Ki*pi,1,ny);	                % valeurs prédites
		normr	= sqrt((yi(:)-yie(:))'*(yi(:)-yie(:))); % normes des résidus
		ddl		= max(1,mi*ny-n-1);                 	% ddl
		[c,R]	= qr(Ai);	                            %#ok<ASGLU> % facteur de Cholesky de la matrice d'interpolation (R'R = Ai'*Ai)
		E		= (pinvs(R')*ki')';	                    % modification de R (équivalent ki/R)
		dn      = size(R,1)-size(Wi,1);
		if dn>0, Wi = [Wi;zeros(dn,mi*ny)]; end %#ok<AGROW>
		E       = E*Wi;
		e		= sqrt(1+sum(E.*E,2)'); % replace sqrt(1+sum((E.*E)')');
		% si tinv disponible
		dydxic	= (normr/sqrt(ddl))*e*tinv(1-(1-proba)/2,ddl);
		% sinon
		%dydxic	= (normr/sqrt(ddl))*e*pnorm_inv(1-(1-proba)/2);
        exponent = repmat(ordreic(:),1,ny);
		if fusion
			dydxic = reshape(dydxic,size(dydxic,1)/nordreic,nordreic);
            exponent = reshape(exponent,size(exponent,1)/nordreic,nordreic);
		else
			dydxic = reshape(dydxic,size(dydxic,1)/nordreic,ny*nordreic);
            exponent = reshape(exponent,size(exponent,1)/nordreic,ny*nordreic);
		end
        dydxic = dydxic ./ (xiamp .^ exponent);
	end
	dydxout(isol,:) = dydx;
	if compute_IC
		dydxicout(isol,:) = dydxic;
	end
end