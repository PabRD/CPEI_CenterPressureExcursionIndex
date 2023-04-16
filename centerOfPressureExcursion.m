%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                      NAME: Traitement CPEI                              %
%                      AUTHOR: PabDawan                                   %
%                      DATE: Octobre 2022                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Description: %% Center of Pressure Excursion Index (CPEI) calculation
% PabRD | Octobre-Novembre 2022
% Outils: PressMat, TexiSense (64*64 capteurs), ~30Hz
tic
clc                                                                         %Clear Command Window
clear                                                                       %Clear Workspace
close all                                                                   %Close figure


%% Chemin d'accès
start_path=pwd;
allFiles = genpath(pwd);
allFiles = split(allFiles,';');

%% Extraction des données
% L'utilisateur sélectionne 5,10 ou 15 fichiers bruts.
initialPath=allFiles{2};
cd(initialPath)
[fileName,pathName,~] = uigetfile(strcat(allFiles{2},'\*.*'),'Select File','MultiSelect', 'on');

if ~iscell(fileName)
    fileName = {fileName};
end

nbFiles = size(fileName,2);

dataRaw=cell(1,15);
CPEI=ones(nbFiles,1);  

for kk = 1 : nbFiles
    str = fileName{kk};
    nameCond = regexp(str,'((?<=pG_).*(?=\.bin))','match');                 %trouver l'expression en sandwich de 'pG_' et '.bin'
    name = nameCond{1}(isstrprop(nameCond{1},'alpha'));
    prepCond = regexp(nameCond,'\d*','Match');
    cond = str2double(prepCond{1}{1});

%% organiser la matrices (64*64)
    fid = fopen(str,'rb');
    C = fread(fid);
    fclose(fid);
    clear fid
    
    
    stock_t=C(1:4097:numel(C));                                             %les intervalles de temps entre chaque frames sont stockés à la premiere case de chaque frame (64*64=4096)
    ncol = size(C, 1);
    C(1:4097:numel(C))=[];                                                  % On supprime les valeurs d'intervalles de temps de chaque frame
    C = reshape(C, [64,numel(C)/64])';                                      % Construction des frames en 64*64
    vals=length(C)-mod(length(C),64);

    data = cell(vals/64,1);                                                 %pre-allocation
    k = 0;
        for i= 1 : vals/64                                                  %nombre de capteurs
            data{i} = C(k+1:k+64,:);
            k = k+64;
        end
        
        data=cellfun(@(x) rot90(x,3),data,'uni',0);                         % Placement de la pressmat

%% Retirer le bruit systematique
    test=arrayfun(@(x) data{x}(end-3:end,:),1:length(data),'uni',0);
    maxArtefact=max(cellfun(@(x) sum(x,'all'),test));
    find(cellfun(@(x) sum(x,'all'),test)==maxArtefact);
    boolNoise=(data{1})==0;                                                	% On considère que sur la première frame, le pied n'est pas en contact: si tous les capteurs ne sont pas égaux à zéro c'est qu'il y a des artefacts
    boolNoise=~boolNoise;                                                 	% Booléen permettant d'avoir les capteurs qui ne sont pas égaux à zéro

        for iNoise = 1:numel(data)
            data{iNoise}(end-3:end,:) = 0;                                  % Je sais qu'il y a du bruit dans cette zone donc je rends toute cette zone égale à zéro
            data{iNoise}(1:3,:) = 0;                                    	% Je sais qu'il y a du bruit dans cette zone donc je rends toute cette zone égale à zéro
            data{iNoise}(:,1:3) = 0;                                       	% Je sais qu'il y a du bruit dans cette zone donc je rends toute cette zone égale à zéro
            data{iNoise}(:,end-3:end) = 0;                                 	% Je sais qu'il y a du bruit dans cette zone donc je rends toute cette zone égale à zéro
        end

%% Virer les cellules ou il ne se passe rien
    sommeCapteurs=cellfun(@(x) sum(x,'all'),data);                          % je fais la somme au sein des matrices 64*64: lorsqu'une matrice est vide, sa somme est égale à zéro.
    boolCellVide=sommeCapteurs==0;                                          % Booléen permettant de trouver les cellules vides
    data=data(~boolCellVide);                                               % Les cellules NON vides sont gardées.
% resume=[data, t_interp(~boolCellVide), num2cell(0:1/fs:(length(data)-1)/fs)'];


%% Visualisation du fichier au cours du temps
% for j=1:length(data_filt)
%     colormap("jet")
%     contourf(data_filt{j},50,'LineStyle','none','ShowText','off')
%     pause(0.1)
% end

dataRaw{kk}=data;
%% Calcul baricentre
    bari_xx=cell(1,length(data));                                               %pre Alloc
    bari_yy=cell(1,length(data));                                               %pre Alloc

        for i = 1:length(data)
            [xx,yy,zz]=find(data{i});
            bari_xx{i}=sum(xx.*zz,'all')/sum(data{i},'all');                        % Calcul baricentre x
            bari_yy{i}=sum(yy.*zz,'all')/sum(data{i},'all');                        % Calcul baricentre y
        end

% plot([bari_xx{:}],[bari_yy{:}],'r')
% hold on
% arrayfun(@(x,y) scatter(x,y,'r','filled'),cell2mat(bari_xx),cell2mat(bari_yy))


%% OPTION 2: Quelle est la frame avec le plus de capteurs allumés ?
    dataBinarize=cell(1,numel(data));
        for j = 1:numel(data)
            boolBinarize=data{j}==0;
            dataBinarize{j} = ~boolBinarize;                                        % Je sais qu'il y a du bruit dans cette zone donc je rends toute cette zone égale à zéro
        end
    sommeCapteursAllumes=cellfun(@(x) sum(x,'all'),dataBinarize);               % je fais la somme au sein des matrices 64*64: lorsqu'une matrice est vide, sa somme est égale à zéro.
    [~,PressionBinMax_Ind] = max(sommeCapteursAllumes);

% figure
% hold on
% plot([bari_xx{:}],[bari_yy{:}],'k')
% arrayfun(@(x,y) scatter(x,y,'k','filled'),cell2mat(bari_xx),cell2mat(bari_yy))
% contourf(rot90(data{PressionBinMax_Ind}',3),1000,'LineStyle','none')
% colormap("jet")
% exportgraphics(gca,"pressionsMat.png",Resolution=600)

%% Préparation à la rotation
    [rows,cols,Zz]=find(data{PressionBinMax_Ind});                              % J'arrete d'utiliser un matrice mais plutot des vecteurs de coordonées x et y
    k=[rows,cols];

% plot(rows, cols, '.','MarkerSize',20)
%% Analyse en Composante Principale (PCA)
% On recentre le pied autour de 0
    A=rows-mean(rows);                                                          % centrer en 0 sur abscises
    B=cols-mean(cols);                                                          % centrer sur 0 sur ordonnées
% On rotate de base le pied pour l'avoir dans le sens de la marche
    r = [cos(pi/2),-sin(pi/2);sin(pi/2),cos(pi/2)];                             % matrice de rotation pour mettre les pieds sur un plan antero-post
    rDataRAW = r*[A';B'];                                                       %on applique la rotation
    A=rDataRAW(1,:);
    B=rDataRAW(2,:);
%     figure
%     plot(A,B,'b.')                                                              %tracer le pied d'origine

%% PCA 1
    C=cov(A,B);                                                                 %covariance 
    [V,D]=eig(C);                                                               %eigen values
    alpha=atan(V(1)/V(2));                                                      %calcul de l'angle alpha

%Deux cas: angle positif ou négatif (selon orientation du pied)
    if alpha < 0
        alpha=alpha+pi/2;
    else
        alpha=alpha-pi/2;
    end

    r = [cos(alpha),-sin(alpha);sin(alpha),cos(alpha)];                         % matrice de rotation avec l'angle alpha
    rData = r*[A;B];                                                            % on applique la rotation

%% Visualisation graphique de PCA1
%     hold on
%     plot(A,B,'k.','MarkerSize',20); hold on                                     % plot raw data
%     plot(rData(1,:),rData(2,:),'r.','MarkerSize',20)                            % plot rotated data
%     axis([-44.2254   39.6416  -30.8524   28.6205])
% Affichage graphique des Eigen vectors 
    D_diag=diag(D);
%     plot([0 3*sqrt(D_diag(1))*V(1)],[0 3*sqrt(D_diag(1))*V(2)],'k','LineWidth',2);hold on   %PC1 RAW
%     plot([0 3*sqrt(D_diag(2))*V(3)],[0 3*sqrt(D_diag(2))*V(4)],'k','LineWidth',2)           %PC2 RAW
    alpha_deg=rad2deg(alpha);

%% PCA 2 pour tracer les nouveaux Eigen Vectors
    C=cov(rData(1,:),rData(2,:));                                               %covariance 
    [V,D]=eig(C);                                                               %eigen values
    alpha2=atan(V(2)/V(1));                                                     %calcul de l'angle alpha
    D_diag=diag(D);

%     plot([0 3*sqrt(D_diag(1))*V(1)],[0 3*sqrt(D_diag(1))*V(2)],'r','LineWidth',2);hold on   %PC1 ROTATED
%     plot([0 3*sqrt(D_diag(2))*V(3)],[0 3*sqrt(D_diag(2))*V(4)],'r','LineWidth',2)           %PC2 ROTATED

%% Maniulation et rotation du COP
% Replacer le COP en x=O et y=O
    bari_xCentre=cell2mat(bari_xx)-mean(rows);
    bari_yCentre=cell2mat(bari_yy)-mean(cols);
    % Rotate COP evolution
    r = [cos(pi/2+alpha),-sin(pi/2+alpha);sin(pi/2+alpha),cos(pi/2+alpha)];     % matrice de rotation avec l'angle alpha
    rDataRAW = r*[bari_xCentre;bari_yCentre];                                   %on applique la rotation
    bari_xCentreRotate=rDataRAW(1,:);
    bari_yCentreRotate=rDataRAW(2,:);
%     scatter(bari_xCentreRotate,bari_yCentreRotate,'filled');hold on

%% Créer un encadrement du pied pour determiner la largeur maximale et le COP Excursion.
% Ici l'objectif c'est de trouver le point le plus a droite et le plus à
% gauche sur le tier antérieur et sur le tier posterieur pour avoir largeur
% max vers les metas et le talon.

% ATTENTION: en fonction pied droit ou pied gauche different
% partie la plus externe vers le meta (a partir de 2/3 du pied)
    rData_x=rData(1,:); %x
    rData_y=rData(2,:); %y

% Trouver la largeur max à partir de 2/3 du pied
    thresholdHaut=min(rData_y)+(2/3)*(max(rData_y)-min(rData_y));               % emplacement du début du tier antérieur du pied
    YrDataThresholdedHaut=rData_y(rData_y>thresholdHaut);
    XrDataThresholdedHaut=rData_x(rData_y>thresholdHaut);
    indPosLargeur_Haut=find(floor(XrDataThresholdedHaut)==max(floor(XrDataThresholdedHaut)));
    indPosLargeur_HautFinal=indPosLargeur_Haut(round(median(1:numel(indPosLargeur_Haut))));
%     scatter(XrDataThresholdedHaut(indPosLargeur_HautFinal),YrDataThresholdedHaut(indPosLargeur_HautFinal),'filled')

    coord_extDroite_haute=[XrDataThresholdedHaut(indPosLargeur_HautFinal) YrDataThresholdedHaut(indPosLargeur_HautFinal)]; %x et y

% partie la plus externe vers le talon (<1/3 du pied)
    thresholdBas=min(rData_y)+(1/3)*(max(rData_y)-min(rData_y)); % emplacement du début du tier antérieur du pied
    YrDataThresholdedBas=rData_y(rData_y<thresholdBas);
    XrDataThresholdedBas=rData_x(rData_y<thresholdBas);
    indPosLargeur_basse=find(floor(XrDataThresholdedBas)==max(floor(XrDataThresholdedBas)));
    indPosLargeur_basseFinal=indPosLargeur_basse(round(median(1:numel(indPosLargeur_basse))));
%     scatter(XrDataThresholdedBas(indPosLargeur_basseFinal),YrDataThresholdedBas(indPosLargeur_basseFinal),'filled')


    coord_extDroite_basse=[XrDataThresholdedBas(indPosLargeur_basseFinal) YrDataThresholdedBas(indPosLargeur_basseFinal)]; %x et y
%     plot([coord_extDroite_basse(1) coord_extDroite_haute(1)],[coord_extDroite_basse(2) coord_extDroite_haute(2)],'--rs','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor',[0.5,0.5,0.5])

%% Calcul d'angle pour replacer la courbe qui passe par les deux extrémités perpendiculaire à l'axe des abscisses

    diffXY=coord_extDroite_basse-coord_extDroite_haute;                         %x2-x1 et y2-y1
    alpha_interieur=atan(diffXY(1)/diffXY(2));                                  %CAHSOHTOA: tan(alpha)=diffX/diffY
    rad2deg(alpha_interieur);                                                   %verif angle
% derniere rota pour aligner l'interieur des pieds
    r = [cos(alpha_interieur),-sin(alpha_interieur);sin(alpha_interieur),cos(alpha_interieur)]; % matrice de rotation avec l'angle alpha
    
    rData_New=r*rData;
%     plot(rData_New(1,:),rData_New(2,:),'b.','MarkerSize',20)                    % plot rotated data

% Je peux re-classer dans l'ordre les données x ou y avec les autres coordonées associées:
    [YSorted,I] = sort(rData_New(2,:),'ascend');
    XSorted=rData_New(1,I);

%rotate COP evolution
    rDataRAW2 = r*[bari_xCentreRotate;bari_yCentreRotate];                      %on applique la rotation
    bari_xCentreRotate2=rDataRAW2(1,:);
    bari_yCentreRotate2=rDataRAW2(2,:);
%     scatter(bari_xCentreRotate2,bari_yCentreRotate2,'filled');hold on



%% Début affichage graphique
    figure
    plot(XSorted,YSorted,'.','MarkerSize',30,'MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor',[0.7 0.7 0.7]); hold on
    axis([-38.6994   39.0647  -25.2727   22.0279])

%% Trouver les points extremes (4 points) automatiquement
%Trouver la valeur la plus extreme à gauche (5 eme meta)
    newX=XSorted(YSorted>0);
    newY=YSorted(YSorted>0);

    ind_minMeta=find(floor(newX)==min(floor(newX)));
    ind_inter_minMeta=ind_minMeta(newY(ind_minMeta)>0);
    medInd_minMeta=floor(median(1:length(ind_inter_minMeta)));

    minMeta_X=newX(ind_inter_minMeta(medInd_minMeta));
    minMeta_Y=newY(ind_inter_minMeta(medInd_minMeta));
    scatter(minMeta_X,minMeta_Y,50,'filled','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5])

% Trouver la valeur la plus extreme à droite (1er meta)
    newX=XSorted(YSorted>0);
    newY=YSorted(YSorted>0);

    ind_maxMeta=find(floor(newX)==max(floor(newX)));
    ind_inter_maxMeta=ind_maxMeta(newY(ind_maxMeta)>0);
    medInd_maxMeta=floor(median(1:length(ind_inter_maxMeta)));
    
    maxMeta_X=newX(ind_inter_maxMeta(medInd_maxMeta));
    maxMeta_Y=newY(ind_inter_maxMeta(medInd_maxMeta));
    scatter(maxMeta_X,maxMeta_Y,50,'filled','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5])

%trouver le point le plus extreme à droite mais en bas (talon)
    thresholdBasVisu=0;
    XSortedTheresholded=XSorted(YSorted<thresholdBasVisu);
    YSortedTheresholded=YSorted(YSorted<thresholdBasVisu);
    inddd=find(round(XSortedTheresholded)==max(round(XSortedTheresholded)));
    medInd_maxTalon=floor(median(1:length(inddd)));
%     inddd(medInd_maxTalon)
        
    maxTalon_X=XSortedTheresholded(inddd(medInd_maxTalon));
    maxTalon_Y=YSortedTheresholded(inddd(medInd_maxTalon));
    scatter(maxTalon_X,maxTalon_Y,50,'filled','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5])

% trouver le point le plus extreme à gauche mais en bas (talon)
    indPos_minTalon=YSorted<maxTalon_Y+0.8 & YSorted>maxTalon_Y-0.8;
    minTalon_Y=min(YSorted(indPos_minTalon));
    minTalon_X=min(XSorted(indPos_minTalon));
    scatter(minTalon_X,minTalon_Y,50,'filled','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5])
%% Visualisation graphique
    largeurX=[maxMeta_X;minMeta_X];
    largeurY=[maxMeta_Y;minMeta_Y];
    plot(largeurX,largeurY,'.-','MarkerSize',50,'LineWidth',4,'Color',"#EDB120")
% posLargeurY=a 2/3 de la hauteur:
    posLargeurY=min(YSorted)+2/3*(abs(max(YSorted)-min(YSorted)));
% plot(largeurX,[largeurY(2),largeurY(2)],'--','LineWidth',4)
    plot(largeurX,[posLargeurY,posLargeurY],'-.','LineWidth',2,'Color',"#A2142F")
    %% anticiper la pronation au talon ou les Morton's toes ET l'attaque mediopied INTERDITE
    % Point COP le plus medial vers les orteils (>2/3)  
    SeuilHautY=min(YSorted)+3/4*(abs(max(YSorted)-min(YSorted)));
    COPtopY=bari_yCentreRotate2(bari_yCentreRotate2>SeuilHautY);
    COPtopX=bari_xCentreRotate2(bari_yCentreRotate2>SeuilHautY);
    [~,indDernierCOP]=max(COPtopX);

    % Point COP le medial vers le talon (<1/3)
    % SeuilBasY a 1/3 de la hauteur:
    SeuilBasY=min(YSorted)+1/3*(abs(max(YSorted)-min(YSorted)));
    COPlowY=bari_yCentreRotate2(bari_yCentreRotate2<SeuilBasY);
    COPlowX=bari_xCentreRotate2(bari_yCentreRotate2<SeuilBasY);
    [~,indPremierCOP]=max(COPlowX);

    % Droite excursion line
%     yDroiteCOPDyna=[bari_yCentreRotate2([1, end])];
%     xDroiteCOPDyna=[bari_xCentreRotate2([1, end])];

    yDroiteCOPDyna=[COPlowY(indPremierCOP),COPtopY(indDernierCOP)];
    xDroiteCOPDyna=[COPlowX(indPremierCOP),COPtopX(indDernierCOP)];
    %% Visualisation graphique
% scatter(bari_xCentreRotate,bari_yCentreRotate,'filled');hold on
    plot(bari_xCentreRotate2,bari_yCentreRotate2,'.-','MarkerSize',30,'LineWidth',2,'Color',"#D95319") %orange
    plot(xDroiteCOPDyna,yDroiteCOPDyna,'.-','MarkerSize',30,'LineWidth',2,'Color',"#EDB120") %yellow
    footWidth=abs(diff(largeurX));

% La coordonée de la droite de la largeur c'est largeurY(2)
% Indice 1 c'est pour trouver le point le plus proche au dessus de ma droite
    ind_1=find(bari_yCentreRotate2>posLargeurY,1,'first');
    scatter(bari_xCentreRotate2(ind_1),bari_yCentreRotate2(ind_1),'filled','MarkerEdgeColor',"#EDB120",'MarkerFaceColor',"#EDB120");hold on
% Indice 2 c'est pour trouver le point le plus proche en dessous de ma droite
    ind_2=find(bari_yCentreRotate2<posLargeurY,1,'last');
    scatter(bari_xCentreRotate2(ind_2),bari_yCentreRotate2(ind_2),'filled','MarkerEdgeColor',"#EDB120",'MarkerFaceColor',"#EDB120");hold on

    %% Interpolation spline ?
% x = bari_xCentreRotate2
% y = bari_yCentreRotate2
% y_interp = min(y):0.1:max(y)
% x_interp_spli = interp1(y,x,y_interp,'spline');
% plot(x_interp_spli,y_interp,'-bo','LineWidth',1); hold on
% % plot(x,y,'ko','MarkerSize',20,'LineWidth',1)
% 

%% Coefficients de droite et intersections
    yDroiteCOP=[bari_yCentreRotate2(ind_2), bari_yCentreRotate2(ind_1)];
    xDroiteCOP=[bari_xCentreRotate2(ind_2), bari_xCentreRotate2(ind_1)];
% a=y2-y1/x2-x1
    coeffDroiteCOP=diff(yDroiteCOP)/diff(xDroiteCOP);
%y=ax+b -> b=y-ax -> y-(y2-y1/x2-x1)*x
    ordOrigineCOP=bari_yCentreRotate2(ind_1)-(coeffDroiteCOP*bari_xCentreRotate2(ind_1));
    x_interceptCOP=(ordOrigineCOP-posLargeurY)/(0-coeffDroiteCOP);

%% Intersection entre excursion et COP trajectory
% a=y2-y1/x2-x1
    coeffDroiteCOPDyna=diff(yDroiteCOPDyna)/diff(xDroiteCOPDyna);
%y=ax+b -> b=y-ax -> y-(y2-y1/x2-x1)*x
    ordOrigineCOPDyna=COPlowY(indPremierCOP)-(coeffDroiteCOPDyna*COPlowX(indPremierCOP));
    x_interceptCOPDyna=(ordOrigineCOPDyna-posLargeurY)/(0-coeffDroiteCOPDyna);


    plot([x_interceptCOP x_interceptCOPDyna],[posLargeurY posLargeurY],'MarkerSize',30,'LineWidth',4,'Color',[0.4350 0.0780 0.1840])

%% Center of Pressure Excursion Index
% Song et al. 1996
%x_interceptCOPDyna %courbe bleue
%x_interceptCOP %courbe verte
%footWidth %largeur globale
    CPE_Index=(x_interceptCOPDyna-x_interceptCOP)/footWidth;

%% Trouver la longueur du pied
%XSorted et YSorted: les coordonées du pied repositionné pour calcul CPEI

% trouver les coeff a et b de la droite distale du pied:
% a=y2-y1/x2-x1
%point 1: MinMeta_X et MinMeta_Y // point 2: minTalon_X et minTalon_Y
    coeffDroiteDistale=(minTalon_Y-minMeta_Y)/(minTalon_X-minMeta_X);
%trouver l'ordonée à l'origine
%y=ax+b -> b=y-ax -> y-(y2-y1/x2-x1)*x
%pour le point 1:
    ordOrigineDistal=minTalon_Y-(coeffDroiteDistale*minTalon_X);

%calculer l'interception en haut du pied
% Premiere interception= je connais deja y (distance entre min et max du pied)
    y_distZeroPiedMedialHAUT=max(YSorted);                                      %distance entre le point le plus haut et le 0
    x_interceptDistalHAUT=(y_distZeroPiedMedialHAUT-ordOrigineDistal)/coeffDroiteDistale;  % x=(y-b)/a
    scatter(x_interceptDistalHAUT,y_distZeroPiedMedialHAUT,'k','filled')

%Deuxième interception= je connais deja x (le point le plus à droite)
    x_interceptDistal=max(XSorted);
    y_interceptDistal=coeffDroiteDistale*x_interceptDistal+ordOrigineDistal;     %y=ax+b
    scatter(x_interceptDistal,y_interceptDistal,'k','filled')

    plot([x_interceptDistal x_interceptDistalHAUT],[y_interceptDistal y_distZeroPiedMedialHAUT],'k')

% Calcul angle entre les deux droites a gauche et droite du pied
% Trouver l'angle pour pouvoir tracer la bisectrice qui permettra de
% trouver la longueur du pied
%x2-x1 et y2-y1: distances
    diffY_distal=abs(y_distZeroPiedMedialHAUT-y_interceptDistal);
    diffX_distal=abs(x_interceptDistal-x_interceptDistalHAUT);
%CAHSOHTOA: tan(alpha)=diffX/diffY
    angleAlpha1=atan(diffX_distal/diffY_distal);
    rad2deg(angleAlpha1);

%Bisectrice
    angleAlpha2=angleAlpha1/2;

%tan(alpha2)=OP/ADJ // OP=tan(alpha)*ADJ
    distBisectrice=max(XSorted)-tan(angleAlpha2)*diffY_distal;
    scatter(distBisectrice,(max(YSorted)),'k','filled')

% line([minMeta_X minTalon_X],[minMeta_Y minTalon_Y])
    yline(min(YSorted))
    yline(max(YSorted))
    xline(max(XSorted))
    plot([distBisectrice max(XSorted)],[max(YSorted) y_interceptDistal],'k')

%% Esthetique graphique
    title(strcat(name," ","Essai ",num2str(cond)))
    subtitle(strcat('CPEI= ',num2str(CPE_Index*100,4),'\%'),'Interpreter','latex')
    hfig= gcf;                                                                  % save the figure handle in a variable
    ax=hfig.Children;
    picturewidth = 20;                                                          % set this parameter and keep it forever
    hw_ratio = 0.65;                                                            % feel free to play with this ratio
    set(findall(hfig,'-property','FontSize'),'FontSize',17)                     % adjust fontsize to your document
    
    set(findall(hfig,'-property','Box'),'Box','on')                             % box on ou off
    set(findall(hfig,'-property','TickLabelInterpreter'), ...                   % Lire les labels en LaTex
        'TickLabelInterpreter','latex')
    ax.Title.Interpreter="latex";
    ax.Subtitle.Interpreter="latex";
    set(hfig,'Units','centimeters','Position', ...                              % Placer la figure sur l'écran
        [3 3 picturewidth hw_ratio*picturewidth])
    pos = get(hfig,'Position');
    set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

%% CPEI et ANGLE a stocker

    CPEI(kk)=CPE_Index;
end





