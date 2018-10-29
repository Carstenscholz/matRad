clear
close all
clc


%load HEAD_AND_NECK
load TG119.mat
%load PROSTATE.mat
%load LIVER.mat
%load BOXPHANTOM.mat

% meta information for treatment plan
pln.numOfFractions  = 30;
pln.radiationMode   = 'carbon';           % either photons / protons / helium / carbon
pln.machine         = 'Generic';

% beam geometry settings
pln.propStf.bixelWidth      = 5; 
pln.propStf.gantryAngles    = [0];
pln.propStf.couchAngles     = [0];
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln.propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

quantityOpt  = 'RBExD';     % options: physicalDose, effect, RBExD
modelName    = 'LEM';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

pln.multScen = matRad_multScen(ct,scenGenType);

stf = matRad_generateStf(ct,cst,pln);

cst{2,5}.alphaX = 0.1;

dij = matRad_calcParticleDose(ct,stf,pln,cst);

resultGUI1  = matRad_fluenceOptimization(dij,cst,pln);

%% 

cst{2,5}.alphaX = 0.5;
dij = matRad_calcParticleDose(ct,stf,pln,cst);
effect2vec  = dij.mAlphaDose{1}*resultGUI1.w + (dij.mSqrtBetaDose{1}*resultGUI1.w).^2;
effect2 = reshape(effect2vec, ct.cubeDim);

%% 

resultGUI2 = matRad_calcCubes(resultGUI1.w,dij,cst)

%%
dvh_lowalphaEffect = matRad_calcDVH(cst,resultGUI1.effect);
dvh_highalphaEffect = matRad_calcDVH(cst,effect2);

dvh_lowalphaRBExD = matRad_calcDVH(cst,resultGUI1.RBExD);
dvh_highalphaRBExD = matRad_calcDVH(cst,resultGUI2.RBExD);

figure
matRad_showDVH(dvh_lowalphaRBExD,cst,pln,1);
hold on
matRad_showDVH(dvh_highalphaRBExD,cst,pln,2);
title('straight line alphaX/betaX = 2, dotted line alphaX/betaX = 10 (target)')
hold off

figure
matRad_showDVH(dvh_lowalphaEffect,cst,pln,1);
hold on
matRad_showDVH(dvh_highalphaEffect,cst,pln,2);
xlabel('Effect [a.u.]');
title('straight line alphaX/betaX = 2, dotted line alphaX/betaX = 10 (target)')
hold off

% plane = 3;
% slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
%doseWindow = [0 max([resultGUI1.RBExD(:); resultGUI2.RBExD(:)])];
%doseIsoLevels = [0:0.1:1.2]*max([resultGUI1.RBExD(:); resultGUI2.RBExD(:)]);
% effectWindow = [0 max([resultGUI1.effect(:); effect2(:)])];
% effectIsoLevels = [0:0.1:1.2]*max([resultGUI1.effect(:); effect2(:)]);

% figure
% matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI1.RBExD,plane,slice,[],0.75,colorcube,[],doseWindow,doseIsoLevels );
% title('alphaX/betaX = 2 (target) RBExD');
% 
% figure
% matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI2.RBExD,plane,slice,[],0.75,colorcube,[],doseWindow,doseIsoLevels );
% title('alphaX/betaX = 10 (target)RBExD');

% figure
% matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI1.effect,plane,slice,[],0.75,colorcube,[],effectWindow,effectIsoLevels);
% title('alphaX/betaX = 2 (target) effect');
% 
% figure
% matRad_plotSliceWrapper(gca,ct,cst,1,effect2,plane,slice,[],0.75,colorcube,[],effectWindow,effectIsoLevels);
% title('alphaX/betaX = 10 (target)effect');

