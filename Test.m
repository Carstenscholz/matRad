load TG119.mat

pln.numOfFractions  = 30;
pln.radiationMode   = 'photons';           % either photons / protons / helium / carbon
pln.machine         = 'Generic';


pln.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles    = [0:72:359]; % [?] ;
pln.propStf.couchAngles     = [0 0 0 0 0]; % [?] ; 
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln.propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,scenGenType);


stf = matRad_generateStf(ct,cst,pln);

if strcmp(pln.radiationMode,'photons')
    dij = matRad_calcPhotonDose(ct,stf,pln,cst);
    %dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'helium') || strcmp(pln.radiationMode,'carbon')
    dij = matRad_calcParticleDose(ct,stf,pln,cst);
end

resultGUI  = matRad_fluenceOptimization(dij,cst,pln);

[dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUI);
%% plot slice wrapper 
 colorAssigned =true;
for i = 1:size(cst,1)
    if ~isfield(cst{i,5},'visibleColor')
        colorAssigned = false;
        break;
    elseif isempty(cst{i,5}.visibleColor)
        colorAssigned = false;
        break;
    end
end

% assign color if color assignment is not already present or inconsistent
if colorAssigned == false
  m         = 64;
  colorStep = ceil(m/size(cst,1));
  colors    = colorcube(colorStep*size(cst,1));
  % spread individual VOI colors in the colorcube color palette
  colors    = colors(1:colorStep:end,:);
  for i = 1:size(cst,1)
    cst{i,5}.visibleColor = colors(i,:);
  end
end

 
slice = round((pln.propStf.isoCenter(1,3)-2)/ct.resolution.z);
c = jet(100);
doseIsoLevels = [0:0.1:1.2]*max(resultGUI.physicalDose(:));
figure;

[hCMap,hDose,hCt,hContour,hIsoDose] =...
    matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDose,3,slice,[],[],[],c,[],doseIsoLevels);


