%% head
clear
close all
clc

%load ct file
load TG119.mat


%personal settings
pln.propStf.gantryAngles_Photon = [0:72:359];
pln.propStf.gantryAngles_Proton = [0 180];
pln.propStf.gantryAngles_Carbon = [0];

sliceSelection = 0;                              % 0 == 162.5 mm; <0 positive z direction

alphaX = 0.1;
betaX = 0.05;

showPhysicalDose = false;                       %shows physikal dose
showEffect = false;                             %shows effect
showMixedPhysicalDose = true;                  %shows mixed plans and compares them
showMixedEffect = true;                         %shows mixed plans effect compares them

fractionDivider = 1/3;                          % percentage of used particle fractions

%general settings
pln.numOfFractions  = 30;
pln.propStf.bixelWidth      = 5;

pln.propOpt.runDAO          = false;   
pln.propOpt.runSequencing   = false; 

scenGenType  = 'nomScen';
pln.multScen = matRad_multScen(ct,scenGenType); 

%% Photon part
pln.radiationMode   = 'photons';          
pln.machine         = 'Generic';

pln.propStf.gantryAngles    = pln.propStf.gantryAngles_Photon;
pln.propStf.couchAngles     = zeros(1,numel(pln.propStf.gantryAngles)); 
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

quantityOpt  = 'physicalDose';    
modelName    = 'none'; 

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

stf = matRad_generateStf(ct,cst,pln);

dij = matRad_calcPhotonDose(ct,stf,pln,cst);

resultGUI_Photon  = matRad_fluenceOptimization(dij,cst,pln);

[dvh_Photon,qi_Photon] = matRad_indicatorWrapper(cst,pln,resultGUI_Photon);
title('Photon DVH')

%% Proton part
pln.radiationMode   = 'protons';           
pln.machine         = 'Generic';

pln.propStf.bixelWidth      = 5;
pln.propStf.gantryAngles    = pln.propStf.gantryAngles_Proton;
pln.propStf.couchAngles     = zeros(1,numel(pln.propStf.gantryAngles)); 
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

quantityOpt = 'RBExD';
modelName = 'constRBE';

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

stf = matRad_generateStf(ct,cst,pln);

dij = matRad_calcParticleDose(ct,stf,pln,cst);

resultGUI_Proton  = matRad_fluenceOptimization(dij,cst,pln);

[dvh_Proton,qi_Proton] = matRad_indicatorWrapper(cst,pln,resultGUI_Proton);
title('Proton DVH');


%% Carbon
pln.radiationMode   = 'carbon';           
pln.machine         = 'Generic';

pln.propStf.bixelWidth      = 5;
pln.propStf.gantryAngles    = pln.propStf.gantryAngles_Carbon;
pln.propStf.couchAngles     = zeros(1,numel(pln.propStf.gantryAngles)); 
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

quantityOpt = 'RBExD';
modelName = 'LEM';

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

stf = matRad_generateStf(ct,cst,pln);

dij = matRad_calcParticleDose(ct,stf,pln,cst);

resultGUI_Carbon  = matRad_fluenceOptimization(dij,cst,pln);

[dvh_Carbon,qi_Carbon] = matRad_indicatorWrapper(cst,pln,resultGUI_Carbon);
title('Carbon DVH');

%% Comparison

%DVH in one picture
figure
matRad_showDVH(dvh_Photon,cst,pln,1);
hold on
matRad_showDVH(dvh_Proton,cst,pln,2);
matRad_showDVH(dvh_Carbon,cst,pln,3);
hold off
title('straight line Photon, dotted line Proton, dashed line Carbon');


%slices
plane = 3;
slice = round((pln.propStf.isoCenter(1,3) - sliceSelection)./ct.resolution.z);


effect_Photon = alphaX .* resultGUI_Photon.physicalDose + betaX .* (resultGUI_Photon.physicalDose.^2);
effect_Proton = alphaX .* resultGUI_Proton.physicalDose * 1.1 + betaX .* ((resultGUI_Proton.physicalDose.^2 )*1.1^2);
effect_Carbon = resultGUI_Carbon.effect;

doseWindow = [0 max([resultGUI_Photon.physicalDose(:); resultGUI_Proton.RBExD(:); resultGUI_Carbon.RBExD(:)])];
effectWindow = [0 max([effect_Photon(:); effect_Proton(:);effect_Carbon(:)])];

doseIsoLevels = [0:0.1:1.2]*max([resultGUI_Photon.physicalDose(:); resultGUI_Proton.physicalDose(:); resultGUI_Carbon.RBExD(:)]);
effectIsoLevels = [0:0.1:1.2]*max([effect_Photon(:); effect_Proton(:);effect_Carbon(:)]);

% absDiffCube_physicalDose = resultGUI_Photon.physicalDose-resultGUI_Proton.physicalDose;
% absDiffCube_effect = effect_Photon-effect_Proton;

if showPhysicalDose == true
    figure
    matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_Photon.physicalDose,plane,slice,[],0.75,colorcube,[],doseWindow,doseIsoLevels );
    title('Photon physical dose');

    figure
    matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_Proton.RBExD,plane,slice,[],0.75,colorcube,[],doseWindow,doseIsoLevels );
    title('Proton physical dose');
    
    figure
    matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_Carbon.RBExD,plane,slice,[],0.75,colorcube,[],doseWindow,doseIsoLevels );
    title('Carbon RBExD');
    
%     figure
%     matRad_plotSliceWrapper(gca,ct,cst,1,absDiffCube_physicalDose,plane,slice,[],[],colorcube);
%     title('Photon dose - Proton dose');
end


if showEffect == true
    figure
    matRad_plotSliceWrapper(gca,ct,cst,1,effect_Photon,plane,slice,[],0.75,colorcube,[],effectWindow,effectIsoLevels);
    title('Photon effect');

    figure
    matRad_plotSliceWrapper(gca,ct,cst,1,effect_Proton,plane,slice,[],0.75,colorcube,[],effectWindow,effectIsoLevels);
    title('Proton effect');

    figure
    matRad_plotSliceWrapper(gca,ct,cst,1,effect_Carbon,plane,slice,[],0.75,colorcube,[],effectWindow,effectIsoLevels);
    title('Carbon effect');
  

%     figure
%     matRad_plotSliceWrapper(gca,ct,cst,1,absDiffCube_effect,plane,slice,[],[],colorcube);
%     title('Photon effect - Proton effect');
end

%% Dose
Dose = [];
j = 0;
for i = 1:size(cst,1)
    if (isempty(cst{i,6}) ~= 1)
        j = j + 1;
        Dose{j,1} = cst{i,2};
        Dose{j,2} = cst{i,6}.dose;
        Dose{j,3} = (alphaX  * cst{i,6}.dose + betaX * cst{i,6}.dose^2/pln.numOfFractions);
        Dose{j,4} = (alphaX * 1.1 * cst{i,6}.dose + betaX * (1.1 * cst{i,6}.dose)^2/pln.numOfFractions);    
    end
end

%% Adding - Mixed Treatment Plans

%absAddCube_physicalDose = resultGUI_Photon.physicalDose + resultGUI_Proton.physicalDose;

addedPhysicalDose_ProtonPhoton = pln.numOfFractions * (resultGUI_Proton.RBExD * fractionDivider + resultGUI_Photon.physicalDose * (1 - fractionDivider));
addedPhysicalDose_CarbonPhoton = pln.numOfFractions * (resultGUI_Carbon.RBExD * fractionDivider + resultGUI_Photon.physicalDose * (1 - fractionDivider));

addedEffect_ProtonPhoton = pln.numOfFractions * (effect_Proton * fractionDivider + effect_Photon * (1 - fractionDivider));
addedEffect_CarbonPhoton = pln.numOfFractions * (effect_Carbon * fractionDivider + effect_Photon * (1 - fractionDivider));

addedDoseIsoLevels = [0:0.1:1.2]*max([addedPhysicalDose_ProtonPhoton(:);addedPhysicalDose_CarbonPhoton(:)]);
addedEffectIsoLevels = [0:0.1:1.2]*max([addedEffect_ProtonPhoton(:); addedEffect_CarbonPhoton(:)]);

addedDoseWindow = [0 max([addedPhysicalDose_ProtonPhoton(:); addedPhysicalDose_CarbonPhoton(:)])];
addedEffectWindow = [0 max([addedEffect_ProtonPhoton(:); addedEffect_CarbonPhoton(:)])];

if showMixedPhysicalDose == true   

%     figure
%     matRad_plotSliceWrapper(gca,ct,cst,1,addedPhysicalDose_ProtonPhoton,plane,slice,[],[],colorcube,[],addedDoseWindow,addedDoseIsoLevels);
%     title('Mixed Plan Proton Photon Physical Dose');
% 
%     figure
%     matRad_plotSliceWrapper(gca,ct,cst,1,addedPhysicalDose_CarbonPhoton,plane,slice,[],[],colorcube,[],addedDoseWindow,addedDoseIsoLevels);
%     title('Mixed Plan Carbon Photon Physical Dose');
   
    dvh_ProtonPhotonDose = matRad_calcDVH(cst,addedPhysicalDose_ProtonPhoton);
    dvh_CarbonPhotonDose = matRad_calcDVH(cst,addedPhysicalDose_CarbonPhoton);
    
    figure
    matRad_showDVH(dvh_ProtonPhotonDose,cst,pln,1);
    hold on
    matRad_showDVH(dvh_CarbonPhotonDose,cst,pln,2);
    hold off
    title('Mixed Plan DVH for Physical Dose, straight line Proton Photon, dotted line Carbon Photon');
%     figure
%     matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_Photon.physicalDose*pln.numOfFractions,plane,slice,[],0.75,colorcube,[],doseWindow_physicalDose*pln.numOfFractions,[]);
%     title('Photon physical dose');
%     
%     figure
%     matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_Proton.physicalDose*pln.numOfFractions,plane,slice,[],0.75,colorcube,[],doseWindow_physicalDose*pln.numOfFractions,[]);
%     title('Proton physical dose');

end  
 
if showMixedEffect == true

%    figure
%    matRad_plotSliceWrapper(gca,ct,cst,1,addedEffect_ProtonPhoton,plane,slice,[],[],colorcube,[],addedEffectWindow,addedEffectIsoLevels);
%    title('Mixed Plan Proton Photon Effect');
% 
%    figure
%    matRad_plotSliceWrapper(gca,ct,cst,1,addedEffect_CarbonPhoton,plane,slice,[],[],colorcube,[],addedEffectWindow,addedEffectIsoLevels);
%    title('Mixed Plan Carbon Photon Effect');
   
   dvh_ProtonPhotonEffect = matRad_calcDVH(cst,addedEffect_ProtonPhoton);
   dvh_CarbonPhotonEffect = matRad_calcDVH(cst,addedEffect_CarbonPhoton);
    
   figure
   matRad_showDVH(dvh_ProtonPhotonEffect,cst,pln,1);
   hold on
   matRad_showDVH(dvh_CarbonPhotonEffect,cst,pln,2);
   hold off
   title('Mixed Plan DVH for Effect, straight line Proton Photon, dotted line Carbon Photon');
   xlabel('Effect [a.u.]');
%     figure
%     matRad_plotSliceWrapper(gca,ct,cst,1,effect_Photon*pln.numOfFractions,plane,slice,[],0.75,colorcube,[],doseWindow_effect*pln.numOfFractions,[]);
%     title('Photon effect');
% 
%     figure
%     matRad_plotSliceWrapper(gca,ct,cst,1,effect_Proton*pln.numOfFractions,plane,slice,[],0.75,colorcube,[],doseWindow_effect*pln.numOfFractions,[]);
%     title('Proton effect');

end



