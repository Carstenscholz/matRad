function [PatStruct] = dicomLoop(Folder)

clear PatStruct;

%% Loop for approved RT plans
%'C:\Users\c129z\Documents\CLEO_Export\'
folderList = dir(Folder);

for i = 3:numel(folderList)                     %Loop trough patients (start by 3 because of hidden folders)
     
    PatientSeries = dir([folderList(i).folder '\' folderList(i).name]);
    
    PatStruct(i-2).PatientName = folderList(i).name;                                                      %Adding patient name
    
    for j = 3:numel(PatientSeries)              %Loop through different therapy series (start by 3 because of hidden folders)
    
        PatientRTPlans = dir([PatientSeries(j).folder '\' PatientSeries(j).name '\RTPLAN']);

        RTPlanNumber = 1;

        for k = 3:numel(PatientRTPlans)             %Loop through RTPlans (start by 3 because of hidden folders)

            dcmInfRTPlan = dicominfo([PatientRTPlans(k).folder '\' PatientRTPlans(k).name]);

            if strcmp(dcmInfRTPlan.ApprovalStatus,'APPROVED')        %Looking for approved RTPlans

                PatStruct(i-2).Patient(j-2).Series(RTPlanNumber).RTPlan.folderlocation = ([PatientRTPlans(k).folder, '\', PatientRTPlans(k).name]);       %Adding RTPlan folderlocation
              
                PatStruct(i-2).Patient(j-2).Series(RTPlanNumber).Info.Patient = folderList(i).name;                                                      %Adding patient name
              
                PatStruct(i-2).Patient(j-2).Series(RTPlanNumber).Info.PatientSeries = [PatientSeries(j).folder '\' PatientSeries(j).name];                %Adding patient series folder location


                if isfield(dcmInfRTPlan,{'IonBeamSequence'})           %Looking for info about ion plans       
                      
                    PatStruct(i-2).Patient(j-2).Series(RTPlanNumber).Info.Particle = dcmInfRTPlan.PrescriptionDescription;        %Adding info about Particle and fractions
                      
                    PatStruct(i-2).Patient(j-2).Series(RTPlanNumber).Info.Information = dcmInfRTPlan.IonBeamSequence.Item_1.TreatmentMachineName;        %Adding infos
                  
                    if isfield(dcmInfRTPlan.IonBeamSequence.Item_1,{'NumberOfRangeShifters'})
                      
                        PatStruct(i-2).Patient(j-2).Series(RTPlanNumber).Info.NumOfRangeShifters = dcmInfRTPlan.IonBeamSequence.Item_1.NumberOfRangeShifters;
                  
                    end    
                      
                elseif isfield(dcmInfRTPlan,{'BeamSequence'})          %Looking for info about other plans (mostly photon)
                  
                    PatStruct(i-2).Patient(j-2).Series(RTPlanNumber).Info.Particle = dcmInfRTPlan.BeamSequence.Item_1.RadiationType;        %Adding info about Particle  
                  
                    PatStruct(i-2).Patient(j-2).Series(RTPlanNumber).Info.Room = dcmInfRTPlan.BeamSequence.Item_1.TreatmentMachineName;        %Adding info about treatment room/photon accelerator
              
                end
              
                PatStruct(i-2).Patient(j-2).Series(RTPlanNumber).Info.NumOfBeams = dcmInfRTPlan.FractionGroupSequence.Item_1.NumberOfBeams;   %Adding info about beam number
              
                if isfield(dcmInfRTPlan,{'ReferencedDoseSequence'})       %Looking for referenced doses in the RTPlan (in that form only ion doses)
                                                                         
                    for l = 1:numel(fieldnames(dcmInfRTPlan.ReferencedDoseSequence))      %Loop through referenced doses

                        refDoseNamesIon = fieldnames(dcmInfRTPlan.ReferencedDoseSequence);

                        PatStruct(i-2).Patient(j-2).Series(RTPlanNumber).Dose.Dosage(l).Dose = [PatientSeries(j).folder '\' PatientSeries(j).name '\RTDOSE\' dcmInfRTPlan.ReferencedDoseSequence.(refDoseNamesIon{l}).ReferencedSOPInstanceUID '.dcm'];
                        %Adding referenced doses

                    end
                  
                else
                  
                    DoseNumber = 1;
                  
                    for m = 3:numel(dir([PatientSeries(j).folder '\' PatientSeries(j).name '\RTDOSE']))       %Loop through all doses
                      
                        RTDoses = dir([PatientSeries(j).folder '\' PatientSeries(j).name '\RTDOSE']);
                      
                        dcmInfRTDoses = dicominfo([PatientSeries(j).folder '\' PatientSeries(j).name '\RTDOSE\' RTDoses(m).name]);
                                           
                      
                        if strcmp([dcmInfRTDoses.ReferencedRTPlanSequence.Item_1.ReferencedSOPInstanceUID '.dcm'], PatientRTPlans(k).name)        %Looking if the photon RTPlan is referenced
                          
                            PatStruct(i-2).Patient(j-2).Series(RTPlanNumber).Dose.Dosage(DoseNumber).Dose = [PatientSeries(j).folder '\' PatientSeries(j).name '\RTDOSE\' RTDoses(m).name];          %Adding referenced doses
                            DoseNumber = DoseNumber + 1;
                                                   
                        end                       
                    end
                end

                PatStruct(i-2).Patient(j-2).Series(RTPlanNumber).RTStruct.Struct = [PatientSeries(j).folder '\' PatientSeries(j).name '\RTSTRUCT\' dcmInfRTPlan.ReferencedStructureSetSequence.Item_1.ReferencedSOPInstanceUID '.dcm'];
                %Adding referenced RTStruct

                check = dir([PatientSeries(j).folder '\' PatientSeries(j).name '\RTSTRUCT']);         %Used for checking if a CT series is linked to the RTStruct

                for n = 3:numel(check)        %Loop through all existing RTStruct

                    if strcmp([check(n).folder '\' check(n).name],[PatStruct(i-2).Patient(j-2).Series(RTPlanNumber).RTStruct.Struct])      %Checks if the referenced RTStruct is in the existing RTStructs

                        RTStructFolderLoc = PatStruct(i-2).Patient(j-2).Series(RTPlanNumber).RTStruct.Struct;
                        dcmInfRTStruct = dicominfo(RTStructFolderLoc);
                        PatStruct(i-2).Patient(j-2).Series(RTPlanNumber).CTSeries.CT = dcmInfRTStruct.(dicomlookup('3006','0010')).Item_1.(dicomlookup('3006','0012')).Item_1.(dicomlookup('3006','0014')).Item_1.(dicomlookup('0020','000E'));     
                        %Adding referenced CTSeries
                    end

                end

                RTPlanNumber = RTPlanNumber + 1;

            end
        end
    end
end
  
clearvars -except PatStruct
end







