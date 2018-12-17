xResolution = 3;

yResolution = 3;

zResolution = 3;

dicomMetaBool = 1;

BoolUseDoseGrid = 0;


%% 

for i = 78:78     %1:numel(PatStruct)                                                      %Loop through patients (i patient index)
       
        
        if isempty(PatStruct(i).Patient(1).Series)                  
            PatientFolder = ['C:\Users\c129z\Documents\Files\' PatStruct(i).Patient(2).Series(1).Info.Patient];           
        else
            PatientFolder = ['C:\Users\c129z\Documents\Files\' PatStruct(i).Patient(1).Series(1).Info.Patient];
        end
    
        mkdir(PatientFolder);
    
        for j = 1:1 %1:numel(PatStruct(i).Patient)                                  %Loop trhough patient series (j patient series index)
        
            SeriesFolder = [PatientFolder '\Series_' int2str(j)];
        
            mkdir(SeriesFolder);
    
            for k = 1:1 %1:numel(PatStruct(i).Patient(j).Series)                           %Loop through RTPlans (k RTPlan index)
                
                %try
                    
    
                    if ((isempty(PatStruct(i).Patient(j).Series(k).CTSeries) == 0) && (isempty(PatStruct(i).Patient(j).Series(k).Dose) == 0))           %Check to not run into the RTPlan without referenced CTSeries and Dose
                    
                        clear files CTinfo dcmInfCT ctLocation;
    
                        files.rtplan{1,1} = PatStruct(i).Patient(j).Series(k).RTPlan.folderlocation;

                        files.rtss{1,1} = PatStruct(i).Patient(j).Series(k).RTStruct.Struct;
     
                        for l = 1:numel(PatStruct(i).Patient(j).Series(k).Dose.Dosage)

                            files.rtdose{l,1} = PatStruct(i).Patient(j).Series(k).Dose.Dosage(l).Dose;
                            files.rtdose{l,2} = 'RTDOSE';
    
                        end

                        files.resx = xResolution;

                        files.resy = yResolution;

                        files.resz = zResolution;
   
                        CTinfo = dir([PatStruct(i).Patient(j).Series(k).Info.PatientSeries '\CT']);

                        CTnumber = 1;

                        for m = 3:numel(CTinfo)
    
                            dcmInfCT = dicominfo([CTinfo(m).folder '\' CTinfo(m).name]);
    
                            if strcmp(dcmInfCT.SeriesInstanceUID,PatStruct(i).Patient(j).Series(k).CTSeries.CT)
    
                                ctLocation.folder{CTnumber,1} = [CTinfo(m).folder '\' CTinfo(m).name];
        
                                CTnumber = CTnumber + 1;
        
                            end
    
                        end
            %% 
    
                        files.ct = ctLocation.folder;
                   
                        files.useDoseGrid = BoolUseDoseGrid;
                
                        files.check = 1;
            
                        if isfield(PatStruct(i).Patient(j).Series(k).Info,{'Information'})

                            files.room = PatStruct(i).Patient(j).Series(k).Info.Information;
                
                        end
               
            
                        [ct, cst, pln, resultGUI,stf] = matRad_importDicom(files, dicomMetaBool);
                               
                
                        FileName = ['RTPlan_' int2str(k) '_' pln.radiationMode '.mat'];
                        PathName = [SeriesFolder '\'];
            
                        save([PathName, FileName],'ct', 'cst', 'pln', 'resultGUI','stf','-v7.3')
                    end
                %end
            end
        end  
end









