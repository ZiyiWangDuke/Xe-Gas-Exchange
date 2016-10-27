rda_file = filepath();

fid = fopen(rda_file );

head_start_text = '>>> Begin of header <<<';
head_end_text   = '>>> End of header <<<';


tline = fgets(fid);

while (isempty(strfind(tline , head_end_text)))
    
    tline = fgets(fid);
    
    if ( isempty(strfind (tline , head_start_text)) + isempty(strfind (tline , head_end_text )) == 2)
        
        
        % Store this data in the appropriate format
        
        occurence_of_colon = findstr(':',tline);
        variable = tline(1:occurence_of_colon-1) ;
        value    = tline(occurence_of_colon+1 : length(tline)) ;
        
        switch variable
        case { 'PatientID' , 'PatientName' , 'StudyDescription' , 'PatientBirthDate' , 'StudyDate' , 'StudyTime' , 'PatientAge' , 'SeriesDate' , ...
                    'SeriesTime' , 'SeriesDescription' , 'ProtocolName' , 'PatientPosition' , 'ModelName' , 'StationName' , 'InstitutionName' , ...
                    'DeviceSerialNumber', 'InstanceDate' , 'InstanceTime' , 'InstanceComments' , 'SequenceName' , 'SequenceDescription' , 'Nucleus' ,...
                    'TransmitCoil' }
            eval(['rda.' , variable , ' = value ']);
            
        case {  'SeriesNumber' , 'InstanceNumber' , 'AcquisitionNumber' , 'NumOfPhaseEncodingSteps' , 'NumberOfRows' , 'NumberOfColumns' , 'VectorSize' }
            %Integers
            eval(['rda.' , variable , ' = str2num(value) ']);
        case { 'PatientWeight' , 'TR' , 'TE' , 'TM' , 'DwellTime' , 'NumberOfAverages' , 'MRFrequency' , 'MagneticFieldStrength' , 'FlipAngle' , ...
                     'SliceThickness' ,  'FoVHeight' , 'FoVWidth' , 'PercentOfRectFoV' , 'PixelSpacingRow' , 'PixelSpacingCol'}
            %Floats 
            eval(['rda.' , variable , ' = str2num(value) ']);
        case {'CSIMatrixSize[0]' }
            rda.CSIMatrix_Size(1) = str2num(value);    
        case {'CSIMatrixSize[1]' }
            rda.CSIMatrix_Size(2) = str2num(value);    
        case {'CSIMatrixSize[2]' }
            rda.CSIMatrix_Size(3) = str2num(value);        
        otherwise
            % We don't know what this variable is.  Report this just to keep things clear
            %disp(['Unrecognised variable ' , variable ]);
        end
        
    else
        % Don't bother storing this bit of the output
    end
    
end

% Read FID data
bytes_per_point = 16;
fid_data = fread(fid , rda.CSIMatrix_Size(1) * rda.CSIMatrix_Size(1) *rda.CSIMatrix_Size(1) *rda.VectorSize * 2 , 'double');  
fclose(fid);

% Make FID data complex
fid_data = reshape(fid_data,  2 , rda.VectorSize , rda.CSIMatrix_Size(1) ,  rda.CSIMatrix_Size(2) ,  rda.CSIMatrix_Size(3) );
fid_data = complex(fid_data(1,:,:,:,:),fid_data(2,:,:,:,:));
t = (0:(rda.VectorSize -1))*rda.DwellTime*1E-6;
%%
% Fit spectrum
chopoff = 5;
fid_data = fid_data(chopoff:end);
t = t(chopoff:end);
% fitObj = NMR_TimeFit(fid_data,t,[1;1;1;1;1],[-5;-4;-3;-2;-1],[30;30;30;30;30],[0;0;0;0;0],[],[]); 
fitObj = NMR_TimeFit(fid_data,t,[1],[1],[30],[0],[],[]); 
% fitObj = NMR_TimeFit(fid_data,t,[1;2],[1;0],[30;29],[0;0],[],[]); 
%% %%% for noise study
noise = abs(fitObj.timeDomainSignal);
% startindex = round(length(noise)*0.75);
startindex =1;
stdnoise = std(noise(startindex:end));
%%
stdnoise = std(abs(fitObj.timeDomainSignal));
figure;
bw = 10000;
data = abs(fitObj.timeDomainSignal);
dwell = 1/bw;
T = length(data)*dwell;
t = dwell:dwell:T;
dwellf = 1/T;
f = dwellf:dwellf:bw;
subplot(2,1,1);plot(t,abs(fitObj.timeDomainSignal));ylim([0,1500]);xlim([0,0.2044])
subplot(2,1,2);plot(f,abs(fft(fitObj.timeDomainSignal)));ylim([0,80000])
std(abs(fft(fitObj.timeDomainSignal)))
%%
fitObj = fitObj.fitTimeDomainSignal();
fitObj.plotTimeAndSpectralFit;
fitObj.area
%% calculate SNR
dwell_time = (fitObj.t(2)-fitObj.t(1));
zeroPaddedTime = min(fitObj.t(:)) + dwell_time*((1:fitObj.zeroPadSize)-1)';
signal = fitObj.calcTimeDomainSignal(zeroPaddedTime);
fitNoise = signal - fitObj.timeDomainSignal;
% take off the first 25% of data points
start_index = round(length(fitNoise)*0.75);
snr_estimate_dB = sum(fitObj.area())/std(abs(fitNoise(start_index:end)));
% snr_estimate_dB = (fitObj.area(1))/std(abs(fitNoise(start_index:end)));
% a = std(abs(fitNoise(start_index:end)))