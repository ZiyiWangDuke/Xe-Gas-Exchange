% function [dFID, gFID, dt, gt] = getFIDFromPfile(phaseCal_pfile, nDownstreamFrames,disToff,gasToff,doPrint)
% Force downstream frames to be shortest Te
nDownstreamFrames = ceil(nDownstreamFrames/4)*4+1;

%% Read Raw Pfile and process pfile
pfile = GE.Pfile.read(phaseCal_pfile);

% Remove baselines
pfile = MRI.DataProcessing.removeBaselineViews(pfile);

% Pull relavent info from header
npts = pfile.rdb.rdb_hdr_frame_size;
bw = pfile.rdb.rdb_hdr_user12;                 % Receiver bandwidth (kHz)
dwell_time = 1/(2*bw*1000);
dwell_time = Math.nearestMultipleOf(dwell_time,0.000002); % Dwell time must be an integer multible of 2us
pfile.rdb.rdb_hdr_user12 = 1/(2*dwell_time*1000);
bw = pfile.rdb.rdb_hdr_user12;
t = dwell_time*(0:(npts-1));%pfile.image.te*1E-6+dwell_time*(0:(npts-1)); %sec

% Print some stuff
if(doPrint)
    se_val = pfile.series.se_desc';
    bw_val = bw;
    tg_val = num2str(pfile.rdb.rdb_hdr_ps_mps_tg);
    r1_val = num2str(pfile.rdb.rdb_hdr_ps_mps_r1);
    r2_val = num2str(pfile.rdb.rdb_hdr_ps_mps_r2);
    freq_val = pfile.rdb.rdb_hdr_ps_mps_freq/10;
    patid_val = regexprep(strtrim(deblank(pfile.exam.patid')),'[^a-zA-Z0-9\-]','_');
    if(pfile.exam.ex_datetime == 0)
        %Create exam and series dates in YYYY_MM_DD format
        exam_timestamp = [pfile.rdb.rdb_hdr_scan_date(1:2)' '/' pfile.rdb.rdb_hdr_scan_date(4:5)' '20' pfile.rdb.rdb_hdr_scan_date(8:9)'        ];
    else
        date_number = pfile.exam.ex_datetime/86400 + datenum(1970,1,1);
        exam_timestamp = datestr(date_number,'mm/dd/yyyy');
    end
    disp(phaseCal_pfile)
    disp(['PatID: ' patid_val '    date:' exam_timestamp  '    se: ' se_val]);
    disp(['bw: ' num2str(bw_val) '   tg: ' tg_val '    r1: ' r1_val '    r2: ' r2_val '    freq: ' num2str(freq_val)]);
end

%% Split disolved data into separate pfile, being careful to handle
% different phase cal file formats/versions
if(strcmp(pfile.series.se_desc(1:5)','CALIB'))
    % Calculate average fid for TE 1
    disolvedPfileData = pfile.data(:,21:220);
%     disolvedPfileData = disolvedPfileData;
    
    gasPfileData = pfile.data(:,221);
elseif(strcmp(pfile.series.se_desc(1:5)','SPECT'))
    % Calculate average fid for TE
    disolvedPfileData = pfile.data(:,nDownstreamFrames:(end-1));
%     disolvedPfileData = disolvedPfileData;
    
    gasPfileData = pfile.data(:,201);
elseif(strcmp(pfile.series.se_desc(1:5)','PHASE'))
    % Calculate average fid for TE 1
    disolvedPfileData = pfile.data(:,nDownstreamFrames:4:200);
%     disolvedPfileData = disolvedPfileData;
    
    gasPfileData = pfile.data(:,201);
elseif(strcmp(pfile.series.se_desc(1:5)','NEW C'))   
    % Calculate average fid for TE 1
    disolvedPfileData = pfile.data(:,nDownstreamFrames:4:200);
%     disolvedPfileData = disolvedPfileData;
    
    gasPfileData = pfile.data(:,201);
else
    % Calculate average fid for TE 1
    disolvedPfileData = pfile.data(:,nDownstreamFrames:4:200);
%     disolvedPfileData = disolvedPfileData;
    
    gasPfileData = pfile.data(:,201);
end

dFID = disolvedPfileData((disToff+1):end,:);
dt = t((disToff+1):end);%-t((disToff+1));
gFID = gasPfileData((gasToff+1):end);
gt = t((gasToff+1):end);%-t((gasToff+1));