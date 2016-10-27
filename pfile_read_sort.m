function mag_data=pfile_read_sort(p_number,view_points,loopfactor,ge_rev,flip_cal);
% This just reads in a p-file of the .split type to display the raw data
% Supply the Pfile number (no 'p' no '.7') residing in the same directory
% differs from pfile_read by displaying real and imag_datainary too
% calculates flip angle from first two fids
% updated 12/1/07 - to undo loopfactor from acquired data
% updated 1/12/08 - for better upclose display, handling file names
% updated 3/3/08 - to use common filename routine and raw read routine.
% updated 3/3/08 - now also reads in direct p-files from scanner
% updated 7/11/08 - remove baseline junk views from direct p-files
% updated 10/31/08 - to do flip angle calibration
% updated 6/19/09 - export mag_datanitude data
% updated 9/1//09 - updated to remove every 2048th frame
% updated 9/1/09 - require passing ge_rev and flip_cal, not nframes
% updated 2/17/10 - deal with files with num_fids<20
% updated 11/29/10 - deal with extra junk view in 15X
% updated 12/1/10 - fixed bug, remove every 2049th rather than 2048th frame
% updated 12/14/10 - Standard display of various magnitude views
% updated 3/14/11 - remove baselines on all data to preserve phase info
% updated 5/23/11 - streamline peak finding in flip angle calibration
% sample function call
% mag_data=pfile_read_sort(09728,256,271,12,0);

% NOTES REGARDING LOOPFACTORS
% common loop factors for 400,800,1600 frames are 131, 271, 661
% Other combos nframes/loop 11520/6733, 5120/2999, 3751/2287, 2400/1063
% note that for 3dpr npts=88 for opxres=128. Same for 2dpr?

%% control variables
plot_baselines=0;                % use to plot baseline views before removing
plot_phases=1;                   % plot real imaginary and magnitude data (useful for spectra)
first_frames=5;                    % number of initial fids to plot
breath_frames=20;                  % number of fids to plot to see breath dynamics nvptrig

a=sprintf('\n');disp(a);         % get first carriage return done

% determine header size based on ge_rev 11,12,14,15
ge_rev=floor(ge_rev);             % in case of 14.5, etc
if ge_rev==11
    header_size=61464                % EXCITE 11.0
elseif ge_rev==12
    header_size=66072;                 % EXCITE 12.0
elseif ge_rev==14
    header_size = 145908;              % EXCITE 14.0
elseif ge_rev==15
    header_size = 145908;              % EXCITE 15.0
    add_junk = 1;                 % EXCITE 15 may have 2 junk views?
else
    header_size = 66072;               % default to CIVM 12.0
    a=sprintf('Unknown software revision, assuming Excite 12.0. BEWARE!\n'); disp(a);
end %if

a=sprintf('GE Excite version %g specified',ge_rev); disp(a);
a=sprintf('Header size = %g',header_size); disp(a);

% Use possible range of filenames to supplied argument into actual filename 
[filename,byte_order,junk_views]=create_filename(p_number);

% now read the raw data into memory
raw=read_pfile(filename,header_size,byte_order);

% report on data found in file
nframes=length(raw)/(2*view_points);         % use length for 1D vectors, not size
a=sprintf('Number of initial frames found in file = %g',nframes);disp(a);

rawR = raw(1:2:end);
rawI = raw(2:2:end);
rdata = complex(rawR,rawI);    % create complex vector of time domain data
% mag_data=abs(rdata);                

% identify and deal with baseline views if not from recon engine
if junk_views==1                 % direct pfile not from recon engine
    base_views=1:2049:nframes;       % 1st and every 2049th frames are baselines
    rdata=reshape(rdata,view_points,nframes); % reshape for accessing baselines
    base_data=rdata(:,base_views);  % this contains only the baselines
    num_bases=length(base_views);       % this is how many baselines there are
    base_data=reshape(base_data,1,[]);   % turn it back into linear array
    
    % plot baseline views if requested
    if plot_baselines==1
        a=sprintf('plotting %g baseline views',num_bases);disp(a);
        figure                    % figure to capture 1st baseline view
        plot(abs(base_data));     % just plot magnitudes of the basline data
        axis tight;
        xlabel('points')
        ylabel('magitude')
        title([filename,'- Baseline Views To Be Removed']);
        box on
        propertyeditor('on');
    end %if

    % now remove baseline views
    a=sprintf('Deleting the following baseline views: %s',int2str(base_views)); disp(a);
    rdata(:,base_views)=[];               % remove baselines from complex data    
    nframes=nframes-length(base_views);   % calculate new nframes
    rdata = reshape(rdata,1,[]);          % turn rdata back into a linear array
    a=sprintf('%g frames remain for processing',nframes); disp(a);
end % baseline views section


% undo loopfactor effects if necessary
if loopfactor~=1
      rdata=undo_loopfactor(rdata,nframes,loopfactor,view_points);
%     mag_data=undo_loopfactor(mag_data,nframes,loopfactor,view_points);
end %if

mag_data=abs(rdata);          % define magnitude data after baseline removal and unlooping

% If there aren't very many fids, adjust display fids to match nframes
if nframes<breath_frames
    breath_frames=nframes;
    if nframes<first_frames
        first_frames=nframes
    end %if
end %if


% separate magnitude plots of various sections of data 
figure   % first number of fids
plot(mag_data(1:first_frames*view_points));    % plot only first_frames
axis tight;
xlabel('points')
ylabel('magnitude')
title([filename,'- Initial FIDs']);
box on
propertyeditor('on');

figure   % first larger number of fids for several breaths
plot(mag_data(1:breath_frames*view_points));    % plot only breath_frames
axis tight;
xlabel('points')
ylabel('magnitude')
title([filename,'- Initial Breaths']);
box on
propertyeditor('on');

figure   % plot all fids
plot(mag_data);    % plot only breath_frames
xlabel('points')
ylabel('magnitude')
title([filename,'- All Magnitude Data']);
box on
propertyeditor('on');

% plot to get both real and imaginary channels
if plot_phases==1
    figure;
    
    subplot(311);
%     plot_section(real(rdata),first_frames,view_points);  % first plot desired frames
    plot(real(rdata));
    grid on;
    xlabel('points');
    ylabel('Real Signal');
    title([filename,'- phase sensitive display']);
    
    subplot(312);
%     plot_section(imag(rdata),first_frames,view_points);
    plot(imag(rdata))
    grid on;
    xlabel('points');
    ylabel('Imaginary Signal');
    
    subplot(313);
%     plot_section(mag_data,first_frames,view_points)   % plot a small section of whole series
    plot(mag_data);
    grid on;
    xlabel('points');
    ylabel('Magnitude');
    propertyeditor('on');
    
end %if

%% calibrate flip angles if specified
if flip_cal==1
    a=sprintf('Calibrating flip angle...');disp(a)
    mag_array=reshape(mag_data,view_points,nframes);  %reshape to find peaks in each view
    peaks=max(mag_array(:,1:nframes));  % find max amplitude in each view
    
    figure
    plot(peaks,'ks')
    box on;
    set(gca,'XMinorTick', 'on')
    set(gca,'YMinorTick', 'on')
    xlabel('frame');
    ylabel('echo peaks');
    title('Flip calibration');
    hold on
    
    %% now curve fit
    %create anonymous function to describe the data%
    fitfunct = @(coefs,xdata)coefs(1)*cos(coefs(2)).^(xdata-1);   % cos theta decay
    guess(1)=peaks(1);
    guess(2)=10*pi/180;       % just guess 10 degrees
    xdata=1:nframes;
    [fitparams,resnorm,residual,exitflag,output,lambda,jacobian]  = lsqcurvefit(fitfunct,guess,xdata,peaks);
    ci = nlparci(fitparams,residual,jacobian);  % returns 95% conf intervals on fitparams by default
    param_err=fitparams-ci(:,1)';
    % focus on flip angle
    flip_angle=fitparams(2)*180/pi;
    flip_err=param_err(2)*180/pi;

    % append fitted curve
    x_fit=linspace(0,nframes,100);    % create
    y_fit=fitfunct(fitparams,x_fit);
    plot(x_fit,y_fit,'-k')               % append the fit
  
    %% DISPLAY FIT COEFFICIENTS AND ERROR BARS
    msg = sprintf('\nFlip angle = %1.3f +-%1.3f',flip_angle,flip_err);
    disp(msg);          % list coefficients as we go

    % display fitted values at 2/3 by 2/3 out in graph
    xloc=round(0.67*nframes);
    yloc=round(0.67*(peaks(1)-peaks(end))+peaks(end));
    text(xloc,yloc,msg); 
end %if

propertyeditor 'on'  % turn on property editor for more axis labelling etc.


