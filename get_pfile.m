function get_pfile(source,run_number);
% function to transfer a p-file from the recon engine to present director
% Created 1/19/08
% Modified 3/4/08 to also get p-files direct from scanner
% Updated 9/17/08 to get p-files from scanner beyond only the latest one
% use as follows
% get_pfile('bullseye','testO');
% get_pfile('lx-ron1','P23516.7');
% get_pfile('lx-ron1','+');  % to get the latest


radish=1;           % set to 1 for new radish recons

% figure out if it's a scanner or a recon engine
if strmatch(source, strvcat('onnes', 'kamy', 'heike', 'lx-ron1'))>0       % it's a scanner
    scanner=1;
    % assign user name, password and path for scanners
    user='sdc';
    pass='adw2.0';
    full_path='/usr/g/mrraw'
    % assign right ftp addresses
    if strcmp(source,'onnes')
        source='152.16.234.8';
    elseif strcmp(source,'kamy')
        source='152.16.234.9';
    elseif strcmp(source,'heike')
        source='152.16.234.10';      
    else
        source='152.16.249.132';      % rad onc scanner
    end % if
elseif strmatch(source, strvcat('bullseye','jessie','manta','tinos'))>0    % it's a recon engine
    scanner=0;
    % assing user name and password for recon engines
    user='omega';
    pass='4.signa!';
    % set paths
    if strcmp(source,'bullseye')
        path='/reconbu';
    elseif strcmp(source,'jessie')
        path='/reconje';
    elseif strcmp(source,'tinos')
        path='/analyzet';
    else
        path='/reconma';
    end % if
    % get filenames right
    if radish==1                % data was reconned with radish, not poodle
        full_path=[path '/' run_number '.work'];
        file_name=[run_number '.rp'];
    else                        % data was reconned with poodle
        full_path=[path '/' run_number '/raw_data'];
    end %if
else
    a=sprintf('Source not recognized. Exiting');disp(a)
    return
end %if

fclose all;

% notify user where we're searching
a=sprintf('searching %s',full_path);disp(a)

% ftp over and get a directory listing in order of most recent
ftp_id=ftp(source,user,pass);
cd(ftp_id,full_path);            

% go either scanner route or recon engine route
if scanner==1   
    if run_number=='+'                   % user just wants latest p-file
        d=dir(ftp_id,'-ltr P*');
        a=sprintf('latest file is');disp(a);
        disp(d(end))
        file_name = d(end).name;
        mget(ftp_id,file_name);
    else                                %user wants specific p-file
        file_name=run_number;
        mget(ftp_id,file_name);
    end %if
else                                  % recon engine case
    if radish==1
        mget(ftp_id,file_name);
    else
        d=dir(ftp_id);
        a=sprintf('latest file is');disp(a);
        disp(d(1))
        file_name = d(1).name;
        mget(ftp_id,file_name);
    end %if
end %if

a=sprintf('%s retrieved. Logging off',file_name);disp(a);
close(ftp_id);
