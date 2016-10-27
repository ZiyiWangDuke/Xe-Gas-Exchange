     % this file calculates RBC:Bar, RBC:Gas, Bar:Gas maps and the average
% ratios. The three inputs are used for comparison
% Ratio maps are calculated by divide them on a pixel base.
% close all;clear all;
%% read in
direct = dir();
nfolder = length(direct)-2;
ratio = struct();

%%
for k = 1:1:nfolder
    folder = direct(k+2).name;
    ratio(k).name = folder;
    cd(folder);
    % read in and preprocess data
    [ gas,dissolved,rbc,barrier,mask,meanRbc2barrier ] = maps_preprocess('parameters');
    
    mask = load_nii('Reg_lungMask_growvent.nii');
    mask = boolean(mask.img);
    % load path data for data saving 
    load('parameters');
    
    TE90 = exp((TE90/1000000)/(2/1000))/exp((TE90/1000000)/0.05);
    flipoff = sin(flipoff/(2.5714*180)*pi)/sin(22/180*pi)*100;

    [pathstr,name,ext] = fileparts(bhute_pfile);
    
    % mask
    dissolved = abs(dissolved).*mask;
    rbc = rbc.*mask;
    barrier = barrier.*mask;
    gas = gas .*mask;
    
    % initial ratio maps
    dis2gas = zeros(size(mask));
    bar2gas = zeros(size(mask));
    rbc2gas = zeros(size(mask));
    rbc2bar = zeros(size(mask));
    
    % small value rule out
    % 0.01 of the average value of the biggest 5% of the data is taken as a
    % threshold
%     [ bmask ] = f_smallruleout( barrier,0.1,0.05);
%     [ gmask ] = f_smallruleout( gas,0.1,0.05);
    % average ratio calculation
    bar2gas(mask) = barrier(mask)*TE90*flipoff./gas(mask);
    rbc2gas(mask) = rbc(mask)*TE90*flipoff./gas(mask);
    rbc2bar(mask) = rbc(mask)./barrier(mask);
    dis2gas(mask) = dissolved(mask)*TE90*flipoff./gas(mask);  
   
%% calculate and save average ratio 
    maskall = sum(mask(:));
    ratio(k).mrbc2bar = meanRbc2barrier;
%     ratio(k).rbc2bar = sum(rbc(:))/sum(barrier(:));
%     ratio(k).rbc2gas = sum(rbc(:))/gasall;
%     ratio(k).bar2gas = sum(barrier(:))/gasall;
%  find(rbc2bar == inf);
    rbc2bar(rbc2bar == inf) = 0;
    ratio(k).rbc2bar = sum(rbc2bar(:))/maskall;
    ratio(k).rbc2gas = sum(rbc2gas(:))/maskall;
    ratio(k).bar2gas = sum(bar2gas(:))/maskall;
    ratio(k).lungvolume = maskall;
    cd('../');
end
