%% using active contour for lung segmentation
% read in data
tic;
uteVol = load_nii('BHUTE_Sub65_P32256_bhute_recon.nii');
uteVol = abs(uteVol.img);

% clustering 
nClusters = 3;
dims = size(uteVol);
disp('Calculating clusters');
[cluster_idx, cluster_center] = kmeans(abs(uteVol(:)),nClusters,'distance','sqEuclidean', ...
    'Replicates',3);
cluster_idx = reshape(cluster_idx,dims);

% User select lung cavity by click at left and right positions
h = imslice(cluster_idx);
ax = gca();
disp('Click on each lung (left/right) in a well insulated slice, then press enter...')
[y_idx, x_idx] = getpts(ax);
x_idx = round(x_idx);
y_idx = round(y_idx);
gui = getappdata(h);
slice_dim = gui.UsedByGUIData_m.dimension_selection.Value;
slice_idx = round(gui.UsedByGUIData_m.slice_slider.Value);
im_val = gui.UsedByGUIData_m.hImage.CData;
lung_idx  = im_val(x_idx(1),y_idx(1));

% acquire initial mask for slice-wise active contour segmentation
% Find lung by connected components, conducted in coronal view
disp('Finding connected lung components');

switch slice_dim  % permute matrix accoridng to user dimension choice
    case 1
        uteVolp = permute(uteVol,[2,3,1]);
    case 2
        uteVolp = permute(uteVol,[3,1,2]);
    case 3
        uteVolp = uteVol;
end

[CC, NUM] = bwlabeln(im_val==lung_idx, 18);
CC_lung_idx1 = CC(x_idx(1),y_idx(1));
CC_lung_idx2 = CC(x_idx(2),y_idx(2));
initial_mask = zeros(size(CC));
initial_mask(CC==CC_lung_idx1|CC==CC_lung_idx2) = 1;

pdims = size(uteVolp);
mask = zeros(pdims);
% active contour for this slice
initial_mask_ac = activecontour(uteVolp(:,:,slice_idx),initial_mask,'Chan-Vese','SmoothFactor',0);
mask(:,:,slice_idx) = initial_mask_ac;

%% forward slice wise active contour segmentation
disp('Conducting forward slice wise active contour segmentation');

switchbuffer = 0; % switch for stop process

current_mask = initial_mask;
current_idx = slice_idx-1;
while current_idx>0 && switchbuffer <2
    current_mask_ac = activecontour(uteVolp(:,:,current_idx),current_mask,'Chan-Vese','SmoothFactor',0.2);
    
    summask = sum(current_mask_ac(:));
    cross_mask = current_mask_ac & current_mask;
    sumcross = sum(cross_mask(:));
    if summask ==0
        current_mask = zeros(size(current_mask));
        switchbuffer = switchbuffer+1;    
    else
        crossratio = sumcross/summask;
        if crossratio>0.4
            mask(:,:,current_idx) = current_mask_ac;
            current_mask = current_mask_ac;
            switchbuffer = 0;
        else
            summask = 0;
            current_mask = zeros(size(current_mask));
            switchbuffer = switchbuffer+1;
        end
    end
    
    current_idx = current_idx-1;
    
%     imshowpair(mask(:,:,current_idx+1),uteVolp(:,:,current_idx));
%     title([num2str(summask),'  ',num2str(crossratio)]);
%     pause(1);
end
%% backward slice wise active contour segmentation 
disp('Conducting backward slice wise active contour segmentation');

switchbuffer = 0; % switch for stop process

current_mask = initial_mask;
current_idx = slice_idx+1;
while current_idx>0 && switchbuffer <2
    current_mask_ac = activecontour(uteVolp(:,:,current_idx),current_mask,'Chan-Vese','SmoothFactor',0.2);
    
    summask = sum(current_mask_ac(:));
    cross_mask = current_mask_ac & current_mask;
    sumcross = sum(cross_mask(:));
    if summask ==0
        current_mask = zeros(size(current_mask));
        switchbuffer = switchbuffer+1;    
    else
        crossratio = sumcross/summask;
        if crossratio>0.4
            mask(:,:,current_idx) = current_mask_ac;
            current_mask = current_mask_ac;
            switchbuffer = 0;
        else
            summask = 0;
            current_mask = zeros(size(current_mask));
            switchbuffer = switchbuffer+1;
        end
    end
    
    current_idx = current_idx+1;
    
%     imshowpair(mask(:,:,current_idx-1),uteVolp(:,:,current_idx));
%     title([num2str(summask),'  ',num2str(crossratio)]);
%     pause(1);
end

switch slice_dim  % permute matrix accoridng to user dimension choice
    case 1
        mask = permute(mask,[3,1,2]);
    case 2
        mask = permute(mask,[2,3,1]);
    case 3
end
toc;