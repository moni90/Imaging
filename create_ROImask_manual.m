function [roiMask_id,roiMask_stack]=create_ROImask_manual(name,opt)
%Function to create binary roi mask from RoiSet.zip files created in FIJI
%ImageJ
%Hirofumi Nakayama 2020

%%
%set up paths
path_original = pwd;
try
    directory=configPath();
    path_imaging = directory.home_2p;
catch
    path_imaging = pwd;
end

%set parameters
bmp = false;
img_size = 512;
strROIArchiveFilename = name;
if exist('opt','var')
    if isfield(opt,'bmp')
        bmp = opt.bmp;
    end
    
    if isfield(opt,'path')
        if contains(name,'RoiSet','IgnoreCase',true)
            %If name of RoiSet is given
            strROIArchiveFilename=fullfile(opt.path,sprintf('%s.zip',name));
        else
            %If name of tiff stack is given
            strROIArchiveFilename=fullfile(opt.path,sprintf('%sRoiSet.zip',name));
        end
    else
        if contains(name,'RoiSet','IgnoreCase',true)
            %If name of RoiSet is given
            strROIArchiveFilename=fullfile(path_imaging,'Segmentation\RoiSet',sprintf('%s.zip',name));
        else
            %If name of tiff stack is given
            strROIArchiveFilename=fullfile(path_imaging,'Segmentation\RoiSet',sprintf('%sRoiSet.zip',name));
        end
    end
    
    if isfield(opt, 'img_size')
        img_size = opt.img_size;
    end
end

%%
%Read RoiSet.zip file
%sROI{i}.vnRectBounds=[Top,Left,Bottom,Right]
[sROI] = ReadImageJROI(strROIArchiveFilename);

%Create binary masks for individual roi
[cc,rr] = meshgrid(1:img_size);
roiMasks=[];
for i=1:numel(sROI)
    if isequal(sROI{i}.strType,'Freehand')
        ind=sROI{i}.mnCoordinates;
        
        %For whatever reasons, ind take values of 0-img_size not 1-img_size.
        %replace 0 in index to 1 to avoid errors
        ind(ind==0)=1;
        
        tmp=false(img_size);tmp(sub2ind([img_size,img_size],ind(:,2),ind(:,1)))=true;
        centOfMass(i,:) = [mean(ind(:,2)), mean(ind(:,1))];
        tmp=imfill(tmp,'holes');
        %         tmp=closeOpenROI(tmp);
        cc= bwconncomp(tmp);
        
        if cc.NumObjects==1
            tmp=closeOpenROI(tmp);
        elseif cc.NumObjects>=2
            tmp=(poly2mask(ind(:,1),ind(:,2),img_size,img_size));
        end
        roiMasks(:,:,i)=tmp;
        % close open loop
        
    elseif isequal(sROI{i}.strType,'Oval')
        r=mean(sROI{i}.vnRectBounds([1,3]));
        c=mean(sROI{i}.vnRectBounds([2,4]));
        rw=(sROI{i}.vnRectBounds(3)-sROI{i}.vnRectBounds(1))/2;
        cw=(sROI{i}.vnRectBounds(4)-sROI{i}.vnRectBounds(2))/2;
        centOfMass(i,:) = [r, c];
        
        %Create ellipse masks from bounding box
        roiMasks(:,:,i) = sqrt(((rr-r)/rw).^2+((cc-c)/cw).^2)<1;
        
    end
end

%Detect duplicates in roiMasks and delete that
ind_exclude = [];
cell_vec = reshape(roiMasks, size(roiMasks,1)^2,[]);
for i=1:size(roiMasks,3)
    tmp = cell_vec - cell_vec(:,i);
    ind = find(~any(tmp));
    if length(ind)>=2
        ind_exclude = [ind_exclude, ind(2:end)];
    end
end
roiMasks(:,:,unique(ind_exclude)) = [];
centOfMass(unique(ind_exclude),:)=[];

%Use this not LabelingBinary.m to deal with adjuscent ROIs
roiMask_id=zeros(size(roiMasks(:,:,1)));

%rank order center of mass in ascending order
[~,p] = sort(centOfMass(:,1),'ascend');
posRank = 1:length(p);
posRank(p) = posRank;

%Create single images containing id of roiMasks
for i=1:length(posRank)
    roiMask_id(logical(roiMasks(:,:,posRank==i)))=i;
end

%Exclulde pixels that are shared by multiple ROIs
roiMask_id(sum(roiMasks,3)~=1)=0;

% figure;imagesc(logical(roiMask));
for i=1:numel(sROI)
    roiMask_stack(:,:,i)=logical(roiMask_id==i);
end

%Create .bmp files for individual ROIs
if bmp
    new_folder = fullfile(path_imaging,'Segmentation','bmp');
    mkdir(new_folder)
    cd(new_folder)
    for i = 1:size(roiMask_stack,3)
        imwrite(1-roiMask_stack(:,:,i),sprintf('%s_roi_bmp_%d.bmp',name,i),'bmp')
    end
end

cd(path_original)