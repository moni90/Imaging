

function Registration_SVD(SM,fname,varargin)
%modified from
% Registration_Segmentation_SVD_021519.m
% Registration_Segmentation_SVD_110518.m
%Do registration only without segementation. Segmentation will be done
%manually.
%
%Hirofumi Nakayama 2020

%Inputs
%SM: cell array containing name of sessions
%fname: header of data file names to be saved
%varargin{1}:taret (index of target session). If not given, target is
%determined based on the alignemnt of blank trials

sn='Registration_SVD.m';
direct=configPath();

%folder path for loading or saving data
path_load_data=direct.local;
path_save_tiff=fullfile(direct.home,'SummaryImages','tiff');
path_save_mat={direct.local,direct.rinberg_data};
path_home = direct.home;
%%
% Check session names using sortSessions.m
% [Sess,SessMouse]=sortSessions();

midline = 128;
target = 0;
if numel(varargin)==1
    target = varargin{1};
end
%%
%Load spatial and temporal component of SVD
% U: Spatial component from all trials
% SV: Temporal component from all trials
% U_empty: Spatial component from blank trials
% SV_empty:Temporal component from blank trials
% Masks: Circular mask of field of view
% trial_info: structure containing information about session

cd(path_load_data)
for sess=1:numel(SM)
    
    load(strcat(SM{sess},'SVD_Compression.mat'),'U','SV','U_empty','SV_empty','Masks','trial_info')
    
    Uall{sess}=U;
    SVall{sess}=SV;
    Ue{sess}=U_empty;
    Ma{sess}=Masks;
    TrialInfo{sess}=trial_info;
end
cd(path_home)
mask=Ma{1}.roi_mask;
clear U SV
SV=SVall;
%%
%Piecewise rigid transformation in frequency domain

%Determine which session to use as a template
if target==0
    %Chose the target session based on the alignment of empty trials
    %Choose the session that have small displacement from other sessiosn
    [Ue_reg,target,Shift]=RegMultipleU(Ue,Ma{1}.roi_mask);
    sprintf('Target calculated from emtpy trials, determined to session %d',target)
else
    sprintf('Target given, session %d', target)
end

%Separate reigid transformation for left and right bulb, using reference
%odors. Use subset of field defined by w x h rectangle for registration.
h=140;
w=100;
mask_L=false(256);mask_L(128-h/2:128+h/2-1,midline-w:midline-1)=true;
mask_R=false(256);mask_R(128-h/2:128+h/2-1,midline+1:midline+w)=true;
num_sess=numel(SM);
%Perform all pairwise alignment between session pairs.
%Perform separate reigstration using each of 4 reference odors to calculate
%amount of shift between each pair
for rs=1:4
    for i=1:num_sess
        for j=setdiff(1:num_sess,i)
            %rs: index of reference odors
            %i: session index of an image to be registered
            %j: session index of tempolate
            s1=TrialInfo{i}.stim_num;
            s2=TrialInfo{j}.stim_num;
            
            f_used=105:150;
            frame_mask1=repmat([f_used]',1,nnz(s1==(max(s1)-4+rs)))+repmat([find(s1==(max(s1)-4+rs))-1]*300,length(f_used),1);
            frame_mask1=frame_mask1(:);
            frame_mask2=repmat([f_used]',1,nnz(s2==(max(s2)-4+rs)))+repmat([find(s2==(max(s2)-4+rs))-1]*300,length(f_used),1);
            frame_mask2=frame_mask2(:);
            
            %Reconstruct average 2d images
            %Excluded 1st and 2nd spatial component for reconstruction
            %because they don't capture glomerular activation.
            im1=reshape(reshape(Uall{i}(:,:,3:end),256^2,[])*mean(SV{i}(frame_mask1,3:end))',256,256);
            im2=reshape(reshape(Uall{j}(:,:,3:end),256^2,[])*mean(SV{j}(frame_mask2,3:end))',256,256);
            
            im1L=im1;im1L(~mask_L)=0;
            im1R=im1;im1R(~mask_R)=0;
            im2L=im2;im2L(~mask_L)=0;
            im2R=im2;im2R(~mask_R)=0;
            
            [~,SL]=dftreg_real(im1L,im2L);
            [~,SR]=dftreg_real(im1R,im2R);
            errorL(i,j,rs)=SL.error;
            errorR(i,j,rs)=SR.error;
            diffphaseL(i,j,rs)=SL.diffphase;
            diffphaseR(i,j,rs)=SR.diffphase;
            row_shiftL(i,j,rs)=SL.row_shift;
            row_shiftR(i,j,rs)=SR.row_shift;
            col_shiftL(i,j,rs)=SL.col_shift;
            col_shiftR(i,j,rs)=SR.col_shift;
        end
    end
end

%Calculate average shift into target sessions across 4 reference odors
for i=1:num_sess
    for j=setdiff(1:num_sess,i)
        %None-zeros diffphse and row_shift>25 col_shift>25 happens when registration is completely wrong
        %avoid above conditions
        od_used=squeeze((diffphaseL(i,j,:)==0)&(diffphaseR(i,j,:)==0)&...
            (row_shiftL(i,j,:)<25)&(row_shiftL(i,j,:)<25)&...
            (col_shiftL(i,j,:)<25)&(col_shiftL(i,j,:)<25));
        Shift_L{i,j}.error=mean(errorL(i,j,od_used),3);
        Shift_L{i,j}.diffphase=mean(diffphaseL(i,j,od_used),3);
        Shift_L{i,j}.row_shift=mean(row_shiftL(i,j,od_used),3);
        Shift_L{i,j}.col_shift=mean(col_shiftL(i,j,od_used),3);
        
        
        od_used=squeeze((diffphaseR(i,j,:)==0)&(diffphaseR(i,j,:)==0)&...
            (row_shiftR(i,j,:)<25)&(row_shiftR(i,j,:)<25)&...
            (col_shiftR(i,j,:)<25)&(col_shiftR(i,j,:)<25));
        Shift_R{i,j}.error=mean(errorR(i,j,od_used),3);
        Shift_R{i,j}.diffphase=mean(diffphaseR(i,j,od_used),3);
        Shift_R{i,j}.row_shift=mean(row_shiftR(i,j,od_used),3);
        Shift_R{i,j}.col_shift=mean(col_shiftR(i,j,od_used),3);
    end
end

% Register spatial components
for sess=setdiff(1:numel(SM),target)
    for s=1:size(Uall{1},3)
        tmp=Uall{sess}(:,:,s);
        tmp(~mask)=mean(tmp(:));
        regL=tmp;
        regL(:,midline:end)=0;
        regL=double(dftreg_real(regL,Shift_L{sess,target}));
        regR=tmp;
        regR(:,1:midline)=0;
        regR=double(dftreg_real(regR,Shift_R{sess,target}));
        
        Ureg{sess}(:,:,s)=regL+regR;
        Ureg{sess}(:,midline-1:midline+1,s)=Uall{sess}(:,midline-1:midline+1,s);
    end
end
Ureg{target}=Uall{target};
%%
%Create tiff stack to check if reference stimuli are correctly aligned
%The first half of frames in the stack is after registration. The second
%half is before registration
for sess=1:numel(SM)
    for rs=1:4
        s1=TrialInfo{sess}.stim_num;
        
        f_used=105:150;
        frame_mask1=repmat([f_used]',1,nnz(s1==(max(s1)-4+rs)))+repmat([find(s1==(max(s1)-4+rs))-1]*300,length(f_used),1);
        frame_mask1=frame_mask1(:);
        
        im1=reshape(reshape(Ureg{sess}(:,:,3:end),256^2,[])*mean(SV{sess}(frame_mask1,3:end))',256,256);
        im2=reshape(reshape(Uall{sess}(:,:,3:end),256^2,[])*mean(SV{sess}(frame_mask1,3:end))',256,256);
        
        img_stack(:,:,sess+numel(SM)*(rs-1))=imadjust(mat2gray(im1));
        img_stack(:,:,4*numel(SM)+sess+numel(SM)*(rs-1))=imadjust(mat2gray(im2));
    end
end

cd(path_save_tiff)
save_tiffstack(img_stack,strcat(fname,'_RefOdors_before_vs_afterRegistration'));
cd(path_home)
%%
%creat tiff stack for manual segmenataion
%Use all odors at multiple timepoints

k=1;
img_stack2=[];
for sess=1:numel(SM)
    str=strsplit(TrialInfo{sess}.tsm_name,'_');
    if contains(str{3},'100hz')
        odt1=str{4}; odt2=str{5};
    else
        odt1=str{3};odt2=str{4};
    end
    
    if TrialInfo{sess}.num_cond==19
        %dilution
        stim_used=[5,10,15];
        
        stim_text=repmat({odt1},1,length(stim_used));
    elseif TrialInfo{sess}.num_cond==12
        %Morphing
        %         stim_used=[2,5,8];
        %         stim_used=2:12;
        stim_used=2:8;
        r=[100,90,71,50,29,10,0];
        stim_text=arrayfun(@(x) sprintf('%s-%d_%s-%d',odt1,x,odt2,100-x),r,'UniformOutput',0);
    elseif TrialInfo{sess}.num_cond==15
        %Odor space mapping
        stim_used=[2:15];
        if iscell(TrialInfo{sess}.OdorNames{1})
            tmp=TrialInfo{sess}.OdorNames(stim_used);
            for i=1:numel(tmp)
                stim_text{i}=strcat(tmp{i}{1},'_',tmp{i}{2});
            end
        else
            stim_text=TrialInfo{sess}.OdorNames(stim_used);
        end
    elseif TrialInfo{sess}.num_cond==10
        stim_used=[2:6];
        stim_text=repmat({odt1},1,length(stim_used));
    elseif TrialInfo{sess}.num_cond==20
        stim_used=2:16;
        stim_text=cellfun(@(x) sprintf('%s_%s',x{1},x{2}),TrialInfo{sess}.OdorNames(stim_used),'UniformOutput',0);
    elseif TrialInfo{sess}.num_cond==29
        stim_used=2:9;
        stim_text=cellfun(@(x) sprintf('%s_%s',x{1},x{2}),TrialInfo{sess}.OdorNames(stim_used),'UniformOutput',0);
    end
    s1=TrialInfo{sess}.stim_num;
    for s=1:length(stim_used)
        
        for phase=1:4
            %inhalation onset is frame 101
            %frame_duration = 10ms
            if phase==1
                f_used=103:106;
            elseif phase==2
                f_used = 107:110;
            elseif phase==3
                f_used=120:140;
            elseif phase==4
                f_used=160:180;
            end
            
            frame_mask1=repmat([f_used]',1,nnz(s1==stim_used(s)))+repmat([find(s1==stim_used(s))-1]*300,length(f_used),1);
            frame_mask1=frame_mask1(:);
            im1=reshape(reshape(Ureg{sess}(:,:,3:end),256^2,[])*mean(SV{sess}(frame_mask1,3:end))',256,256);
            img=imadjust(mat2gray(im1));
            img=rgb2gray(insertText(img,[1,1],stim_text{s}));
            obj=imshow(img);
            ffprintf('sess=%d,s=%d,phase=%d',sess,s,phase)
            img_stack2(:,:,k)=obj.CData;
            k=k+1;
        end
    end
    
    %visualize spatial components of SVD
    %SVD spatial components are not very informative for segmentation.
    %So, commented out so far.
    %     for s=3:6
    %     img = imadjust(mat2gray(Ureg{sess}(:,:,s)));
    %     img=rgb2gray(insertText(img,[1,1],sprintf('SVD_%d',s)));
    %     obj=imshow(img);
    %     img_stack2(:,:,k)=obj.CData;
    %     k=k+1;
    %     end
end

cd(path_save_tiff)
save_tiffstack(img_stack2,strcat(fname,'_ForManualSegmentation'));
cd(path_home)

%%
%Save SVD data alone
U=Ureg(:);
SV=SV(:);
TrialInfo=TrialInfo(:);
SM=SM(:);
reg_params.target_sess = SM{target};
reg_params.target = target;
reg_params.midline = midline;
reg_params.sn = sn;

for p=1:numel(path_save_mat)
    cd(path_save_mat{p})
    save(strcat(fname,'_SVD_registration.mat'),'U','SV','mask','TrialInfo','SM','reg_params','-v7.3')
end

cd(path_home)
end
