function [U,SV]=splitSVD(frames,trial_info,mask,varargin)
%This function is to conduct SVD cleaning to video frames
%Based on
% Spontaneous behaviors drive multidimensional, brain-wide neural activity
% Carsen Stringer*1;2, Marius Pachitariu*1;3;4, Nicholas Steinmetz5, Charu Reddy5, Matteo Carandiniy5 and
% Kenneth D. Harrisy3;4
%
%a session is split into multiple blocks and SVD are performed for each
%block.
%Concatenate all spatial factors U from each block
%Apply SVD into concatenated U to obtain the estimate of U to original
%frame
%S*V' is given by the projection of frames into U
%Input
%frames:pixels x concatenated frames
%(reshape to 2d if frames are 256x256xframes
%Outpu
%U:pixel x num_sval
%V:temporal factors weighted by svals
%frames x num_sval (S*V')' of svd)
%2019 Hirofumi Nakayama

if iscell(frames)
    num_trials=numel(frames);
    
    if ndims(frames{1})==3
        if size(frames{1},1)==256&&size(frames{1},2)==256
            %each component in a cell is 256x256xframes
            frames=cell2mat(cellfun(@(x) reshape(x,[],size(x,3)),frames,'UniformOutput',0)');
        else
            error('check frames, ndims=%d,size(frames{1})=%d,%d,%d',ndims(frames{1}),size(frames{1},1),size(frames{1},2),size(frames{1},3))
        end
    elseif ndims(frames{1})==2
        %each component in a cell is pixels x frames
        frames=cat(2,frames{:});
    end
elseif ndims(frames)==3
    %This will fail if a subset of trials in a session is given
    num_trials=length(trial_info.inh_onset);%This may need to be written in a different way
    
    if size(frames,1)==256&&size(frames,2)==256
        frames=reshape(frames,[],size(frames,3));
    else
        error(sprintf('check frames, ndims=%d,size(frames{1})=%d,%d,%d',ndims(frames{1}),size(frames{1},1),size(frames{1},2),size(frames{1},3)))
    end
end


block_size=30;
%block_size=25;
if numel(varargin)==1
    num_svals=varargin{1};
elseif numel(varargin)==2
    [num_svals,bwf]=varargin{1};%Butterworth filter to V
elseif numel(varargin)==3
    [num_svals,bwf,fig_plot]=varargin{:};
elseif numel(varargin)==3
    [num_svals,bwf,fig_plot,block_size]=varargin{:};
end

fps=trial_info.fps;
frame_dur=floor(1000/fps);

highpassCutoff = 0.01;% hz
%num_svals=300;
fnum=size(frames,2)/num_trials;


%Todo: Need to handle trials_unread!!!!!!!!!!!!!!!!!!!!!!!!
if num_trials<=block_size
    num_block=1;
    trials=[0,num_trials];
else
    num_block=floor(num_trials/block_size);
    trials=0:block_size:num_trials;trials(end)=num_trials;
end

trials_read = ~logical(trial_info.trials_unread);
frames_read = repmat(trials_read,fnum,1);frames_read=frames_read(:);
for i=1:num_block
    tic
    subF=single(frames(:,1+trials(i)*fnum:trials(i+1)*fnum));
    %Remove unread trials from svd
    subF=subF(:,frames_read(1+trials(i)*fnum:trials(i+1)*fnum));
    subF=subF(:,(max(subF)-min(subF))~=0);%Remove all 0 frames which correspond to skipped trials
    subF=subF(mask(:),:);
    [Ub,Sb,~] = svd(subF,'econ');
    G{i}=Ub(:,1:num_svals)*Sb(1:num_svals,1:num_svals);
    sprintf('svd for block %d /%d = %.2f sec',i,num_block,toc)
end

clear subF

G_all=cell2mat(G);

[U,~,~] = svd(G_all,'econ');
SV=single(frames(mask(:),:))'*U;

U=U(:,1:num_svals);U=roimask2full(U,mask);
SV=SV(:,1:num_svals);

%U(:,1):black pigment,
%U(:,2):baseline edges(may be because of motion),
%U(:,3):Midline blood vessel,
%U(:,4:5):Odors
%Therefore, reconstruct movies from U(:,4:5)

%Butterworth filter if necessary (Slghtly helps to reduce speckle noies but no significantly)
if bwf==1
    [b100s, a100s] = butter(2, highpassCutoff/(fps/2), 'high');
    dV = detrend(SV, 'linear')';
    SV = filter(b100s,a100s,dV,[],2)';
end

if fig_plot==1
    figure2;
    for i=1:20
        subplot(4,5,i)
        g1=U(:,:,i);
        g1(~mask)=mean(g1(mask));
        imshow(mat2gray(g1));title(num2str(i))
    end
end
