function [trialframes]=read_tsm_int16(filename,range)
%Modified from read_tms.m trialframes are stored in int16 not double
%This file is to read .tsm file of imaging data from RedShirt NeuroCCD
%INPUT
%filename: name of binary file read
%range: range of frames in .tsm file read
%
% header of .tsm file compose of 2880 characters (80 characters x 36 lines)
% you need to extract width (NAXIS1), hight (NAXIS2) and # of frames (NAXIS3) from the header
% you can also extract other info for imaging
% frame rate (EXPOSURE x 1000 ms), 
% gain setting(High: SATURATE = 16383, Med: SATURATE = 8191, Low: SATURATE = 4095)
%
% example header
% --------------------------------------
% SIMPLE  =                    T
% BITPIX  =                   16
% NAXIS   =                    2
% NAXIS1  =                  256
% NAXIS2  =                  256
% NAXIS3  =                  400
% DATE    = 2013-04-21T10:20:40
% ORIGIN  = SciMeasure Analytical Systems, Inc.
% CREATOR = Turbo-SM
% INSTRUME= CCD67-VER1
% EXPOSURE=         0.0100000000
% SATURATE=                16383
% DATAMAX =                16383
% DATAMIN =                    0
% END
% ----------------------------------------
% Hirofumi Nakayama 2017
%

%% set path
path_original = pwd;
try
    [direct]=configPath();
    path_dir = fullfile(direct.rinberg_smdata,filename(1:end-7));
catch
    path_dir = fullfile(path_original, filename(1:end-7));
end

if length(range) == 2
    range = range(1):range(2);
end

%% extracting file content
fullpath=fullfile(path_dir,filename);
fid= fopen(fullpath,'rb');
h_line = fread(fid, 2880, '*char')';
pixels = fread(fid, 'int16');
fclose(fid);

%% reading header
max_head=2880/80;
header=cell(max_head,1);
for i=1:max_head
    header{i}=strrep(h_line(1+(i-1)*80:i*80),' ','');
    
    if ~isempty(strfind(header{i},'NAXIS1'))
        width=str2double(header{i}(8:end));
    elseif ~isempty(strfind(header{i},'NAXIS2'))
        height=str2double(header{i}(8:end));
    elseif ~isempty(strfind(header{i},'NAXIS3'))
        nFrames=str2double(header{i}(8:end));
    end
end

%% reshaping pixel values into frames
nPixels = width * height;

% Last frame: pixels((nFrames)*nPixels+1:(nFrames+1)*nPixels) is the dark frame
darkframe=int16(pixels((nFrames)*nPixels+1:(nFrames+1)*nPixels));
darkframe=reshape(darkframe,[width,height])';

if ~exist('range', 'var')
    range = 1:nFrames;
end

trialframes=zeros(height,width,length(range),'int16');

for frameno=1:length(range)
    frame=int16(pixels((range(frameno)-1)*nPixels+1:range(frameno)*nPixels));
    frame=reshape(frame,[width,height])';
    
%     subtraction of darkframe
    trialframes(:,:,frameno)=frame-darkframe;
end
