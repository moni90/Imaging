function [img_reg,varargout]=dftreg_real(img,varargin)
%This function performs dftregistration

% %[output, Greg] = dftregistration(buf1ft,buf2ft,usfac)
% % Copyright (c) 2016, Manuel Guizar Sicairos, James R. Fienup, University of Rochester
% % All rights reserved.
% %
% % Redistribution and use in source and binary forms, with or without
% % modification, are permitted provided that the following conditions are
% % met:
% %
% %     * Redistributions of source code must retain the above copyright
% %       notice, this list of conditions and the following disclaimer.
% %     * Redistributions in binary form must reproduce the above copyright
% %       notice, this list of conditions and the following disclaimer in
% %       the documentation and/or other materials provided with the distribution
% %     * Neither the name of the University of Rochester nor the names
% %       of its contributors may be used to endorse or promote products derived
% %       from this software without specific prior written permission.
% %
% % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% % ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% % LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% % CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% % SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% % INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% % CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% % ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% % POSSIBILITY OF SUCH DAMAGE.

%Modified above function to make both input and output as real image not
%fft
%Also return variable 'shifts' that can be used to transform other images
%Inputs
%img: image to be registered
%varargin: target image or shifts imformation
%Outputs
%img_reg: reisgeterd image
%varargout{1}:shift
%
%Hirofumi Nakayama 2019

if numel(varargin)>=2
    error('Wrong number of input arguments %d',numel(varargin))
end

if isequal(size(img),size(varargin{1}))
    target=varargin{1};
    usfac=1;
    [output, Greg] = dftregistration(fft2(target),fft2(img),usfac);
    shifts.error=output(1);
    shifts.diffphase=output(2);
    shifts.row_shift=output(3);
    shifts.col_shift=output(4);
    varargout{1}=shifts;
else
    
    shifts=varargin{1};
    [nr,nc]=size(img);
    Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
    Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);
    
    diffphase=shifts.diffphase;
    row_shift=shifts.row_shift;
    col_shift=shifts.col_shift;
    
    [Nc,Nr] = meshgrid(Nc,Nr);
    Greg = fft2(img).*exp(1i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
    Greg = Greg*exp(1i*diffphase);
end

img_reg=real(ifft2(Greg));


