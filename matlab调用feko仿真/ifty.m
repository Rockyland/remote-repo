
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse FFT w.r.t. the second variable %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s=ifty(fs,n)
if nargin<2
    n=0;
end

n0=size(fs,2);

Npad=n-n0;
if Npad<=0
    Npad1=0;
    Npad2=0;
else
    if mod(n0,2)==0
        Npad1=fix(Npad/2);
        Npad2=ceil(Npad/2);
    else
        Npad1=ceil(Npad/2);
        Npad2=fix(Npad/2);
    end
end
fs=padarray(fs,[0,Npad1],0,'pre');
fs=padarray(fs,[0,Npad2],0,'post');

%  s=fftshift(ifft(ifftshift(fs.'))).';
 s=fftshift(ifft(ifftshift(fs,2),[],2),2);
