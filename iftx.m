
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse FFT w.r.t. the first variable %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s=iftx(fs,n)
if nargin<2
    n=0;
end

n0=size(fs,1);

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

fs=padarray(fs,[Npad1,0],0,'pre');
fs=padarray(fs,[Npad2,0],0,'post');

 s=fftshift(ifft(ifftshift(fs)));

