
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward FFT w.r.t. the second variable %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fs=fty(s,n)
if nargin<2
    n=0;
end

n0=size(s,2);

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
s=padarray(s,[0,Npad1],0,'pre');
s=padarray(s,[0,Npad2],0,'post');

%  fs=fftshift(fft(ifftshift(s.'))).';
 fs=fftshift(fft(ifftshift(s,2),[],2),2);
