%Keystone法实现越距离单元徙动校正
function S_key = Keystone(Sf,B,fc)
%Sr：原始回波
%B：带宽
%fc：载频
%S_key:经过keystone变换的校正的回波
[N,M] = size(Sf);
posa = B/(fc*N);

hWaitbar = waitbar(0, '回波计算，wait...', 'CreateCancelBtn', 'delete(gcbf)');
set(hWaitbar, 'Color', [0.9, 0.9, 0.9]);
%
hBtn = findall(hWaitbar, 'type', 'Uicontrol');
set(hBtn, 'String', 'cancle', 'FontSize', 10);

for l=-M/2:M/2-1
    stemp=zeros(N,1);
    for m = -M/2:M/2-1
        stemp = stemp+Sf(1:N,m+M/2+1).*exp(-j*(2*pi*(1+posa*[-N/2:N/2-1].')/M)*l*m);
    end
    DFTtemp(1:N,l+M/2+1) = stemp;

    compT = (l+M/2+1) / M;
    waitbar(compT, hWaitbar, ['Keystone进度...', num2str(round(compT, 2) * 100), '%']);
end
close(hWaitbar)
S_key = ifty(DFTtemp);