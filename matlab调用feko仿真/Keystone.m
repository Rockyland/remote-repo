%Keystone��ʵ��Խ���뵥Ԫ�㶯У��
function S_key = Keystone(Sf,B,fc)
%Sr��ԭʼ�ز�
%B������
%fc����Ƶ
%S_key:����keystone�任��У���Ļز�
[N,M] = size(Sf);
posa = B/(fc*N);

hWaitbar = waitbar(0, '�ز����㣬wait...', 'CreateCancelBtn', 'delete(gcbf)');
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
    waitbar(compT, hWaitbar, ['Keystone����...', num2str(round(compT, 2) * 100), '%']);
end
close(hWaitbar)
S_key = ifty(DFTtemp);