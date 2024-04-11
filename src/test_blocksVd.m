h = 1;
N = 3;

thet = rand(N,1)*(pi/180)*1; % random between 0 and 1 degree

[V,d,DV,Dd,D2V,D2d] = blocksVd(h, N, thet);

% compute finite difference (FD) approx of first and second derivatives
DVFD = zeros(1,N);
DdFD = zeros(1,N);
D2VFD = zeros(N);
D2dFD = zeros(N);

for k = 1:N
    dthet = zeros(N,1);
    dthet(k) = 1e-6;
    [V1,d1,DV1,Dd1,~,~] = blocksVd(h, N, thet+dthet);
    DVFD(k) = (V1-V)/1e-6;
    DdFD(k) = (d1-d)/1e-6;
    D2VFD(k,:) = (DV1-DV)/1e-6;
    D2dFD(k,:) = (Dd1-Dd)/1e-6;
end

disp([DV; DVFD])

disp([Dd; DdFD])

disp([D2V D2VFD])

disp([D2d D2dFD])