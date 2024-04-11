function [V, d, DV, Dd, D2V, D2d] = blocksVd(h, N, thet)

% A stack of N blocks, each of width 2 units, and height 2h units
% Rotation of each block k relative to the one below it, k-1, is thet(k)
% V is the resulting potenial energy, d is the horizontal displacement of
% the center of the top block. The first and second derivatives are with
% respect to thet. Thus DV and Dd are 1xN and D2V and D2d are NxN

p0 = [1; -h];
R0 = eye(2);
V = 0;
DV = zeros(1, N);
D2V = zeros(N);

Rhat = @(thet_)([cos(thet_) -sin(thet_); sin(thet_) cos(thet_)]);
Rhatp = @(thet_)([-sin(thet_) -cos(thet_); cos(thet_) -sin(thet_)]);

e1 = [1; 0];
e2 = [0; 1];

for k = 0:N-1
    Rhat1 = Rhat(thet(k+1));
    Rhatp1 = Rhatp(thet(k+1));
    Rhatpp1 = -Rhat1;

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Function, V
    %%%%%%%%%%%%%%%%%%%%%%%%%
    R1 = R0*Rhat1;
    p1 = p0 + R0*[-1; h] + R1*[1; h];
    V = V + e2'*p1;

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % First derivative, DV
    %%%%%%%%%%%%%%%%%%%%%%%%%
    DR1 = cell(k+1,1);
    Dp1 = cell(k+1,1);
    for l = 1:k
        DR1{l} = DR0{l}*Rhat1;
        Dp1{l} = Dp0{l} + DR0{l}*[-1; h] + DR1{l}*[1; h];
        DV(l) = DV(l) + e2'*Dp1{l};
    end
    DR1{k+1} = R0*Rhatp1;
    Dp1{k+1} = DR1{k+1}*[1; h];
    DV(k+1) = DV(k+1) + e2'*Dp1{k+1};

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Second derivative, D2V
    %%%%%%%%%%%%%%%%%%%%%%%%%
    D2R1 = cell(k+1,k+1);
    D2p1 = cell(k+1,k+1);
    % m ~= k+1
    for m = 1:k
        for l = 1:k
            D2R1{m,l} = D2R0{m,l}*Rhat1;
            D2p1{m,l} = D2p0{m,l} + D2R0{m,l}*[-1; h] + D2R1{m,l}*[1; h];
            D2V(m,l) = D2V(m,l) + e2'*D2p1{m,l};
        end
        D2R1{m,k+1} = DR0{m}*Rhatp1; % l = k+1
        D2p1{m,k+1} = D2R1{m,k+1}*[1; h];
        D2V(m,k+1) = D2V(m,k+1) + e2'*D2p1{m,k+1};
    end
    % m = k+1, l~=k+1
    for l = 1:k
        D2R1{k+1,l} = DR0{l}*Rhatp1;
        D2p1{k+1,l} = D2R1{k+1,l}*[1; h];
        D2V(k+1,l) = D2V(k+1,l) + e2'*D2p1{k+1,l};
    end
    % m = k+1, l=k+1
    D2R1{k+1,k+1} = R0*Rhatpp1;
    D2p1{k+1,k+1} = D2R1{k+1,k+1}*[1; h];
    D2V(k+1,k+1) = D2V(k+1,k+1) + e2'*D2p1{k+1,k+1};

    R0 = R1;
    p0 = p1;
    DR0 = DR1;
    Dp0 = Dp1;
    D2R0 = D2R1;
    D2p0 = D2p1;
end

% At this stage, p1=pN, Dp1=DpN and D2p1=D2pN
d = 1 - e1'*p1;
Dd = cellfun(@(x_)(-e1'*x_),Dp1)';
D2d = cellfun(@(x_)(-e1'*x_),D2p1);

end
