
%% constellations 16 QAM

% preperation of complex plane points
N = 16; % Number of the points 
z_16QAM = constellationGen(N); % Constellation points
A = sqrt(sum(sum(abs(z_16QAM).^2))/N);
z_16QAM = z_16QAM/A;
z_16QAM = z_16QAM(:);

% Getting ready for maximization problem
p0 = ones(N,1)/N;
Aeq = ones(1,N);
beq = 1;
lb = zeros(N,1);
ub = ones(N,1);
SNR_vec = 0:1:30;
c_16QAM = zeros(1,length(SNR_vec));
optim_p_16QAM = zeros(N,length(SNR_vec));
options = optimoptions(@fmincon,'ConstraintTolerance',1e-4,'MaxIterations',500,'ObjectiveLimit',-1e10,...
    'MaxFunctionEvaluations',800,'StepTolerance',1e-6); % Setting parameters to reduce the convergence time

% Maximizing for each snr
for i = 1:length(SNR_vec)
    fun = @(p) awgn_channel_capacity(SNR_vec(i), z_16QAM, p );
    [optim_p_16QAM(:,i),capacity] = fmincon(fun,p0,[],[],Aeq,beq,lb,ub,[],options);
    c_16QAM(i) = -capacity;
end

figure(1)
plot(SNR_vec,c_16QAM,"r", SNR_vec,log2(1+10.^(SNR_vec/10)),"k--");
grid on;legend('16QAM','Shannon');xlabel('SNR');ylabel('Capacity');title('16QAM vs. Shannon limit');


%% constellation 64 QAM

% preperation of complex plane points
N = 64; % Number of the points 
z_64QAM = constellationGen(N); % Constellation points
A = sqrt(sum(sum(abs(z_64QAM).^2))/N);
z_64QAM = z_64QAM/A;
z_64QAM = z_64QAM(:);

% Getting ready for maximization problem
p0 = ones(N,1)/N;
Aeq = ones(1,N);
beq = 1;
lb = zeros(N,1);
ub = ones(N,1);
SNR_vec = 0:1:30;
c_64QAM = zeros(1,length(SNR_vec));
optim_p_64QAM = zeros(N,length(SNR_vec));
options = optimoptions(@fmincon,'ConstraintTolerance',1e-4,'MaxIterations',500,'ObjectiveLimit',-1e10,...
    'MaxFunctionEvaluations',800,'StepTolerance',1e-6); % Setting parameters to reduce the convergence time

% Maximizing for each snr
for i = 1:length(SNR_vec)
    fun = @(p) awgn_channel_capacity(SNR_vec(i), z_64QAM, p );
    [optim_p_64QAM(:,i),capacity] = fmincon(fun,p0,[],[],Aeq,beq,lb,ub,[],options);
    c_64QAM(i) = -capacity;
end

figure(2)
plot(SNR_vec,c_64QAM,"b", SNR_vec,log2(1+10.^(SNR_vec/10)),"k--");
grid on;legend('64QAM','Shannon');xlabel('SNR');ylabel('Capacity');title('64QAM vs. Shannon limit');

figure(3)
plot(SNR_vec, c_16QAM, 'r', SNR_vec,c_64QAM,"b", SNR_vec,log2(1+10.^(SNR_vec/10)),"k--");
grid on;legend('16QAM','64QAM','Shannon');xlabel('SNR');ylabel('Capacity');title('16QAM vs 64QAM vs. Shannon limit');



