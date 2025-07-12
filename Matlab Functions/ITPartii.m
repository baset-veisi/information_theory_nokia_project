
%% constellations 16 QAM MC AND SC

% preperation of complex plane points
N = 16; % Number of the points 
n = ceil(sqrt(N)); % Number of the degrees of freedom
x_16QAM = constellationGen(N); % Constellation points
x_16QAM = x_16QAM(:);

%Getting ready for maximization problem
p0 = ones(n,1)/n;
Aeq = ones(1,n);
beq = 1;
lb = zeros(n,1);
ub = ones(n,1);
SNR_vec = 0:0.1:30;
c_16QAM = zeros(1,length(SNR_vec));
optim_p_16QAM = zeros(n,length(SNR_vec));
%options = optimoptions(@fmincon,'ConstraintTolerance',1e-6,'MaxIterations',1000,'ObjectiveLimit',-1e10,...
  %  'MaxFunctionEvaluations',3000,'StepTolerance',1e-8); % Setting parameters to reduce the convergence time

% Maximizing for each snr
for i = 1:length(SNR_vec)
    fun = @(p) awgn_channel_capacity(SNR_vec(i), x_16QAM, p );
    [optim_p_16QAM(:,i),capacity] = fmincon(fun,p0,[],[],Aeq,beq,lb,ub);
    c_16QAM(i) = -2*capacity - log2(2*pi*exp(1));
end


% Getting ready for maximization problem
p0 = ones(2*n+1,1)/(n);
p0(1) = 0;
Aeq = ones(2,2*n+1);Aeq(:,1) = 0;Aeq(1,n+2:end)=0;Aeq(2,2:n+1) = 0;
beq = ones(2,1);
lb = zeros(2*n+1,1);
ub = ones(2*n+1,1);
c_16QAM_MC = zeros(1,length(SNR_vec));
optim_p_16QAM_MC = zeros(2*n+1,length(SNR_vec));
alpha = 1;
% Maximizing for each snr
for i = 1:length(SNR_vec)
    fun = @(p) (my_awgn_channel_capacity( p(1), SNR_vec(i), x_16QAM, p(2:n+1) )...
        + my_awgn_channel_capacity( 1-p(1),SNR_vec(i)+10*log10(alpha^2), x_16QAM, p(n+2:end) ));
    [optim_p_16QAM_MC(:,i),capacity] = fmincon(fun,p0,[],[],Aeq,beq,lb,ub);
    c_16QAM_MC(i) = -4*capacity - 4*log2(2*pi*exp(1));
end


% figure(1)
% plot(SNR_vec,c_16QAM_MC,"b",SNR_vec,2*ShannonLimitAttenuated(SNR_vec,alpha),"k--");
% grid on;legend('Multi-carrier 16QAM damped');xlabel('SNR');ylabel('Capacity');title('Multi-carrier 16QAM  vs. Shannon limit');
% 

%% constellation MC 64-QAM

% preperation of complex plane points
N = 64; % Number of the points 
n = ceil(sqrt(N)); % Number of the degrees of freedom
x_64QAM = constellationGen(N); % Constellation points
x_64QAM = x_64QAM(:);


% Getting ready for maximization problem
p0 = ones(n,1)/n;
Aeq = ones(1,n);
beq = 1;
lb = zeros(n,1);
ub = ones(n,1);
SNR_vec = 0:0.1:30;
c_64QAM = zeros(1,length(SNR_vec));
optim_p_64QAM = zeros(n,length(SNR_vec));

% % Maximizing for each snr
for i = 1:length(SNR_vec)
    fun = @(p) awgn_channel_capacity(SNR_vec(i), x_64QAM, p );
    [optim_p_64QAM(:,i),capacity] = fmincon(fun,p0,[],[],Aeq,beq,lb,ub);
    c_64QAM(i) = -2*capacity - log2(2*pi*exp(1));
end


% Getting ready for maximization problem
p0 = ones(2*n+1,1)/(n);
p0(1) = 0;
Aeq = ones(2,2*n+1);Aeq(:,1) = 0;Aeq(1,n+2:end)=0;Aeq(2,2:n+1) = 0;
beq = ones(2,1) ;
lb = zeros(2*n+1,1);
ub = ones(2*n+1,1);
c_64QAM_MC = zeros(1,length(SNR_vec));
optim_p_64QAM_MC = zeros(2*n+1,length(SNR_vec));
% Maximizing for each snr
alpha = 1;
% Maximizing for each snr
for i = 1:length(SNR_vec)
    fun = @(p) (my_awgn_channel_capacity( p(1), SNR_vec(i), x_64QAM, p(2:n+1) )...
        + my_awgn_channel_capacity( 1-p(1),SNR_vec(i)+10*log10(alpha^2), x_64QAM, p(n+2:end) ));
    [optim_p_64QAM_MC(:,i),capacity] = fmincon(fun,p0,[],[],Aeq,beq,lb,ub);
    c_64QAM_MC(i) = -4*capacity - 4*log2(2*pi*exp(1));
end





% figure(2)
% plot(SNR_vec, c_64QAM_MC,SNR_vec,2*ShannonLimitAttenuated(SNR_vec,alpha),"k--");
% grid on;legend('64QAM Multi-Carrier','Shannon');xlabel('SNR');ylabel('Capacity');title('Multi-Carrier 64QAM att. ch. vs. its Shannon attenuated limit');
% 

plot(SNR_vec, 0.25*c_64QAM_MC,'r',SNR_vec,c_64QAM,'b');grid on;title('Normal ch, 64 QAM vs. 1 quarter of 64 QAM alpha 1');xlabel('SNR');ylabel('Capacity');legend('Single Carrier','Multi-Carrier')
