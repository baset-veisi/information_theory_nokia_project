function c = my_awgn_channel_capacity(  beta, SNR, x, p )
% CAPACITY_AWGN computes the capacity of the Gaussian channel with power 
% constraint, where the signal-to-noise ratio is SNR (in dB) and the input 
% is distributed over x (the constellation) with probability density function
% p, i.e., the probability to observe x(i) is p(i).
%
% By construction, x and p are vectors that should have the same length.
%
% The code works in dimension 1 and 2. 
%
% The channel model is:
% Y = sqrt(snr)X+N
% where E[X^2]<1 and N ~ N(0,1). Here snr is in linear scale. The capacity 
% is equal to I(X;Y) where I is the mutual information. The mutual
% information computation gives (H is the entropy function):
% I(X;Y) = H(Y) - H(Y|X)
%        = H(Y) - H(sqrt(snr)X+N|X)
%        = H(Y) - H(N)
% For a one dimensional Gaussian N(mu, sigma^2), the entropy is equal to
% 1/2*log2(2*pi*e*sigma^2). To compute H(Y), we compute the pdf of Y as 
% follow:
% f(y) = sum(a in x) Pr(x=a) Pr(y|x=a).
%
% Arguments (input):
% SNR - Signal-to-noise ratio [dB]
% x   - Input support, in other word the constellation [vector] 
% p   - Input distribution [vector]
%
% Arguments (output):
% c - Capacity of the Gaussian channel
%
% Examples:
% For the uniformly distributed BPSK constellation, the capacity at 10 dB is 
% c = capacity_awgn(10, [-1, 1], [0.5 0.5])
% 
% For the uniformly distributed M-QAM at 5 dB:
% x = (-sqrt(M)+1) : 2 : (sqrt(M)-1)    
% x = ones(sqrt(M), 1) * t              
% x = x - 1i*x'                         
% x = x / sqrt(1/M*norm(x, 'fro')^2)  -> Normalise the constellation power  
% x = reshape(X, 1, M)                  
% c = capacity_awgn(5, x, 1/M*ones(1,M))
%
% Author      : Hugo MERIC
% Homepage    : http://hugo.meric.perso.sfr.fr/index.html
% Release     : 1.0
% Release date: 2015-04-16
%alpha = p(1);


% --------------------------
% ----- Initialization -----
% --------------------------

snr = 10^(SNR/10) ;
if beta < 0.01
   c = 0;
else
    snr = beta*snr;
    E_0 = sum(abs(x').^2*p);
    snr = snr/E_0;
    % --------------------
    % ----- Capacity -----
    % --------------------
    
    
    func = int(snr,x,p);
    entY = integral(func, -50*max(x), 50*max(x));
    c =  entY;



    
    
    % ---------------------
    % ----- Functions -----
    % ---------------------
        
end
function z = pdf_channel_output(y, snr, x, p)
            z = 0 ;
            for k=1:length(x)
                z = z + p(k) .* sqrt(1/(2*pi)).*exp( -((y-sqrt(snr).*x(k)).^2)/2 );
            end
        end
    
function z = int(snr, x, p)
            z = @(y) pdf_channel_output(y, snr, x, p) .* log2(pdf_channel_output(y, snr, x, p)+1e-50) ;
end

end




