function [out, rates, b] = print_summary_chain(chain, epsigma, show, exchange)
% Function to print a summary of chain composition
% out:  K-by-3 or K-by-5 table 
if ~exist('exchange', 'var')
    exchange=0;
end
Rej = chain{1}; K = size(Rej, 1); 
n = chain{2}; ne = chain{3};
b=n-ne+1; 
epsilon = epsigma{1}; SIGMA = epsigma{2};

if(exchange)
    nsw = chain{4}; sw = chain{5};
    % Exchange move acceptance rate across all chains
    fprintf('%f percent exchange moves accepted \n', sw/nsw*100)
end

if(exchange)
    out = zeros(K, 6);
    rates = zeros(K, 2);
else
    out = zeros(K, 3);
    rates = zeros(K, 1);
end

for k=1:K
    % Rejection due to prior
    out(k, 1) = sum(Rej(k, 1:n(k))==2)/n(k)*100;
    % Rejection due to race
    out(k, 2) = sum(Rej(k, 1:n(k))==1)/n(k)*100;
    % Accepted local move
    out(k, 3) = sum(Rej(k, 1:n(k))==0)/n(k)*100;
    % Local move acceptance rate
    rates(k, 1) = sum(Rej(k, 1:n(k))==0)/sum(Rej(k, 1:n(k))<3)*100;
    
    if(exchange)
        % Rejected exchange move
        out(k, 4) = sum(Rej(k, 1:n(k))==3)/n(k)*100;
        % Accepted exchange move
        out(k, 5) = sum(Rej(k, 1:n(k))==4)/n(k)*100;
        % Total percentage of exchange move
        out(k, 6) = out(k, 4)+ out(k, 5);
        % Exchange move acceptance rate
        rates(k, 2) = sum(Rej(k, 1:n(k))==4)/sum(Rej(k, 1:n(k))>=3)*100;        
    end
    
    if(show)
        fprintf(' \n --- CHAIN %d, epsilon = %.2f --- \n \n', k, epsilon(k))
        fprintf('SIGMA:  ');     disp(SIGMA(:,k)')
        fprintf('Chain is made up of: \n')
        fprintf('Rejection due to prior: %f percent |', out(k, 1))
        fprintf(' due to race: %f percent \n', out(k, 2))
        fprintf('Accepted local move: %f percent \n',out(k, 3))
        if(exchange)
            fprintf('Rejected swap: %f percent |',out(k, 4))
            fprintf(' Accepted swap: %f percent \n',out(k, 5))
        end
        fprintf('Acceptance rates: \n')
        fprintf('Local moves rate: %f percent ', rates(k, 1))
        if(exchange)
            fprintf('| Exchange moves rate: %f percent \n', rates(k, 2))
        end
    end
end

fprintf('\nSample sizes:  ')
disp(b')
if(exchange)
    % Exchange moves across all chains
    fprintf('Percentage exchange moves: %f \n \n',nsw/sum(n)*100)
end
end