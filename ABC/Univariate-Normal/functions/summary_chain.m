function [out, rates] = summary_chain(chain, eprho, show)
% Function to print a summary of chain composition
% out:  K-by-5 table 

R_in = chain{1}; K = size(R_in, 2);
b = chain{2}; nsw = chain{3};
sw = chain{4};
epsilon = eprho{1}; rho = eprho{2};

% Exchange move acceptance rate across all chains
fprintf('%f percent exchange moves accepted \n', sw./nsw*100)


out = zeros(K, 5);
rates = zeros(K, 2);


for k=1:K
    % Rejection due to prior
    out(k, 1) = sum(R_in{k}==2)/b(k)*100;
    % Rejection due to race
    out(k, 2) = sum(R_in{k}==1)/b(k)*100;
    % Accepted local move
    out(k, 3) = sum(R_in{k}==0)/b(k)*100;
    % Local move acceptance rate
    rates(k, 1) = sum(R_in{k}==0)/sum(R_in{k}<3)*100;
    
    % Rejected exchange move
    out(k, 4) = sum(R_in{k}==3)/b(k)*100;
    % Accepted exchange move
    out(k, 5) = sum(R_in{k}==4)/b(k)*100;
    % Exchange move acceptance rate
    rates(k, 2) = sum(R_in{k}==4)/sum(R_in{k}>=3)*100;
    
    if(show)
        fprintf(' \n --- CHAIN %d, epsilon = %.2f --- \n \n', k, epsilon(k))
        fprintf('rho:  ');     disp(rho)
        fprintf('Chain is made up of: \n')
        fprintf('Rejection due to prior: %f percent |', out(k, 1))
        fprintf(' due to race: %f percent \n', out(k, 2))
        fprintf('Accepted local move: %f percent \n',out(k, 3))
        fprintf('Rejected swap: %f percent |',out(k, 4))
        fprintf(' Accepted swap: %f percent \n',out(k, 5))
        fprintf('Acceptance rates: \n')
        fprintf('Local moves rate: %f percent ', rates(k, 1))
        fprintf('| Exchange moves rate: %f percent \n', rates(k, 2))
    end
end

fprintf('\nSample sizes:  ')
disp(b')
% Exchange moves across all chains
fprintf('Percentage exchange moves: %f \n \n',nsw'/sum(b)*100)

end