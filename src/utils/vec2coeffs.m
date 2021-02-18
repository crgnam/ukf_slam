function [Cnm,Snm] = vec2coeffs(Cnm_Snm_vec)
    % Calculate size of the Cnm coefficient matrix:
    L = size(Cnm_Snm_vec,1);
    for ii = 1:2000
        if L == eval(ii)
            N = ii+1;
            break
        end
    end
    Cnm = zeros(N);
    Snm = zeros(N);
    
    % Split appropriately:
    num_Snm = (L-N)/2 + 1;
    num_Cnm = L - num_Snm;
    Cnm_vec = Cnm_Snm_vec(1:num_Cnm);
    Snm_vec = Cnm_Snm_vec(num_Cnm+1:end);

    % Scale factor to try to maintain numerical stability:
    scale = 1e-7;
    
    kk = 1;
    for ii = 3:N
       for jj = 1:ii
          Cnm(ii,jj) = Cnm_vec(kk)*scale;
          kk = kk+1;
       end
    end
    Cnm(1) = 1;
    
    kk = 1;
    for ii = 3:N
        for jj = 2:ii
            Snm(ii,jj) = Snm_vec(kk)*scale;
            kk = kk+1;
        end
    end
end

% Function to get the appropriate size (guess and check is best approach
% for now):
function [val] = eval(n)
    step = 5;
    val = 0;
    for ii = 1:n-1
        val = val+step;
        step = step+2;
    end
end