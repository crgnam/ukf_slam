% pcov_pos = posSemiDefCov(pcov)
%
% This function reads in a covariance matrix and checks if it is positive
% definite.  If it is not, but non-positive eigenvalues are very small in
% magnitude, those eigenvalues are set to a small positive number.

function pcov_pos = posSemiDefCov(pcov)

[~, err] = cholcov_codegen(pcov);
% err is zero if the matrix is positive semidefinite

if err ~= 0 %the covariance matrix is not positive semidefinite
    [v, d]=eig(pcov);
    
    % set any eigenvalues that are not positive to a small positive number
    for ii=1:size(d,1)
        if d(ii,ii)<=1e-14
            if abs(d(ii,ii)/max(max(d)))>1e-3 %large negative eigenvalues
                seven = 7;
                %error('Matrix is negative (semi)definite')
            end
            d(ii,ii)=1e-14;
        end
    end
    
    % Reconstruct the covariance matrix, which should now be positive
    % definite
    pcov_pos=real(v*d*v');
    %disp('Not Pos')
else %already positive definite
    pcov_pos=real(pcov);
    %disp('Pos')
end
