function [xhat,ck]=ffs(xt,t,n,T)
    %Inputs: 
    % xt- original signal (time-limited)
    % t- time vector 
    % n- number of harmonics (positive)
    % T- period to assume 
    dt=t(2)-t(1);
    L=length(t);
    ck= zeros (2*n+1,1); %coefficients
    k_vals=-n:n;
    for idx=1:length(k_vals)
        k=k_vals(idx);
        ck(idx)=(1/T)*sum(xt .* exp(-1j*2*pi*k*t/T))*dt; %numerical integration
    end
%Now we need to reconstruct the signal from the coefficients 
    xhat=zeros(size(t));
    for idx=1:length(k_vals)
        k=k_vals(idx);
        xhat=xhat+ck(idx)*exp(1j*2*pi*k*t/T);
    end
    xhat=real(xhat);
end

