function V_RF = change_V_RF(beta,K,P,sigma2)
global N N_RF
global H_t V_RF

% A = cell(1,9);
% B = cell(1,9);
% D = cell(1,9);
% zeta_B = zeros(N,N_RF);
% zeta_D = zeros(N,N_RF);
% eta_B = zeros(N,N_RF);
% eta_D = zeros(N,N_RF);
% c = zeros(N,N_RF);
% z = zeros(N,N_RF);
% phi = zeros(N,N_RF);

for j = 1:N_RF
    % 计算A_j
    V_RF_bar = V_RF;
    V_RF_bar(:,j) = [];
    A{j} = H_t * V_RF_bar*(V_RF_bar')* H_t';
    % 计算B_j,D_j,用于计算zeta和eta
    A_inv = A{j}^(-1);
    B{j} = H_t'* A_inv^2 * H_t;
    D{j} = H_t'* A_inv * H_t;
    for i = 1:N
        % 计算zeta和eta
%         temp = 0;
%         for m = 1:1:N
%             if( m~=i )
%                 for n = 1:1:N
%                     if(n~=i)
%                         temp = temp + conj(V_RF(m,j))*B{j}(m,n)*V_RF(n,j);
%                     end
%                 end
%             end
%         end
        % 计算zeta_B
        V_RF_i = V_RF;
        V_RF_i(i,:) = []; % 去掉第i行
        B_temp = B{j};
        B_temp(i,:) = [];
        B_temp(:,i) = []; % 去掉第i行，第i列
        VV_B = V_RF_i(:,j)' * B_temp * V_RF_i(:,j);
        zeta_B(i,j)= B{j}(i,i) + 2 * real(VV_B);
%         temp = 0;
%         for m = 1:1:N
%             if( m~=i )
%                 for n = 1:1:N
%                     if(n~=i)
%                         temp = temp + conj(V_RF(m,j))*D{j}(m,n)*V_RF(n,j);
%                     end
%                 end
%             end
%         end
        % 计算zeta_D
        D_temp = D{j};
        D_temp(i,:) = [];
        D_temp(:,i) = [];
        VV_D = V_RF_i(:,j)' * D_temp * V_RF_i(:,j);
        zeta_D(i,j)= D{j}(i,i) + 2 * real(VV_D);
        
        % 计算eta_B, eta_D
        V_B = B{j}(i,:) * V_RF(:,j);
        V_D = D{j}(i,:) * V_RF(:,j);
        eta_B(i,j) = V_B - B{j}(i,i)*V_RF(i,j);
        eta_D(i,j) = V_D - D{j}(i,i)*V_RF(i,j);
        
        % 计算theta_1，theta_2
        c(i,j) = (1 + zeta_D(i,j)) * eta_B(i,j) - zeta_B(i,j) * eta_D(i,j);
        z(i,j) = imag( 2 * eta_B(i,j) * eta_D(i,j));
        phi(i,j) = 0;
        tt = asin(imag(c(i,j))/abs(c(i,j)));
        if(real(c(i,j))>=0) phi(i,j) = tt; 
        else phi(i,j) = pi - tt; 
        end
        theta_1 = -phi(i,j) + asin(z(i,j)/abs(c(i,j)));
        theta_2 = pi - phi(i,j) - asin(z(i,j)/abs(c(i,j)));
        
        % 判断最优theta
        V_RF1 = exp(-1j * theta_1);
        V_RF2 = exp(-1j * theta_2);
        trace_A_inv = trace(A_inv);
        f1 = trace_A_inv - ...
            ( zeta_B(i,j) + 2 * real(conj(V_RF1)*eta_B(i,j)))/...
            (1+zeta_D(i,j)+ 2 * real(conj(V_RF1)*eta_D(i,j)));
        f2 = trace_A_inv - ...
            ( zeta_B(i,j) + 2 * real(conj(V_RF2)*eta_B(i,j)))/...
            (1+zeta_D(i,j)+ 2 * real(conj(V_RF2)*eta_D(i,j)));
        if(f1 < f2) theta_opt = theta_1; 
        else theta_opt = theta_2;
        end
        
        V_RF(i,j) = exp(-1j*theta_opt); % 更新V_RF
    end
end

% R = zeros(1,K);
for k = 1:1:K
    R(k) = beta(k) * log2(1+(P(k,k)/sigma2));
end
sum_r = sum(R)