function [lambda,temp_p] = get_lambda(Q_t,P_t,sigma2,K,beta)
lambda = 1/sigma2;
while 1
    initPower = 0;
    posi = 0;
    for k = 1:1:K
        temp_p(k) = (beta(k)/lambda) - Q_t(k,k)*sigma2;
        if(temp_p(k) > 0)
            initPower = initPower + temp_p(k);
            posi = posi + 1;
        end
    end
%             initPower / P_t -1;
%             initPower;
    if( abs(initPower / P_t -1) <= 0.001 )
        disp("found P");
        break;
    end
%             if(posi > 0) lamda = lamda + 0.01*(initPower - P_t)/posi;
%             else lamda = lamda/5;
%             end
    lambda = lambda + 0.001*(initPower - P_t);
    if(lambda <= 0) 
        pause;
    end
end