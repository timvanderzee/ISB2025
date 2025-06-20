function [InitGuess] = CreateInitialGuess(bool_guess,params_OS, data_exp, muscles, coeff_LMT_ma)
%UNTITLED4 Summary of this function  here
N      = data_exp.Nspline; 
offset = data_exp.offset;  
x      = data_exp.qspline; 

% vMtilde guess
vMGuess = ones(length(muscles),N)*-0.3;
    
% lMtilde guess
if bool_guess == 0; 
    lMtildaGuess = ones(length(muscles),N)*0.9; % 0.9 of 1.1? 
else
    for i = 1:N
    lMT_ext(i)  = coeff_LMT_ma(1,1)  + coeff_LMT_ma(2,1)*(x(i)+offset(i))  + coeff_LMT_ma(3,1)*(x(i)+offset(i)).^2  + coeff_LMT_ma(4,1)*(x(i)+offset(i)).^3;
    ma_ext(i)   = -coeff_LMT_ma(2,1) + -coeff_LMT_ma(3,1)*(x(i)+offset(i)) + -coeff_LMT_ma(4,1)*(x(i)+offset(i)).^2;
    lMT_flex(i) = coeff_LMT_ma(1,2)  + coeff_LMT_ma(2,2)*(x(i)+offset(i))  + coeff_LMT_ma(3,2)*(x(i)+offset(i)).^2  + coeff_LMT_ma(4,2)*(x(i)+offset(i)).^3;
    ma_flex(i)  = -coeff_LMT_ma(2,2) + -coeff_LMT_ma(3,2)*(x(i)+offset(i)) + -coeff_LMT_ma(4,2)*(x(i)+offset(i)).^2;
    end
    lMtildaGuess_ext = (lMT_ext - params_OS.MT(3,1))./params_OS.MT(2,1);  
    lMtildaGuess_flex = (lMT_flex - params_OS.MT(3,2)) ./params_OS.MT(2,2);
    lMtildaGuess = [lMtildaGuess_ext' lMtildaGuess_flex']; 
end
    
% lM projected guess
params.MTparams = params_OS.MT;
lMo      = params.MTparams(2,:);
alphao   = params.MTparams(4,:)';
lMGuess  = lMtildaGuess.*lMo;
w        = lMo'.*sin(alphao);
lM_projectedGuess = sqrt((lMGuess.^2 - w'.^2));

%dLdMT guess
dlMdtGuess = vMGuess'.*params.MTparams(5,:)./lMo;

% Guesses
InitGuess.vM           = vMGuess; 
InitGuess.lMtilda      = lMtildaGuess; 
InitGuess.lM_projected = lM_projectedGuess; 
InitGuess.dlMdt        = dlMdtGuess; 

end

