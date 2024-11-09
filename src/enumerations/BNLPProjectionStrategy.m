
classdef BNLPProjectionStrategy
    enumeration
        LPEC  % solve lpec
        LPECFlipped  % solve lpec with flipped sign
        LPECRandom % solve Lpec with random objective (-1,1)
        Simple % xk1>=xk2
        ActiveSetFunction % use fancy function for active set (todo)
    end
   
end

