classdef LpecSolver
    enumeration
        Gurobi
        Highs
        Matlab
        Projected_Gradient  %cf. Kirches 2022 % only if g = []; and G and H subvectors of x
        Reg % use scholtes Reg solver 
        Ell_1 % use scholtes solver 
        Ell_inf % use scholtes solver/costum implement inside lpec not used;  
        Nlp % solve with IPOPT directly
    end   
end

