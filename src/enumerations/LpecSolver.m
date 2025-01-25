classdef LpecSolver
    enumeration
        Gurobi
        Highs % HiGHS solver solved via matlab's intlinprog function 
        Highs_casadi % HiGHS solver solved via casadi's conic function 
        Matlab
        Projected_Gradient  %cf. Kirches 2022 % only if g = []; and G and H subvectors of x
        Reg % use scholtes Reg solver 
        Ell_1 % use scholtes solver 
        Ell_inf % use scholtes solver/costum implement inside lpec not used;  
        Nlp % solve with IPOPT directly
    end   
end

