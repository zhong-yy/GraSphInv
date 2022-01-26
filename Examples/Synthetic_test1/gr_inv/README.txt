Each folder is a workspace for inversion using different parameters and constraints

├── Adap_Inver_wide_bounds     :Inversion using an adaptively refined mesh with an ridiculously wide upper/lower bound
├── Adaptive_Inversion         :an adaptively refined mesh, without a-priori constraints except upper/lower bounds    
├── Adaptive_Inversion_crg     :an adaptively refined mesh, constrained by an a-priori velocity model using cross-gradient coupling
├── Adaptive_Inversion_pet     :an adaptively refined mesh, constrained by an a-priori velocity model using direct-parameter coupling
├── Non-adaptive_Inversion     :a fixed uniform mesh, without a-priori constraints except upper/lower bounds
├── Non-adaptive_Inversion_crg :a fixed uniform mesh, constrained by an a-priori velocity model using cross-gradient coupling
├── Non-adaptive_Inversion_pet :a fixed uniform mesh, constrained by an a-priori velocity model using direct-parameter coupling
└── README.txt                 :This file
