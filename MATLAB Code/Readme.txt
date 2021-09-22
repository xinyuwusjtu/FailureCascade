This folder contains the MATLAB source files for the work: "Predicting Failure Cascades in Large Scale Power
Systems via the Influence Model Framework". 

1. Main Part: Training the Influence Model to Predict Failure Cascade Sequence

There are 4 code files that need to be run in order as follows:
- "FlowCal_MatPower.m": generate a certain number of failure cascade sequences for training and testing based on flow calculattion over DC and AC models.
- "HybridLearning_IPOPT.m": learn the pairwise influence matrices based on Monte Carlo method, and the weighted matrix based on quadratic optimization (implemented by IPOPT, a software package for large-scale ​nonlinear optimization)
- "ThresholdEstimation.m": learn the threshold values
- "FailureCascadePrediction.m": predict the failure cascade sequences for the testing initial failures which are evaluated under different metrics

Remark: 
- Files "fCascadeNew_Matpower_Revised_Slack.m" and "fCascadeNew_Matpower_Revised_AC_New_Slack.m" are the codes for calculating the failure cascade sequences given initial link failures under DC and AC flow model respectively, 
which are called in file "FlowCal_MatPower.m". 
- Files "runpf_me.m", "rundcpf_me.m", "runacpf_me.m", "bustype_me.m" are the codes for calculating the power flow, which are called in "fCascadeNew_Matpower_Revised_Slack.m" and "fCascadeNew_Matpower_Revised_AC_New_Slack.m". 
All these files are mildly modified from the corresponding original files in IEEE Matpower Toolbox to fit into the failure cascade generation code, including "runpf.m", "rundcpf.m", "bustype.m";
- Files "fun_full_D_weighted_function_vector.m" and "fun_full_D_weighted_gradient_vector.m" are for quadratic optimization in learning the weighted matrix, which are called in "HybridLearning_IPOPT.m".

2. Other Code Files

- "FailureBasicInformation.m": output the basic information of the given power system, including number of buses, generators, transmission lines, power generation and loading conditions, topology, etc.
- "Property_Sparsity.m": study the sparsity of the influence (corresponding to Section V.F in the paper)
- "Critical_Component_Identification.m": identifies the critical components that lead to severe failure cascade based on the trained influence model.
- "Critical_Component_FigureDraw.m": plot the result for "Critical_Component_Identification.m".