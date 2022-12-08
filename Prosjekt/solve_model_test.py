import BeamModels as models
import SolverAlgorithms as algs

num_nodes = 9
#model = models.SimplySupportedBeamModel(num_nodes)
model = models.CantileverWithEndMoment(num_nodes)
#model = models.SimplySupportedBeamModel(num_nodes)
# model 3

#algs.solveLinearSteps(model,load_steps=0.01,max_steps=100)
#algs.solveNonlinLoadControl(model,load_steps=0.01,max_steps=100)
#algs.solveArchLength(model,archLength=0.01,max_steps=100)
algs.solveArchLength_test(model,archLength=0.01,max_steps=100)
