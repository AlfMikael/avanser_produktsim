import BeamModels as models
import SolverAlgorithms as algs

num_nodes = 9
model = models.SimplySupportedBeamModel(num_nodes)
#model = models.CantileverWithEndMoment(num_nodes)
# model 3

#algs.solveLinearSteps(model,load_steps=0.01,max_steps=100)
#algs.solveNonlinLoadControl(model,load_steps=0.01,max_steps=100)
#algs.solveArchLength(model,archLength=0.01,max_steps=100)
algs.solveArchLength_copy(model,arcLength=0.01,max_steps=2)
#algs.solveArchLength_test(model,archLength=0.1,max_steps=10)

num_steps = len(model.load_history)

for iStep in range(num_steps):
   print("LoadFactor= {:12.3e}".format(model.load_history[iStep]))
   print("dispVec={:}".format(iStep))
   print(model.disp_history[iStep])

step_inc = (num_steps // 10)
for iStep in range(0,len(model.load_history), step_inc):
    model.plotDispState(iStep)

print("End")
