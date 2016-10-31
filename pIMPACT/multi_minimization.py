'''
require gnu parallel
require objFuncName.py which is the cost function for optimization
'''
import numpy as np
import os

def multi_minimization(objFuncName, bounds, tol=1E99, population=None, popsize=192, ncore=1):
    
    limits = np.array(bounds, dtype="float").T
    center = 0.5 * (limits[0] + limits[1])
    scale = np.fabs(limits[0] - limits[1])
    np.savetxt('generated_center.data',center)
    np.savetxt('generated_scale.data',scale)
    
    if popsize > ncore*100:
        microPopsize = popsize/ncore/10
    else:
        microPopsize = popsize/ncore
    nTask = popsize/microPopsize
    print 'nTask = ',nTask
        
    gen_lines=['import numpy as np',
                 'import sys',
                 'from scipy import optimize',
                 'from '+objFuncName+' import '+objFuncName,
                 'center = np.loadtxt("generated_center.data")',
                 'scale = np.loadtxt("generated_scale.data")',    
                 'population = center + scale*(np.random.rand('+\
                     str(microPopsize)+','+str(len(bounds))+')-0.5)',
                 'result=[]',
                 'for i in range('+str(microPopsize)+'):',
                 '  temp_result = optimize.minimize('+objFuncName+\
                     ', population[i], method = "L-BFGS-B")',
                 '  if temp_result.fun < ' + str(tol) + ':',
                 '    result.append( temp_result.x )',
                 '  if result:',
                 '    np.savetxt("generated_result"+sys.argv[1]+".data", np.array(result))']

    with open('generated_Optim.py','w') as fp:
        fp.write('\n'.join(gen_lines) + '\n')

    arg = [str(i) for i in range(nTask)]
    with open('generated_arg.data','w') as fp:
        fp.write('\n'.join(arg) + '\n')
    os.system('parallel -a generated_arg.data python generated_Optim.py')
    #os.system('parallel python generated_Optim.py ::: '+str(np.arange(nTask))[1:-1])
    
    result=[]
    for i in range(nTask):
        filename = 'generated_result'+str(i)+'.data'
        if os.path.isfile(filename):
            temp = np.loadtxt(filename)
            if temp.ndim > 1:
                for j in range(len(temp)):
                    result.append(temp[j])
            else:
                result.append(temp)
            
    os.system('rm generated_*')
    
    result=np.array(result)
    np.savetxt('pop.data',result)

    return result
