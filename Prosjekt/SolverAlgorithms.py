# example Beam_models# ----------------------------------------------------------------# PURPOSE#  Starting point for a couple of beam models#import mathimport numpy as npimport matplotlib.pyplot as pltimport CorotBeamimport matplotlib.animation as anmfrom copy import deepcopy# ----- Topology -------------------------------------------------def solveArchLength(model, archLength=0.02, max_steps=50, max_iter=30):    """ Initial step, outside step-loop """    num_dofs = model.get_num_dofs()    u_vec = np.zeros(num_dofs)    residual_vec = np.zeros(num_dofs)    d_q_prev = np.zeros(num_dofs)    Lambda = 0    # Performing a single predictor step    #for iStep in range(max_steps):    for iStep in range(0):        res_Vec = np.zeros(num_dofs)        # Getting stiffness matrix        K_mat = model.get_K_sys(uVec)        print("K_mat", K_mat)        Lambda = 0.0        q_Vec = model.get_incremental_load(Lambda)        print("q_vec når lambda 0", q_Vec)        print("Current q-vector: ", q_Vec)        #        d_q = np.linalg.solve(K_mat, q_Vec)        print("d_q", d_q)        uVec = d_q * Lambda        # Unit tangent along equilibrium path        w = np.array([0.5, 0.5])  # for example        f = math.sqrt(1 + w.T @ w)        """Predictor step"""   #     res_Vec = problem.get_residual(uVec,Lambda)  #      q_Vec   = problem.get_vecq(Lambda) #       K_mat   = problem.get_Kmat(uVec)        d_r = np.linalg.solve(K_mat,res_Vec)        d_q = np.linalg.solve(K_mat,q_Vec)        delta_v = 0        """ Corrector iteration """        for iIter in range(1):            # TODO: Implement this            # Finn et uttrykk for delta_v            # regn ut residual            res_Vec = model.get_residual(Lambda, uVec)            K_mat = model.get_K_sys(uVec)            # Displacement som resultat av residual            d_r = np.linalg.solve(K_mat, res_Vec)            # Displacement som resultat av q            d_q = np.linalg.solve(K_mat, q_Vec)            #uVec += delta_uVec            #K =            #d_r = np.linalg.solve(K, res_Vec)            # Hvis residual liten nok -> delta_v OK            if (res_Vec.dot(res_Vec) < 1.0e-15):                break        model.append_solution(Lambda, uVec)        print(" ")def solveArchLength_copy(model, arcLength=0.02, max_steps=50, max_iter=30):    """ Initial step, outside step-loop """    num_dofs = model.get_num_dofs()    uVec = np.zeros(num_dofs)    residual_vec = np.zeros(num_dofs)    d_q_prev = np.zeros(num_dofs)    Lambda = 0.0    # Performing a single predictor step    #for iStep in range(max_steps):    for iStep in range(max_steps):        for iIter in range(max_iter):            res_Vec = model.get_residual(Lambda, uVec)            q_Vec = model.get_incremental_load()            K_mat = model.get_K_sys(uVec)            d_r = np.linalg.solve(K_mat, res_Vec)            d_q = np.linalg.solve(K_mat, q_Vec)            dLambda = 0            if (iIter == 0):  # Predictor step                f = math.sqrt(1.0 + d_q.dot(d_q))                dLambda = arcLength / f                dotV = d_q.dot(d_q_prev)                if (dotV < 0.0):                    dLambda = -dLambda                d_q_prev = dLambda * d_q            else:  # Corrector step                f = (1.0 + d_q.dot(d_q))                dLambda = -(d_r.dot(d_q)) / f            delta_uVec = d_r + dLambda * d_q            uVec = uVec + delta_uVec            Lambda += dLambda            print("i = {:2d}{:3d} dotV={:8.1e} dLambda={:10.1e} Lambda={:12.3e} u =[{:12.3e} {:12.3e}],  r = [{:12.3e} {:12.3e}]".format(                    iStep, iIter, dotV, dLambda, Lambda, uVec[-1], uVec[-2], res_Vec[-1], res_Vec[-2]))            res_Vec = model.get_residual(Lambda, uVec)            if (res_Vec.dot(res_Vec) < 1.0e-15):                break    model.append_solution(Lambda, uVec)    print(" ")def solveArchLength_copy(model, arcLength=0.02, max_steps=50, max_iter=30):    """ Initial step, outside step-loop """    num_dofs = model.get_num_dofs()    uVec = np.zeros(num_dofs)    residual_vec = np.zeros(num_dofs)    d_q_prev = np.zeros(num_dofs)    Lambda = 0.0    # Performing a single predictor step    #for iStep in range(max_steps):    for iStep in range(max_steps):        for iIter in range(max_iter):            res_Vec = model.get_residual(Lambda, uVec)            q_Vec = model.get_incremental_load()            K_mat = model.get_K_sys(uVec)            d_r = np.linalg.solve(K_mat, res_Vec)            d_q = np.linalg.solve(K_mat, q_Vec)            dLambda = 0            if (iIter == 0):  # Predictor step                f = math.sqrt(1.0 + d_q.dot(d_q))                dLambda = arcLength / f                dotV = d_q.dot(d_q_prev)                if (dotV < 0.0):                    dLambda = -dLambda                d_q_prev = dLambda * d_q            else:  # Corrector step                f = (1.0 + d_q.dot(d_q))                dLambda = -(d_r.dot(d_q)) / f            delta_uVec = d_r + dLambda * d_q            uVec = uVec + delta_uVec            Lambda += dLambda            print("i = {:2d}{:3d} dotV={:8.1e} dLambda={:10.1e} Lambda={:12.3e} u =[{:12.3e} {:12.3e}],  r = [{:12.3e} {:12.3e}]".format(                    iStep, iIter, dotV, dLambda, Lambda, uVec[-1], uVec[-2], res_Vec[-1], res_Vec[-2]))            res_Vec = model.get_residual(Lambda, uVec)            if (res_Vec.dot(res_Vec) < 1.0e-15):                break    model.append_solution(Lambda, uVec)    print(" ")def predictor_only(model, archLength=0.02, max_steps=50, max_iter=30):    """ Initial step, outside step-loop """    num_dofs = model.get_num_dofs()    u_vec = np.zeros(num_dofs)    q_vec = model.get_incremental_load()    #residual_vec = np.zeros(num_dofs)    d_q_prev = np.zeros(num_dofs)    Lambda = 0.0    # Performing a single predictor step    residual_vec = model.get_residual(Lambda, u_vec)    K = model.get_K_sys(u_vec)    # Finn displacements fra q    d_u = np.linalg.solve(K, q_vec)    print("q-vektor:", q_vec)    print("d_u:", d_u)    print("Stiffness, last element, K:", K[-1])    for iStep in range(max_steps):        model.append_solution(Lambda, u_vec)def solveNonlinLoadControl(model, load_steps=0.01, max_steps=100, max_iter=30):    num_dofs = model.get_num_dofs()    uVec   = np.zeros(num_dofs)    d_uVec = np.zeros(num_dofs)    Lambda = 0    for iStep in range(max_steps):        Lambda = load_steps * iStep        #TODO: Implement this        for iIter in range(max_iter):            # TODO: Implement this            res_Vec = model.get_residual(Lambda, uVec)            if (res_Vec.dot(res_Vec) < 1.0e-15):                break        model.append_solution(Lambda, uVec)        print(" ")        model.append_solution(Lambda, uVec)        print("Non-Linear load step {:}  load_factor= {:12.3e}".format(iStep, Lambda))def solveLinearSteps(model, load_steps=0.01, max_steps=100):    num_dofs = model.get_num_dofs()    uVec = np.zeros(num_dofs)    for iStep in range(max_steps):        Lambda = load_steps * iStep        #q_Vec = model.get_incremental_load(Lambda)        q_Vec = model.get_incremental_load()        print("q_Vec", q_Vec)        K_mat = model.get_K_sys(uVec)        d_q = np.linalg.solve(K_mat, q_Vec)        print("d_q", d_q)        uVec = d_q * Lambda        model.append_solution(Lambda, uVec)        print("Linear load step {:}  load_factor= {:12.3e}".format(iStep, Lambda))