<<<<<<< HEAD
=======
from pymatching import Matching
import matplotlib.pyplot as plt
from numpy.random import rand
from hwbsc.Col_union_find import Alm_lin_union_find
from ldpc.mod2 import row_echelon, rank, nullspace, row_basis
from hwbsc.stabilizers_utils import find_log_op_color
from hwbsc.color_code import hexa_color_pcm_generator, hexa_color_graph
import numpy as np
import random
import pickle
import os

def error_generator(phys_err_rate: float, n: int):
    
    error= np.zeros(n).astype(int)
    
    for i in range(n):
        r = random.random() 
        # error:
        if r < phys_err_rate:
            error[i] = 1
    return error

#print(error_generator(H))

def log_err_rate(H, colG, facegraph, logicals, num =20000, p= 0.006):
    # creates matching using pymatching
    n =H.shape[1]
    m = H.shape[0]
    
    err_counter = 0
    alt_counter = 0
    for y in range(num):
        if y%1000 == 0:
            print(y)
        decode_fail = False
        # x error:
        error = error_generator(p, n)
        
        errornodes = []
        for idx,bin in enumerate(error):
            if bin == 1:
                errornodes.append(idx+facegraph.number_of_nodes())

        syndrome = (H @ error) % 2
        
        # print("error",error)
        # print("syndrome",syndrome)
        # print("xerror\n",xerror,"\nxsyndrome\n",xsyndrome)
        #xdecode_prediction = xmatching.decode(xsyndrome).T
        # colG, facegraph, _, _ = color_graph(H)
        decodealmlinunfind = Alm_lin_union_find(colG, facegraph.copy(),syndrome).decode
        # print("hier2",H@np.array(decodealmlinunfind)%2)
        residual = (error + decodealmlinunfind) % 2
        # print(H @ residual % 2)
        # print(logicals @ residual %2)
        if np.array_equal(H @ residual % 2, np.zeros(H.shape[0])) and not np.array_equal(logicals @ residual %2, np.zeros(logicals.shape[0])):
            # print('decode fail')
            decode_fail = True
                
        if decode_fail == True:
           
            alt_counter += 1
        
        





        #print(check_commutativity(logicals[1],xresidual))
        # if np.array_equal(hx @ xresidual%2, np.zeros(hx.shape[0])) and check_commutativity(logicals[1],xresidual) != 0:
        #     alt_counter += 1/2

        # if np.array_equal(hz @ zresidual%2, np.zeros(hz.shape[0])) and check_commutativity(logicals[0],zresidual) != 0:
        #     alt_counter += 1/2
            
        # predicted_flips = (logicals @ decode_prediction %2)
         
        # if not np.array_equal(actual_flip, predicted_flips):
        #     err_counter += 1
    return  alt_counter/num #err_counter/num 

def threshold_plotter(distances:list, p: list, n:int, text:bool= False, plot:bool = True):
    np.random.seed(2)

    log_err_distances = []

    for dist in distances:
        print(f"d = {dist}")
        pcm = hexa_color_pcm_generator(dist,int(dist-6))[0]
        
        
        

        if not os.path.exists(f"graphs/graphdist{dist}"):
            colG, facegraph, node_colors, edge_colors = hexa_color_graph(pcm)
            print(f"erstelle {dist}")
            with open(f"graphs/graphdist{dist}","wb") as f:
                pickle.dump(colG,f)
            with open(f"graphs/facegraphdist{dist}","wb") as f:
                pickle.dump(facegraph,f)
        
        with open(f"graphs/graphdist{dist}","rb") as f:
            colG = pickle.load(f)
        with open(f"graphs/facegraphdist{dist}","rb") as f:
            facegraph = pickle.load(f)



        logicals = find_log_op_color(pcm)
        log_errs = []

        for phys_err_rate in p:
            print(f"per={phys_err_rate}")
            log_errs.append(log_err_rate(pcm,colG,facegraph, logicals, n, phys_err_rate))

        log_err_distances.append(np.array(log_errs))
    
    plt.figure()



    for distance, logical_errors in zip(distances, log_err_distances):
        if text == True:
            data = np.column_stack((p, logical_errors))
            np.savetxt(f"thresholddata/coloruniondistance{distance}data.txt", data,header='x,y', delimiter=',')
            print("distance",dist,data.tolist())
        elif plot == True:
            plt.plot(p, logical_errors, label = f"distance {distance}")
            plt.grid(True)

        # std_err = (logical_errors *(1-logical_errors)/n)**0.5
        # plt.errorbar(p, logical_errors, yerr = std_err, label = f"distance {dist}")
    
    if plot == True:
        plt.xscale('log') 
        plt.yscale('log')
        #plt.plot(0.1629,0.3,'ro')
        plt.xlabel("phys err rate")
        plt.ylabel("logical err rate")
        plt.legend(loc = 0)
        plt.savefig("col threshold")


import multiprocessing
# run threshold plotter from here: output is threshold.png
def task1():
    threshold_plotter([36],np.logspace(-3,-1.44,9),n = 10000, plot = True, text = True)
    print("task one is done")

def task2():
    threshold_plotter([48],np.logspace(-3,-1.44,9),n = 10000, plot = True, text = True)
    print("task two is done")

def task3():
    threshold_plotter([36],np.logspace(-1.44,-0.4,7),n = 500, plot = True, text = True)
    print("task three is done")

def task4():
    threshold_plotter([48],np.logspace(-1.44,-0.4,7),n = 500, plot = True, text = True)
    print("task four is done")

process1 = multiprocessing.Process(target=task1)
process2 = multiprocessing.Process(target=task2)
# process3 = multiprocessing.Process(target=task3)
# process4 = multiprocessing.Process(target=task4)

process1.start()
process2.start()
# process3.start()
# process4.start()

process1.join()
process2.join()
# process3.join()
# process4.join()


>>>>>>> 6cb60cdd1f8d65da7aed06ab833d7f64cbf0d911
