from pymatching import Matching
import matplotlib.pyplot as plt
from numpy.random import rand
from hwbsc.Tri_col_union_find import Tri_alm_lin_union_find
from ldpc.mod2 import row_echelon, rank, nullspace, row_basis
from hwbsc.stabilizers_utils import find_log_op_color
from hwbsc.color_code import tri_color_graph, tri_color_pcm_generator
import numpy as np
import random
import os
import pickle
import itertools

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
        # error = np.array([0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0])
        # error = np.array([0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0])
        # error = np.array([0, 0, , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1])
        errornodes = [] 
        for idx,bin in enumerate(error):
            if bin == 1:
                errornodes.append(idx+facegraph.number_of_nodes())

        syndrome = (H @ error) % 2
        syndrome2 = syndrome.copy()
        
        # print("error",error)
        # print("syndrome",syndrome)
        # print("xerror\n",xerror,"\nxsyndrome\n",xsyndrome)
        #xdecode_prediction = xmatching.decode(xsyndrome).T
        # colG, facegraph, _, _ = color_graph(H)
        decodealmlinunfind = Tri_alm_lin_union_find(colG, facegraph.copy(),syndrome).decode
        # print("hier2",H@np.array(decodealmlinunfind)%2)
        residual = (error + decodealmlinunfind) % 2
        # print(H @ residual % 2)
        # print(logicals @ residual %2)
        if np.array_equal(H @ residual % 2, np.zeros(H.shape[0])) and not np.array_equal(logicals @ residual %2, np.zeros(logicals.shape[0])):
            # print('decode fail')
            decode_fail = True
                
        if decode_fail == True:
        
            alt_counter += 1
        
        if not np.array_equal(H@decodealmlinunfind%2,syndrome2):
            print("FAIL")
            print(error.tolist())
            break
        else:
            pass
            # print("YAY")
            
        





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
        pcm = tri_color_pcm_generator(dist,dist)[0]
        print(dist,"pcm erstellt")
        if not os.path.exists(f"graphs/trigraphdist{dist}"):
            colG, facegraph, node_colors, edge_colors = tri_color_graph(pcm)
            print(f"erstelle {dist}")
            with open(f"graphs/trigraphdist{dist}","wb") as f:
                pickle.dump(colG,f)
            with open(f"graphs/trifacegraphdist{dist}","wb") as f:
                pickle.dump(facegraph,f)
        
        with open(f"graphs/trigraphdist{dist}","rb") as f:
            colG = pickle.load(f)
        with open(f"graphs/trifacegraphdist{dist}","rb") as f:
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
            np.savetxt(f"multithreading/tricoloruniondistance{distance}data.txt", data,header='x,y', delimiter=',')
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

# run threshold plotter from here: output is threshold.png
# print(np.logspace(-1.44,-0.4,7))
#np.logspace(-2,-0.7,11)


#print(log_err_rate(creatematrix(3),find_log_op(creatematrix(3)), num = 1, p= 0.3))


import multiprocessing
# run threshold plotter from here: output is threshold.png
# def task1():
#     threshold_plotter([10],np.logspace(-3,-1.44,7),n = 10000, plot = True, text = True)
#     print("task one is done")

# def task2():
#     threshold_plotter([12],np.logspace(-3,-1.44,7),n = 10000, plot = True, text = True)
#     print("task two is done")

# def task3():
#     threshold_plotter([14],np.logspace(-3,-1.4,7),n = 10000, plot = True, text = True)
#     print("task three is done")

# def task4():
#     threshold_plotter([16],np.logspace(-1.44,-0.4,7),n = 200, plot = True, text = True)
#     print("task four is done")

# process1 = multiprocessing.Process(target=task1)
# process2 = multiprocessing.Process(target=task2)
# process3 = multiprocessing.Process(target=task3)
# # process4 = multiprocessing.Process(target=task4)

# process1.start()
# process2.start()
# process3.start()
# # process4.start()

# process1.join()
# process2.join()
# process3.join()
# # process4.join()


import threading

def process_item(item):
    dist = 12
    pcm = tri_color_pcm_generator(dist,dist)[0]
    print(dist,"pcm erstellt")
    if not os.path.exists(f"graphs/trigraphdist{dist}"):
        colG, facegraph, node_colors, edge_colors = tri_color_graph(pcm)
        print(f"erstelle {dist}")
        with open(f"graphs/trigraphdist{dist}","wb") as f:
            pickle.dump(colG,f)
        with open(f"graphs/trifacegraphdist{dist}","wb") as f:
            pickle.dump(facegraph,f)
    with open(f"graphs/trigraphdist{dist}","rb") as f:
        colG = pickle.load(f)
    with open(f"graphs/trifacegraphdist{dist}","rb") as f:
        facegraph = pickle.load(f)
    logicals = find_log_op_color(pcm)
    bla = log_err_rate(pcm,colG,facegraph, logicals, 20000, item)  
    print(item,bla)  

items_to_process = np.logspace(-3,-1.44,9)

def thread_worker(items):
    for item in items:
        process_item(item)

# Determine the number of threads to use
num_threads = 9  # Adjust as needed

# Split the list of items into equal portions for each thread
split_items = [items_to_process[i::num_threads] for i in range(num_threads)]

# Create and start the threads
threads = []
for i in range(num_threads):
    thread = threading.Thread(target=thread_worker, args=(split_items[i],))
    threads.append(thread)
    thread.start()

# Wait for all threads to complete
for thread in threads:
    thread.join()

