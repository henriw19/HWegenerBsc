from pymatching import Matching
import matplotlib.pyplot as plt
from numpy.random import rand
from hwbsc.stab_codes import *
from hwbsc.stabilizers_utils import hypergraph_product
from hwbsc.almlinunionfind import Alm_lin_union_find
import random


def createlinmatrix(n):
    # create linear repetition code of len n
    rep = np.zeros((n-1,n),dtype=int)
    for i in range(n-1):
        for j in range(n):
            rep[i,i] = 1
            rep[i,i+1] = 1
    H = rep
    hz = np.vstack((np.kron(np.eye(H.shape[0],dtype = int),H),np.kron(H.T,np.eye(H.shape[1],dtype = int))))
    hx = np.vstack((np.kron(H,np.eye(H.shape[0],dtype = int)),np.kron(np.eye(H.shape[1],dtype = int),H.T)))


    hyprep = hypergraph_prod(rep).H
    return hyprep, hx, hz

def createringmatrix(n):
    # create ring repetition code of len n
    rep = np.zeros((n,n),dtype=int)
    for i in range(n):
        for j in range(n):
            rep[i,i] = 1
            rep[i,(i+1)%n] = 1
    # print(rep)
    H = rep
    hz = np.vstack((np.kron(np.eye(H.shape[0],dtype = int),H),np.kron(H.T,np.eye(H.shape[1],dtype = int))))
    hx = np.vstack((np.kron(H,np.eye(H.shape[0],dtype = int)),np.kron(np.eye(H.shape[1],dtype = int),H.T)))

    hyprep = hypergraph_prod(rep).H
    # print(hz.T)
    return hyprep, hx, hz

def error_generator(phys_err_rate: float, n: int):
    

    
    xerror = np.zeros(n).astype(int)
    zerror= np.zeros(n).astype(int)
    
    for i in range(n):
        r = random.random() 
        # x error:
        if r < phys_err_rate/3:
            xerror[i] = 1
        # z error:
        elif phys_err_rate/3 < r < 2*phys_err_rate/3:
            xerror[i] = 1
            zerror[i] = 1
        # y error:
        elif 2*phys_err_rate/3 < r < phys_err_rate:
            zerror[i] = 1
        
    return xerror, zerror

#print(error_generator(H))

def log_err_rate(H, hx, hz, logicals, num =20000, p= 0.006):
    # creates matching using pymatching
    n =H.shape[1]//2
    m = H.shape[0]
    

    lx = logicals[0,:n]
    lz = logicals[1,n:]
    xmatching = Matching.from_check_matrix(hz) 
    zmatching = Matching.from_check_matrix(hx)
    # matching = Matching.from_check_matrix(H)
    err_counter = 0
    alt_counter = 0
     
    for _ in range(num):
        # error = error_generator(H,p)
        # actual_flip = logicals @ error % 2
        # syndrome = (H @ error.T) % 2
        # decode_prediction = matching.decode(syndrome).T
        # residual = (error + decode_prediction) % 2
        
        decode_fail = False

        # x error:
        xerror, zerror = error_generator(p, n)
        #print(xerror, zerror)
        # x syndrome:
        #print(hx @ np.array([0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0]))
        xsyndrome = (hz @ xerror.T) % 2
        # print("xerror\n",xerror,"\nxsyndrome\n",xsyndrome)
        #xdecode_prediction = xmatching.decode(xsyndrome).T

        xdecodealmlinunfind = Alm_lin_union_find(hz.T, xsyndrome).decodebin
        xresidual = (xerror + xdecodealmlinunfind) % 2
        
        if np.array_equal(hz @ xresidual % 2, np.zeros(hx.shape[0])) and lz @ xresidual %2 != 0:

            decode_fail = True
                

        # z error:
        
        zerror = np.zeros(hz.shape[1],dtype=int)

        for i in range(hz.shape[1]):
            if rand() < p/2:
                zerror[i] = (zerror[i] +1) % 2
        # z syndrome:
        zsyndrome = (hx @ zerror.T) % 2
        
        # print("zerror\n",zerror,"\nzsyndrome\n",zsyndrome)
        #zdecode_prediction = zmatching.decode(zsyndrome).T
        zdecodealmlinunfind = Alm_lin_union_find(hx.T, zsyndrome).decodebin

        zresidual = (zerror + zdecodealmlinunfind) % 2
        if np.array_equal(hx @ zresidual%2, np.zeros(hz.shape[0])) and (lx @ zresidual) %2 != 0:
         
            decode_fail = True 
      

        # y error:
        # for i in range(2):
        #     if xerror[i] == 1 and xerror[i] == zerror[i]:
        #         pass

        
        

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
        H, hx, hz = createringmatrix(dist)
        logicals = find_log_op(H)
        log_errs = []

        for phys_err_rate in p:
            print(f"per={phys_err_rate}")
            log_errs.append(log_err_rate(H, hx.T, hz.T, logicals, n, phys_err_rate))

        log_err_distances.append(np.array(log_errs))
    
    plt.figure()



    for distance, logical_errors in zip(distances, log_err_distances):
        if text == True:
            data = np.column_stack((p, logical_errors))
            np.savetxt(f"uniondistance{distance}data.txt", data,header='x,y', delimiter=',')
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
        plt.savefig("threshold")

# run threshold plotter from here: output is threshold.png
threshold_plotter([3,5,7,9],np.logspace(-1,-0.71,4),n = 1000, plot = True, text = True)
#np.logspace(-2,-0.7,11)


#print(log_err_rate(creatematrix(3),find_log_op(creatematrix(3)), num = 1, p= 0.3))

