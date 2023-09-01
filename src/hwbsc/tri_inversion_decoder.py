import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import itertools
import galois
import sympy as sp
from ldpc.mod2 import row_echelon, rank, nullspace, row_basis
# from hwbsc.codes import hamming_code
from scipy.linalg import solve
from numpy.linalg import matrix_rank
import cProfile
from ldpc2.gf2sparse import LuDecomposition





def create_linear_comb_with_zero(vec1: np.ndarray, vec2: np.ndarray, index: int) -> np.ndarray:
    """ Creates a linear combination of vec1 and vec2;  such that the component at the 
    given index is zero by subtracting a multiple of vec1 from vec2. """
    assert vec1.ndim == 1
    assert vec1.shape == vec2.shape
    ### BEGIN SOLUTION
    new1 = vec2 - (vec2[index] / vec1[index]) * vec1
    # new2 = vec4 + (vec4[index] / vec3[index]) * vec3
    return new1, -vec2[index] / vec1[index]
    ### END SOLUTION

def pivot_matrix_inplace(ext_coeff_mat: np.ndarray, column: int,unit):
    """ Pivots extended coefficient matrix such that the diagonal element in
    the given column has a larger absolute value than any element below it.
    If column >= ext_coeff_mat.shape[1] - 1; do not change matrix"""
    
    
    ### BEGIN SOLUTION
    
    if column >= ext_coeff_mat.shape[1] - 1:
        
        return 
    # find argmax of absolute value for all entries below diagonal elem
    # add colum to it to get index in the colum
    pivot_pos = column + np.argmax(np.abs(ext_coeff_mat[column : , column]))
    # unit_pivot_pos = column + np.argmax(np.abs(unit[column : , column]))
    
    # swap rows
    temp = np.copy(ext_coeff_mat[column])
    ext_coeff_mat[column] = np.copy(ext_coeff_mat[pivot_pos])
    ext_coeff_mat[pivot_pos] = temp

    temp2 = np.copy(unit[pivot_pos])
    unit[pivot_pos] = unit[column]
    unit[column] = temp2
    

    # actions.append((column, pivot_pos))

    # temp2 = np.copy(unit[column])
    # unit[column] = np.copy(unit[pivot_pos])
    # unit[pivot_pos] = temp2

    # ext_coeff_mat[[column, pivot_pos]] = ext_coeff_mat[[pivot_pos, column]]


def upper_triangle(ext_coeff_mat: np.ndarray, pivoting: bool, actions):
    unit = np.eye(ext_coeff_mat.shape[0])
    actions = []
    for colum_idx in range(ext_coeff_mat.shape[1] ):
        for i, row in enumerate(ext_coeff_mat[colum_idx + 1 : ], colum_idx):
            
            # row_ref is the row we use to form linear combinations with
            # this has to be the row with the same index as the colum we are 
            # working with
            row_ref = ext_coeff_mat[colum_idx]
            if pivoting:
                pivot_matrix_inplace(ext_coeff_mat, colum_idx, unit)
            
            # ##########print(unit)
            # ##########print(colum_idx)
            new_row, factor = create_linear_comb_with_zero(row_ref, row, colum_idx)
            # ##########print(new_row, new_unit_row)
            
            if not np.array_equal(ext_coeff_mat[i+1], new_row):
                unit[i+1] += -factor * unit[colum_idx]
              
                actions.append((i+1,colum_idx, factor))
            
            ext_coeff_mat[i + 1] = new_row  
              
            # ##########print(unit,"\n")
    

    
    return ext_coeff_mat,unit

# upper = upper_triangle(ldpc,True,[]) #,"\n",upper_triangle(A_test,True)[1])



# ##########print(upper[0],"\n", upper[1])
A_test = np.array([[1,2,-3,1],
                   [2,4,0,7],
                   [-1,3,2,0]]).astype(float)

hamming = np.array([[0,1,1,0,1,1,0],
                       [0,0,1,1,1,0,1],
                       [1,1,1,1,0,0,0]])

# ldpc = GF([[1,0,0,0,0,0,1],
#             [1,1,1,0,0,0,0],
#             [0,0,0,1,0,0,0],
#             [0,0,1,1,1,0,0],
#             [1,0,0,0,1,1,1]])

# error = GF([0,0,1,0,0,1,1])
# syndrome= GF([1, 1, 0, 1, 0])
# ##########print("syndrome:",ldpc @ error)


def subpcm(colG, subcluster, facegraph):
    inside = []
    translbtwinsidenodesandpcm = {}

    complement = set(facegraph.nodes()) - set(subcluster)

    for l,node in enumerate(colG):
        # if (1 < len(list(set(colG.neighbors(node)).intersection(subcluster))) < 4 and len(list(colG.neighbors(node))) < 4) or (len(list(set(colG.neighbors(node)).intersection(subcluster))) == 1 and len(list(colG.neighbors(node))) ==1):
        if len(list(set(colG.neighbors(node)).intersection(complement))) == 0 and not node in facegraph.nodes():
            inside.append(node)
            translbtwinsidenodesandpcm[node] = len(inside)-1

    # suberror = []
    # if len(err) > 1:
    #     for node in list(colG.nodes()):
    #         if node in inside:
    #             suberror.append(err[int(node-facegraph.number_of_nodes())])
    translbtwancillasandpcm = {}
    for k,ancilla in enumerate(subcluster):
        translbtwancillasandpcm[ancilla] = k
    subpcmm = np.zeros((len(subcluster),len(inside))).astype(int)
 
    for i,insidenode in enumerate(inside):
        for j,neighbor in enumerate(list(colG.neighbors(insidenode))):
            if neighbor in subcluster:
                subpcmm[translbtwancillasandpcm[neighbor],translbtwinsidenodesandpcm[insidenode]] = 1
    return subpcmm, translbtwancillasandpcm, translbtwinsidenodesandpcm, inside

def findaustauschbar(U):
    austauschen = {}
    delcols = []
    for colidx,col in enumerate(U.T):
        if sum(np.array(col)) == 0:
            delcols.append(colidx)
        else:
            lastzero = np.nonzero(col)[-1][-1]
            if lastzero not in austauschen:
                austauschen[lastzero]= {colidx}
            else:
                austauschen[lastzero]|={colidx}
    austauschbar = []
    for col,tauschis in austauschen.items():
        if len(tauschis) > 1:
            austauschbar.append(list(tauschis))
    #####print("delcols",austauschbar,delcols)
    return austauschbar, delcols

def solver(merged,U,y,Ushape1, combinations,delcols,subpcmm, inside):
    lowestweight = U.shape[1]
    subsolution =0
    length = len(combinations)
    ###print(length)
    c=0
    GF = galois.GF(2)
    rang = 0
    y = list(y)
    U = np.array(U)

    for i,delete in enumerate(combinations):
        combinations[i] = set(delete)
    if delcols:
        merged += delcols
    merged = set(merged)
    for delete in combinations[100000:]:
        mergednow = merged.copy()
        # if delcols:
        #     mergednow += delcols
        # mergednow = list(mergednow)
        # mergednow = [x for x in mergednow if x not in delete]
        mergednow = mergednow - delete


        # mergednow = sorted(mergednow)  
        ynow = y.copy()
        Unow = U.copy()
        
        for index in mergednow:
            addline = np.zeros(Ushape1)
            addline[index] = 1
            Unow = np.vstack((Unow,addline)).astype(int)
            ynow.append(0)
        # Unow = np.array(Unow).astype(int)
        # ynow = np.array(ynow).astype(int)
        # Unow = GF(Unow)
        # ynow = GF(ynow)
      
            
        if Unow.shape[0] == Unow.shape[1]:
            # solut = np.linalg.solve(Unow,ynow)
            solut = solve(Unow,ynow)
            
            solut = np.abs(np.array(solut).astype(int)) %2
            predsyn = np.array(GF(subpcmm) @ GF(solut))
            if np.array_equal(predsyn,np.array(insidesubsyndrome)):
                subsolution = solut
    
                ###print('cluster successfully solved')
                # weight = sum(np.array(solut))
                # if weight < lowestweight:
                #     lowestweight = weight
                #     subsolution = solut
                break
                # solut = list(np.array(solut).astype(int))
            else:
                c += 1
                #####print("ERROR")
                if c%10 == 0:
                    ###print(c/length)
                    pass
        if c >= len(combinations)-1:
            pass
            ####print("decode FAIL")
            fullsolution = np.zeros(facegraph.number_of_nodes()*2).astype(int)
        else:
            fullsolution = np.zeros(facegraph.number_of_nodes()*2).astype(int)
            for node in range(len(subsolution)):
                if subsolution[node] == 1:
                    fullsolution[inside[node]-facegraph.number_of_nodes()] = 1

        errnodes = []
        for idx,i in enumerate(fullsolution):
            if i == 1:
                errnodes.append(idx)
    return subsolution


def check_solvable(colG, syndrome, facegraph, subcluster, letztercluster):
    subcluster =sorted(list(subcluster))
    subpcmm, translbtwancillasandpcm, translbtwinsidenodesandpcm, inside = subpcm(colG,subcluster, facegraph)
    insidesubsyndrome = []
    for num in subcluster:
        insidesubsyndrome.append(int(syndrome[num]))
    plud = LuDecomposition(subpcmm, lower_triangular=True)
    subsolution = plud.solve(insidesubsyndrome)


    #print("letzt",letztercluster)
    insidesubsyndromeletzte = []
    for num0 in letztercluster:
        insidesubsyndromeletzte.append(int(syndrome[num0]))
   
    for i in subcluster:
        if i in letztercluster and syndrome[i] == 1:
            insidesubsyndromeletzte[i] = 0 
    if letztercluster:
        #print("jaaa letztercluster",letztercluster)
        subletzte = subpcm(colG,letztercluster, facegraph)[0]
        #print("vorletzteiteration PCM:",subletzte)
        #print("vorletzteiteration syndrome:",insidesubsyndromeletzte)
        
        plud0 = LuDecomposition(subletzte, lower_triangular=True)
        subsolution0 = plud0.solve(insidesubsyndromeletzte)

    # check if cluster is solvable, return False if not

    if letztercluster:
        if np.array_equal(np.array(insidesubsyndrome), np.array(subpcmm) @ np.array(subsolution)%2) and np.array_equal(np.array(insidesubsyndromeletzte), np.array(subletzte) @ np.array(subsolution0)%2):
            #print("beide lösbar")
            return True
        else:
            #print("subcluster not solvable")
            return False
    elif np.array_equal(np.array(insidesubsyndrome), np.array(subpcmm) @ np.array(subsolution)%2):
        return True
    else:
        #print("nicht lösbar")
        return False

def mattosym(numeric_matrix,i,rechts = []):
    # Get the dimensions of the matrix
    num_rows = len(numeric_matrix)
    num_cols = len(numeric_matrix[0])

    # Create the symbolic variables for the unknowns
    variables = sp.symbols(f'{i}0:%d' % num_cols)

    # Create a list to store the symbolic equations
    symbolic_equations = []

    # Iterate over each row of the numeric matrix
    for i in range(num_rows):
        # Build the left-hand side of the equation
        lhs = sum(variables[j] * numeric_matrix[i][j] for j in range(num_cols))

        # Build the right-hand side of the equation (assume it is 0, change if needed)
        if len(rechts) > 0:
            rhs = rechts[i]
        else:
            rhs = 0

        # Create the symbolic equation and add it to the list
        equation = sp.Eq(lhs, rhs)
        symbolic_equations.append(equation)
    return symbolic_equations, variables

def lines_not_in_matrix(matrix_a, matrix_b):
    # Convert the columns of matrix_a and matrix_b to sets
    matrix_a = matrix_a.T
    matrix_b = matrix_b.T
    set_a = set(tuple(col) for col in matrix_a.T)
    set_b = set(tuple(col) for col in matrix_b.T)

    # Find the columns in matrix_a that are not in matrix_b
    difference = set_a.difference(set_b)

    # Convert the resulting set back to NumPy arrays
    indices = []
    for col in difference:
        indices.append(np.argwhere(np.all(matrix_a.T == col, axis=1)))
    return np.sort(np.concatenate(indices).ravel())

def gfinversion(colG, subcluster, facegraph, syndrome): 
    
    subcluster =sorted(list(subcluster))
    subpcmm, translbtwancillasandpcm, translbtwinsidenodesandpcm, inside = subpcm(colG,subcluster, facegraph)
    #print(subcluster)
    #print(inside)
    #print(subpcmm)
    insidesubsyndrome = []
    for num in subcluster:
        insidesubsyndrome.append(int(syndrome[num]))

    plud = LuDecomposition(subpcmm, lower_triangular=True)
    subsolution = plud.solve(insidesubsyndrome)
    
    # check if cluster is solvable, return False if not
    if not np.array_equal(np.array(insidesubsyndrome), np.array(subpcmm) @ np.array(subsolution)%2):
        #print("subcluster not solvable")
        return False

    ##print("subs",list(subsolution))
    ###print(len(subsolution),subsolution)
    ###print(subpcmm@subsolution%2)
    ###print(np.array(insidesubsyndrome))
    #print("subcluster solved")
    fullsolution = np.zeros(int(len(list(colG.nodes()))-len(syndrome))).astype(int)
    for node in range(len(subsolution)):
        if subsolution[node] == 1:
            fullsolution[inside[node]-facegraph.number_of_nodes()] = 1
    
    subchecknodes = []

    for idx,ii in enumerate(insidesubsyndrome):
        if ii == 1:
            subchecknodes.append(subcluster[idx])

    subsolutionnodes = []
    for idx2,iii in enumerate(subsolution):
        if iii == 1:
            subsolutionnodes.append(inside[idx2])


    # ##print("suberrornodes",suberrornodes)
    # ##print("hier",subchecknodes)
    # ##print(subsolutionnodes)
    
    # ##print("suberror",np.array(suberror))
    # ##print("1",np.array(insidesubsyndrome))
    # ##print("2",subpcmm @ suberror %2)
    if np.array_equal(subpcmm @ subsolution%2, insidesubsyndrome):
        #print("Success")  
        pass 
    else:
        #print("Failure")
        pass
    
    
    
    return fullsolution





#     GF = galois.GF(2)
#     ###print("test",subpcmm@suberror%2)    
#     ###print("insidesub\n",np.array(insidesubsyndrome))
#     ####print("subcluster",subclustr)

#     subsym, variables = mattosym(subpcmm,'x',insidesubsyndrome)
#     solution_space = list(sp.solve(subsym, variables, check=False))
#     ###print("subsym",solution_space)

#     subgraph = colG.subgraph(inside)
    
#     GF = galois.GF(2)
#     insidesubsyndrome = GF(insidesubsyndrome)
#     A = GF(subpcmm) ; A
#     P, L, U = A.plu_decompose()
#     # U = row_echelon(subpcmm)[0]
#     y = np.linalg.solve(L, P.T@insidesubsyndrome)
#     U = np.array(U)

#     # eliminate redundant rows:
#     rowU = np.array(row_basis(U))
#     deleted_row = lines_not_in_matrix(U,rowU)
#     ##print(np.array_equal(rowU,np.delete(U,deleted_row,axis=0)))
#     ##print(deleted_row)
#     ##print(rowU)

#     # eliminate redundant columns:
#     ##print(U)
#     ##print(U.shape)
#     rowcolU = row_basis(rowU.T).T
#     ##print(rowcolU)
#     ##print(rowcolU.shape)
#     deleted_col = lines_not_in_matrix(rowU.T,rowcolU.T)
#     ##print(deleted_col)
#     ##print(np.array_equal(rowcolU,np.delete(rowU,deleted_col,axis=1)))

    
#     newsyn = np.array(np.delete(insidesubsyndrome,deleted_row))
#     # newH = np.delete(subpcmm,deleted_row,axis=0)
#     newH = np.delete(subpcmm,deleted_col,axis=1)

#     solu = np.linalg.solve(rowcolU,newsyn).astype(int)%2
#     ##print("solu",newH @solu%2)
#     ##print(newsyn)
#     fullsolu = solu
#     for ind,num in enumerate(deleted_col):
#         fullsolu = fullsolu.copy()
#         fullsolu = np.insert(fullsolu,num,0)
#     ##print(subpcmm@fullsolu%2)
#     ##print(insidesubsyndrome)



#     for i,col in enumerate(U.T):
#         if not np.any(np.all(colU == col, axis=0),axis=1):
#             ##print(i)
#         else:
#             ##print(col)

    
#     for i,col in enumerate(U.T):
#         if col not in row_basis(U.T):
#             ##print(i)
#     ##print(np.delete(U,[9,14,15,17,20,21,22],axis = 1))
#     # ##print(row_basis(U.T).T)
#     ##print(row_basis(U.T).T.shape, U.shape)
   
    
#     ding = np.linalg.solve(row_basis(row_basis(U.T).T), y[:-2])
#     ##print(ding)
#     ding2= [0 ,0, 0, 0, 1 ,0 ,1 ,1 ,0 ,0 ,0 ,1 ,1 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,0 ,0 ,0, 1]
#     ##print(subpcmm@ding2%2)
#     ##print(insidesubsyndrome)
#     ##print("newinside",np.insert(ding,[9,14,15,17,20,21,22],0))
        




    
#     blaa1, variablen2 = mattosym(U,'x',y)
#     conditionlist = list(list(sp.solve(blaa1, variablen2,check = False,  manual = True))[0])
#     ###print(conditionlist)
#     # for condition in conditionlist:
#     #     ###print(condition[-1])

#     rel_gleichungen = [str(gleichung) for gleichung in conditionlist if (len(str(gleichung)) > 2)\
#                         and not ((str(gleichung)[-2] == 'x') \
#                        or (str(gleichung)[-3] == 'x')) and not (int(str(gleichung)[-1]) % 2 ==0 and \
#                                                              str(gleichung)[-2]!= '/')]
    
#     # rel_gleichungen2 = rel_gleichungen
#     for i,gleichung in enumerate(rel_gleichungen):
#         # ###print(sp.Eq(gleichung,0, evaluate = False))
#         gleichung2 = gleichung
#         # .replace('/2','')
#         splitgl = gleichung.split()
#         for element in splitgl:
#             if '2*' in element or '4*' in element:
#                 gleichung2 = gleichung2.replace(element, "")
#             if '/2' in element:
#                 for element2 in splitgl:
#                     if not '/2' in element2:
#                         gleichung2 = gleichung2.replace(element2, "")
#                 gleichung2 = gleichung2.replace('/2', "")
#         rel_gleichungen[i] = gleichung2
    

#     # ###print(rel_gleichungen)
#     def maker(rel_gleichungen):

#         solucion = dict()
#         solucion['+'] = 0
#         solucion['-'] = 0
#         for gleichung in rel_gleichungen:
        
#             global_break = False
#             counter = 1
#             for element in gleichung.split():
#                 if element in solucion:
#                     counter += 1
#             if counter % 2 == 0:
#                 break 
#             for element in gleichung.split():
#                 if element not in solucion:
#                     solucion[element] = 1
#                     break

            
#         return solucion
#     # ###print(maker(rel_gleichungen))
#     # for gleichung in conditionlist:
#     #     gleichung = list(gleichung)
#         # for in gleichung:
#         #     ###print(thing)
#             # if type(thing) is int:
#             #     for variable in gleichung:
#             #         if variable.notset():
#             #             variable = -thing
#             #             break to naechstegleichung






    
# # # Delete the free variables of the subpcm
# #     free_variables = set()
# #     for row in range(U.shape[0]):
# #         for col in range(U.shape[1]):
# #             if row == col and U[row,col] != 1 or sum(np.array(U[row])) == 0:
# #                 free_variables.add(row)
# #     free_variables = list(free_variables)
# #     # U = np.delete(U, free_variables,axis=0)
# #     # y = np.delete(y, free_variables)
# #     ###print("test15",(np.array(U))@np.array(suberror)%2)
# #     ###print(check_solvable(subpcmm,insidesubsyndrome))
# #     # ###print("free_variables",free_variables)
# #     ###print("y",y)
# #     ###print("U",U)
    
# #     # syndrome_equations, variables = mattosym((np.array(P))@np.array(subpcmm),'x',P@insidesubsyndrome)
# #     symbolic_equations, variables = mattosym(U,'x',y)
# #     # ###print(variables)
# #     # ###print("syndrome euqaions",syndrome_equations)
# #     # # for i,equation in enumerate(symbolic_equations):
# #     # #     symbolic_equations[i] = sp.Eq(symbolic_equations[i].lhs,y[i])
# #     # ###print("symbolic_equations",symbolic_equations)
# #     # U = sp.Matrix(U)
# #     # y = sp.Matrix(y)
    
# #     # all_equations = set(syndrome_equations) | set(symbolic_equations)
# #     # ###print("all_equations",list(all_equations))
# #     # allmatrix = sp.linear_eq_to_matrix(list(all_equations), variables)[0]
# #     # ###print(allmatrix.shape,rank(np.array(allmatrix)))
    
# #     # solution_space = sp.solve(list(all_equations), variables, check=False, manual = True)  #,exclude=['x0','x2','x7','x32','x33','x50'])
# #     # ###print(solution_space)

# #     U = sp.Matrix(U)
# #     y = sp.Matrix(y)




# #     # actual_error = solution_space.subs({'x0': 1, 'x2': 1, 'x7':1,'x32':1,'x33':1,'x50':1,'x1':0,'x3':0,'x4':0,'x5':0,'x6':0,'x11':0,'x8':0,'x9':0,'x10':0,'x12':0,'x13':0,'x14':0,'x15':0,'x16':0,'x17':0,'x18':0,'x19':0,'x20':0,'x21':0,'x22':0,'x23':0,'x24':0})
# #     # ###print("actual_error",actual_error)
# #     solution_space = sp.linsolve(symbolic_equations, variables)
# #     loesungsraum = []
# #     # ###print("solution_space",solution_space)


# #     # finde Konstanten und definiere Lösungsraum als liste an Equations
# #     # constants = [0]*len(list(solution_space)[0])
# #     # for j,solution_eq in enumerate(list(solution_space)[0]):
# #     #     numbers = str((solution_space[0][j].args))[1]
# #     #     if len(str(solution_space[0][j])) > 3:
# #     #         if str(solution_space[0][j])[-2] == ' ' or str(solution_space[0][j])[-3] == ' ' and str(solution_space[0][j])[-2] != 'x':
# #     #             constants[j] = int(str(solution_space[0][j])[-1])%2
# #     #     loesungsraum.append(sp.Eq(solution_space[0][j],0))
# #     # ###print("constants",constants)
    
# #     # constants = [0,1,1,0,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,1,0,0]
# #     loeschen = []
# #     for i,equation in enumerate(loesungsraum):
# #         loesungsraum[i] = sp.Eq(loesungsraum[i].lhs - variables[i], loesungsraum[i].rhs)
# #         if loesungsraum[i] == True:
# #             loeschen.append(i)


# #     # ###print("loesungsraum",loesungsraum)
# #     U = np.array(U)
# #     syndrome_equations, variables = mattosym((np.array(P))@np.array(subpcmm),'x',P@insidesubsyndrome)

# #     symbolic_equations, variables = mattosym(U,'x',y)
# #     # ###print(variables)
# #     # ###print("syndrome euqaions",syndrome_equations)
# #     # for i,equation in enumerate(symbolic_equations):
# #     #     symbolic_equations[i] = sp.Eq(symbolic_equations[i].lhs,y[i])
# #     # ###print("symbolic_equations",symbolic_equations)
# #     U = sp.Matrix(U)
# #     y = sp.Matrix(y)
    
# #     all_equations = set(syndrome_equations) | set(loesungsraum)
# #     # ###print("all_equations",list(all_equations))
# #     # allmatrix = sp.linear_eq_to_matrix(list(all_equations), variables)[0]
# #     # ###print(allmatrix.shape,rank(np.array(allmatrix)))
# #     solution_space = sp.solve(list(all_equations), variables, check=False, manual = True)  #,exclude=['x0','x2','x7','x32','x33','x50'])
# #     # ###print("sol space",solution_space)











# #     # ###print("loesungsraum",loesungsraum)
# #     # loesungsraum = np.delete(loesungsraum,loeschen)
# #     # maybe = np.array(list(sp.linsolve(loesungsraum, variables))[0])%2
# #     # ###print(np.array(insidesubsyndrome)%2)
# #     # ###print("hier?",subpcmm@maybe%2)
    
    
# #     # erstelle Matrixfrom des loesungsraums
# #     # coeff_matrix = sp.linear_eq_to_matrix(solution_space, variables)[0]
# #     # coeff_matrix = np.array(coeff_matrix).astype(int)%2
# #     # ###print("coeff_matrix\n", coeff_matrix)
# #     # # finde kleinste weight lösung im fehlerraum
# #     # smallest = coeff_matrix.shape[0]
# #     # for k,cols in enumerate(coeff_matrix.T):
# #     #     weight = sum(cols)
# #     #     errorrr = cols
# #     #     subsolution = (errorrr + np.array(constants))%2
# #     #     # ###print(np.array(P)@np.array(subpcmm)@subsolution%2,np.array(P)@np.array(insidesubsyndrome)%2)
# #     #     if np.array_equal(subpcmm@subsolution%2,np.array(P)@np.array(insidesubsyndrome)%2):# weight < smallest and weight >= 1:
# #     #         ###print("SUBERFOLG")
# #     #         break
# #     #         smallest = weight
    



# #     # ###print("err",errorrr)
# #     # subsolution = (errorrr + np.array(constants))%2
    
# #     # ###print("test2",subpcmm@subsolution%2,"\n","syn",insidesubsyndrome,"\n")


    
# #     # Assuming you already have the eigenvalues in a list

# #     # Assuming you have the matrix A (replace this with your actual matrix)
# #     # A = coeff_matrix
# #     # eigenvalues = [1,1,1]
# #     # Iterate over each eigenvalue and find the corresponding eigenvectors
# #     # eigenvectors = []
# #     # for eigenvalue in eigenvalues:
# #     #     # Form the matrix (A - λI)
# #     #     eigen_matrix = A - eigenvalue * np.eye(A.shape[0])

# #     #     # Find the null space (kernel) of the matrix
# #     #     # _, null_space = np.linalg.qr(eigen_matrix.T, mode='complete')
# #     #     null_spacee = nullspace(eigen_matrix)
# #     #     # Append the null space vectors to the eigenvectors list
# #     #     for el in null_spacee:
# #     #         if sum(el) > 0:
# #     #             eigenvectors.append(el)

# #     # # ###print the eigenvectors
# #     # ###print("Eigenvectors:")
# #     # for eigenvector in eigenvectors:
# #     #     blu = (eigenvector + constants )%2
# #     #     ikl = subpcmm@blu%2
# #     #     if np.array_equal(np.array(ikl), np.array(P)@ np.array(insidesubsyndrome)%2):
# #     #         ###print("YAY")
# #     #         break














# #     # free_variables = [10,17,18]

# #     # U[30,:] = U[5,:]
# #     # y[30] = y[5]
# #     # U[22,:] = U[13,:]
# #     # y[22] = y[13]
# #     # U[19,:] = U[17,:]
# #     # y[19] = y[17]
# #     # U[24,:] = U[21,:]
# #     # y[24] = y[21]
# #     # U[28,:] = U[15,:]
# #     # y[28] = y[15]
# #     # U[32,:] = U[23,:]
# #     # y[32] = y[23]
# #     # U = np.delete(U, [5,11,13,15,17,20,21,23,31],axis=0)
# #     # y = np.delete(y, [5,11,13,15,17,20,21,23,31])
# #     # anzahl = U.shape[1]
# #     # U = list(U)
# #     # y = list(y)
# #     # deletethecolums = [5,11,13,15,17,20,21,23,31,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53]
# #     # for colu in deletethecolums:
# #     #     newline = np.zeros(anzahl,dtype=int)
# #     #     newline[colu] = 1
# #     #     U.append(newline)
# #     #     y.append(0)
# #     # U = np.array(U)        
# #     # y = np.array(y)
# #     # ###print(np.linalg.solve(U,y))
# #     import sys

# #     # Open a file in write mode to save the ###print output
# #     with open('subpcmoutputx.txt', 'w') as f:
# #         # Redirect the standard output to the file
# #         sys.stdout = f

# #         # Your ###print statements here
# #         for row in np.array(subpcmm):
# #             ###print(row)
        

# #     # Reset the standard output back to its original value (console)
# #     sys.stdout = sys.__stdout__
    

# #     # #delrows = 5, 11, 13, 15, 17, 20,21, 23
    
# #     # # U = np.delete(U, free_variables, axis= 0)
# #     # # y = np.delete(y, free_variables)
    
# #     # ###print("free variables",free_variables)


# # # [0,2,7,32,33,50]

















 
# # trial and error:
#     U = np.array(U).astype(int)
#     austauschbar, delcols = findaustauschbar(U)
#     ###print("delcols",delcols)
#     ###print(U)
#     U = GF(U)
#     ###print(austauschbar,delcols)
#     combinations = list(itertools.product(*austauschbar))#austauschbar[0], austauschbar[1], austauschbar[2],austauschbar[3],austauschbar[4],austauschbar[5],austauschbar[6]))#,austauschbar[7],austauschbar[8]))
#     merged = []
#     c =0
#     length = len(combinations)

#     ####print("hier4")

#     ####print(austauschbar, delcols)
# # merge alle austauschbar Listen zsm
#     for sublist in austauschbar:
#         merged.extend(sublist)

#     lowestweight = U.shape[1]
#     subsolution =0
#     length = len(combinations)
    
#     # ###print("combinations",combinations)
#     Ushape1 = U.shape[1]

#     subsolution =0
#     length = len(combinations)
#     ###print(length)
#     c=0
#     GF = galois.GF(2)
#     rang = 0
#     ###print("y",y)
#     y = list(y)

#     U = np.array(U)



#     for i,col in enumerate(U.T):
#         if rank(np.delete(U,col).T) == rank(U.T):
#             loesch.append(i)


#     for i,delete in enumerate(combinations):
#         combinations[i] = set(delete)
#     if delcols:
#         merged += delcols
#     merged = set(merged)
#     for delete in combinations:
#         mergednow = merged.copy()
#         # if delcols:
#         #     mergednow += delcols
#         # mergednow = list(mergednow)
#         # mergednow = [x for x in mergednow if x not in delete]
#         mergednow = mergednow - delete


#         # mergednow = sorted(mergednow)  
#         ynow = y.copy()
#         Unow = U.copy()
        
#         for index in mergednow:
#             addline = np.zeros(Ushape1)
#             addline[index] = 1
#             Unow = np.vstack((Unow,addline)).astype(int)
#             ynow.append(0)
#         # Unow = np.array(Unow).astype(int)
#         # ynow = np.array(ynow).astype(int)
#         # Unow = GF(Unow)
#         # ynow = GF(ynow)
#             ###print(rank(U.T),rank(Unow.T))

            


#         if Unow.shape[0] == Unow.shape[1]:
#             # solut = np.linalg.solve(Unow,ynow)
#             solut = solve(Unow,ynow)
            
#             solut = np.abs(np.array(solut).astype(int)) %2
#             predsyn = np.array(GF(subpcmm) @ GF(solut))
#             if np.array_equal(predsyn,np.array(insidesubsyndrome)):
#                 subsolution = solut
#                 # ###print("ynow",np.array(ynow))
#                 ###print("solut",solut)
#                 # ###print(insidesubsyndrome)
#                 # ###print(suberror)
#                 # ###print("delcols",mergednow)
#                 ###print('cluster successfully solved')
#                 # weight = sum(np.array(solut))
#                 # if weight < lowestweight:
#                 #     lowestweight = weight
#                 #     subsolution = solut
                
#                 # solut = list(np.array(solut).astype(int))
#             else:
#                 c += 1
#                 #####print("ERROR")
#                 if c%10 == 0:
#                     # ###print(c/length)
#                     pass
#     if c >= len(combinations)-1:
        
#         ####print("decode FAIL")
#         fullsolution = np.zeros(facegraph.number_of_nodes()*2).astype(int)
#     else:
#         fullsolution = np.zeros(facegraph.number_of_nodes()*2).astype(int)
#     for node in range(len(subsolution)):
#         if subsolution[node] == 1:
#             fullsolution[inside[node]-facegraph.number_of_nodes()] = 1






# #     errnodes = []
# #     for idx,i in enumerate(fullsolution):
# #         if i == 1:
# #             errnodes.append(idx)
    
# #     ###print(subsolution)
# #     ###print(inside)
# #     ###print(fullsolution)
#     return fullsolution



# # decodepred = gfinversion([[1,0,0,0,0,0,1],
# #             [1,1,1,0,0,0,0],
# #             [0,0,0,1,0,0,0],
# #             [0,0,1,1,1,0,0],
# #             [1,0,0,0,1,1,1]], [1,0,1,1,1])
# # ##########print(decodepred)

# cluster = np.zeros((6,22)).astype(int)
# cluster[0][0] =1
# cluster[0][1] =1
# cluster[0][2] =1
# cluster[0][6] =1
# cluster[0][7] =1
# cluster[0][8] =1 

# cluster[1][2] =1
# cluster[1][3] =1
# cluster[1][4] =1
# cluster[1][8] =1
# cluster[1][9] =1
# cluster[1][10] =1 

# cluster[2][5] =1
# cluster[2][6] =1
# cluster[2][7] =1
# cluster[2][12] =1
# cluster[2][13] =1
# cluster[2][14] =1 

# cluster[3][7] =1
# cluster[3][8] =1
# cluster[3][9] =1
# cluster[3][14] =1
# cluster[3][15] =1
# cluster[3][16] =1 

# cluster[4][11] =1
# cluster[4][12] =1
# cluster[4][13] =1
# cluster[4][17] =1
# cluster[4][18] =1
# cluster[4][19] =1

# cluster[5][13] =1
# cluster[5][14] =1
# cluster[5][15] =1
# cluster[5][19] =1
# cluster[5][20] =1
# cluster[5][21] =1 

# syn = np.array([1,1,1,1,1,1])
# # ##########print(gfinversion(cluster,syn))

