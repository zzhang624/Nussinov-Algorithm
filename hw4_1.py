import numpy as np

min_loop_length = 4



def pair_check(tup):
    if tup in [('A', 'U'), ('U', 'A'), ('C', 'G'), ('G', 'C')]:
        return True
    return False


def cal_OPT(sequence):
    """\
    Description:
    ------------
        Implement Nussinov algorithm(dynamic programming) and construct the OPT matrix that stores the optimal score. You may use bottom up approach or top down approach.
    Parameter:
    ------------
        sequence: RNA sequence of length n, type `list`.
    Return:
    ------------
        OPT: OPT chart, an np.array of the shape (n,n)
    """
    length = len(sequence)
    OPT = np.zeros((length,length), dtype=int) 
    for j in range(5, length): #OPT[i][j]
        for i in range (0, length - j ): # OPT[i][j+i], j is [j+i] here.
            maxs = [] #store max values
            maxs.append(OPT[i][j+i-1]) #none-match case
            
            for t in range(i,j+i-4): # t runs from i to j-4, [j+i] pairs to [t].
                tup = (sequence[t],sequence[j+i])
                if pair_check(tup): # check match
                    if t == i:
                        opt_value = 1 + OPT[t+1][j+i-1]
                    else:    
                        opt_value = 1 + OPT[i][t-1] + OPT[t+1][j+i-1]
                    maxs.append(opt_value)
                    
            OPT[i][j+i] = max(maxs)
    return OPT


def traceback(OPT, sequence, i, j, structure):
    """\
    Description:
    ------------
        Backtracking algorithm, to find the pairing of bases in RNA sequence according to DP chart.    
    Parameter:
    ------------
        OPT: OPT matrix from cal_OPT(sequence)
        sequence: RNA sequence of length n, type `list`.
        fill feel to include more parameters into the function aside from `OPT` and `sequence`...
    Return:
    ------------
        structure: 
            a list of tuples that stores the pairing of the optimal solution, 
            e.g. 
            For a sequence with connection of bases 1 and 2, 3 and 4(1, 2, 3, 4 are indices, counting from 0 to n-1 for sequence of length n),
            the returned structure should be [(1, 2), (3, 4)]
    """
    # Implement your algorithm here
    if OPT[i][j] == 0:
        return
    elif OPT[i][j] == OPT[i][j-1]:
        traceback(OPT, sequence, i, j-1, structure)
    elif OPT[i][j] == OPT[i+1][j-1] +1 and pair_check((sequence[i],sequence[j])):
        structure.append((i, j))
        traceback(OPT, sequence, i+1, j-1, structure)
    else:
        for t in range(i+1,j-4):
            if OPT[i][j] == 1 + OPT[i][t-1] + OPT[t+1][j-1] and pair_check((sequence[t],sequence[j])):
                structure.append((t, j))
                traceback(OPT, sequence, i, t-1, structure)
                traceback(OPT, sequence, t+1, j-1, structure)
                break


sequences = []
with open("test_data", "r") as fp:
    for seq in fp:
        sequences.append(seq.strip("\n"))

for i, seq in enumerate(sequences):
    OPT = cal_OPT(seq)
    
    # fill free to include more parameters according to the function that you defined.
    structure = []
    traceback(OPT, seq, 0, len(seq)-1, structure)
    print(structure)
    with open("result_" + str(i) + ".txt", "w") as fp:
        for pair in structure:
            fp.write(str(pair) + "\n")
    
