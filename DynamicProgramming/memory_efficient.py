import sys
import time
import psutil

def process_memory():
    process = psutil.Process()
    memory_info = process.memory_info()
    memory_consumed = int(memory_info.rss/1024)
    return memory_consumed

def time_wrapper(func,x,y):
    start_time = float(time.time())
    x,y,z,mem = func(x,y)
    end_time = float(time.time())
    end_time = (end_time - start_time)*1000
    return x,y,z,end_time,mem

alpha = {'AA':0,'AC':110,'AG':48,'AT':94,'CA':110,'CC':0,'CG':118,'CT':48,'GA':48,'GC':118,'GG':0,'GT':110,'TA':94,'TC':48,'TG':110,'TT':0}
delta = 30
def gen(txt):
    str1 = None
    str2 = None
    ind1 = []
    ind2 = []

    #read file and get values
    with open(txt, 'r') as file:
        current = None
        for line in file:
            line = line.strip()
            try:
                current.append(int(line))
            except:
                if current == None:
                    str1 = line
                    current = ind1
                else:
                    str2 = line
                    current = ind2
        file.close()
    len1 = len(str1)
    len2 = len(str2)
    #generate strings
    for ind in ind1:
        str1 = str1[:ind+1] + str1 + str1[ind+1:]

    for ind in ind2:
        str2 = str2[:ind+1] + str2 + str2[ind+1:]

    if len(str1) != 2**len(ind1)*len1:
        return 'str1 error'
    if len(str2) != 2**len(ind2)*len2:
        return 'str2 error'
    return (str1,str2)

def find_sequence(dpMatrix, str1, str2):
    seq1 = ''
    seq2 = ''
    row = len(str1)
    col = len(str2)
    # Creating the sequence from the dp matrix
    while (row > 0 or col > 0):
        tempMin = min(delta + dpMatrix[row][col-1], delta + dpMatrix[row-1][col], alpha[str1[row-1]+str2[col-1]] + dpMatrix[row-1][col-1])
        if tempMin == alpha[str1[row-1]+str2[col-1]] + dpMatrix[row-1][col-1]:
            seq1 += str1[row-1]
            seq2 +=str2[col-1]
            row -= 1
            col -= 1
        elif tempMin == delta + dpMatrix[row][col-1]:
            seq1 += '_'
            seq2 +=str2[col-1]
            col -= 1
        else:
            seq1 += str1[row-1]
            seq2 +='_'
            row -= 1

    # return (Sequence 1, Sequence 2, Cost of Alignment)
    return (seq1[::-1], seq2[::-1], dpMatrix[-1][-1],None)


def match(x,y):
    #Def pentalties
    #Base Cases for Recusrion. Stop when X has len < 2
    if len(x) <=2 and len(y) != 0:

        dp = [[0]*(len(y)+1) for _ in range(len(x)+1)]

        #Base Case
        for i in range(1,len(x)+1):
            dp[i][0] = delta*i
        for j in range(1,len(y)+1):
            dp[0][j] = delta*j

        # Find optimal solution for x of len <= 2
        for i in range(1,len(x)+1):
            for j in range(1,len(y)+1):
                dp[i][j] = min((alpha[x[i-1]+y[j-1]]+dp[i-1][j-1]),
                                delta + dp[i-1][j],
                                delta + dp[i][j-1]
                                )
        return find_sequence(dp, x, y) # Calling the function to find sequence
    elif len(x) <=2 and len(y) == 0:
        return(x,'_'*len(x), delta*len(x),None)


    #Split X into 
    half1 = x[:int(len(x)/2)]
    half2 = x[int(len(x)/2):]

    #Solve First Half
    x = half1

    #Inorder to save memory and not deal with the complexities of python garbage collection.
    #I will initiate a dp matrix of size half2 since this will always be bigger than or equal to half1
    dp = [[0]*(len(y)+1) for _ in range(2)]

    #Base Case for DP program
    for j in range(1,len(y)+1):
        dp[0][j] = delta*j


    #DP loop to fill matrix
    for i in range(1,len(x)+1):
        if i%2:
            dp[1][0] = delta + dp[0][0]
            for j in range(1,len(y)+1):
                    dp[1][j] = min((alpha[x[i-1]+y[j-1]]+dp[0][j-1]),
                                    delta + dp[0][j],
                                    delta + dp[1][j-1]
                                    )
        else:
            dp[0][0] = delta + dp[1][0]
            for j in range(1,len(y)+1):
                    dp[0][j] = min((alpha[x[i-1]+y[j-1]]+dp[1][j-1]),
                                    delta + dp[1][j],
                                    delta + dp[0][j-1]
                                    )

    #Save important row information give new pointer.
    if len(x)%2:
        opt1 = dp[1].copy()
    else:
        opt1 = dp[0].copy()
    
    x = half2[::-1]
    #Reassign the same values so that we are memory efficient and we do not have to wait on the garbage collector to clear memory.
    for i in range(2):
        for j in range(len(dp[i])):
            dp[i][j] = 0 

    #Base Case
    for j in range(1,len(y)+1):
        dp[0][j] = delta*j

    #DP loop to fill matrix
    for i in range(1,len(x)+1):
        if i%2:
            dp[1][0] = delta + dp[0][0]
            for j in range(1,len(y)+1):
                    dp[1][j] = min((alpha[x[i-1]+y[-j]]+dp[0][j-1]),
                                    delta + dp[0][j],
                                    delta + dp[1][j-1]
                                    )
        else:
            dp[0][0] = delta + dp[1][0]
            for j in range(1,len(y)+1):
                    dp[0][j] = min((alpha[x[i-1]+y[-j]]+dp[1][j-1]),
                                    delta + dp[1][j],
                                    delta + dp[0][j-1]
                                    )

    if len(x)%2:
        opt2 = dp[1].copy()
    else:
        opt2 = dp[0].copy()

    # Calculate alignment at split bestI
    bestMatch = opt1[0] + opt2[-1]
    bestI = 0
    for i in range(len(opt1)):
        if opt1[i]+opt2[-i-1] < bestMatch:
            bestMatch = opt1[i]+opt2[-i-1]
            bestI = i
    # Recursively call function again
    x1,y1,z1,zz = match(half1,y[:bestI])
    x2,y2,z2,zz = match(half2,y[bestI:])

    return ((x1+x2),(y1+y2),(z1+z2),process_memory())

# Main function enables calling from the command line
if __name__ == '__main__':
    x, y = gen(sys.argv[1])

    # CALL DP FUNCTION, time calculations
    soln = time_wrapper(match,x,y)
    
    # PRINT STATS
    lines = [str(soln[2]), "\n"+ str(soln[0]), "\n"+ str(soln[1]), "\n"+ str(soln[3]), "\n"+ str(soln[4])]
    with open(sys.argv[2], 'w') as f:
        f.writelines(lines)
