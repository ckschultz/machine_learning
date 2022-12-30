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

# Gloabl variables/hardcoded alpha/delta values
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

# Find the string sequence
def find_sequence(dpMatrix, str1, str2):
    seq1 = ''
    seq2 = ''
    row = len(str1)
    col = len(str2)

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
    return (seq1[::-1], seq2[::-1], dpMatrix[-1][-1])

# Constructing the dp matrix
def match(str1,str2):
    alpha = {'AA':0,'AC':110,'AG':48,'AT':94,'CA':110,'CC':0,'CG':118,'CT':48,'GA':48,'GC':118,'GG':0,'GT':110,'TA':94,'TC':48,'TG':110,'TT':0}
    delta = 30

    #dp[i][j]
    dp = [[0]*(len(str2)+1) for _ in range(len(str1)+1)]

    #Base Case
    for i in range(1,len(str1)+1):
        dp[i][0] = delta*i
    for j in range(1,len(str2)+1):
        dp[0][j] = delta*j

    for i in range(1,len(str1)+1):
        for j in range(1,len(str2)+1):
            dp[i][j] = min((alpha[str1[i-1]+str2[j-1]]+dp[i-1][j-1]),
                            delta + dp[i-1][j],
                            delta + dp[i][j-1]
                            )
    x,y,z = find_sequence(dp,str1,str2)

    return x,y,z,process_memory()

# Write main function to be called from the command line
if __name__ == '__main__':
    x, y = gen(sys.argv[1])

    # CALL DP FUNCTION 
    soln = time_wrapper(match,x,y)
    
    # PRINT STATS
    lines = [str(soln[2]), "\n"+ str(soln[0]), "\n"+ str(soln[1]), "\n"+ str(soln[3]), "\n"+ str(soln[4])]
    with open(sys.argv[2], 'w') as f:
        f.writelines(lines)
    

