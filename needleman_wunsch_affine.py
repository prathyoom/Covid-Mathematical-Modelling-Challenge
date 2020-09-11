from matplotlib import pyplot as plt

def count_freq(arr):
    freq = {}
    for items in arr:
        freq[items] = arr.count(items)
    a1 = []
    a2 = []
    for key, value in freq.items():
        a1.append(key)
        a2.append(value)
    return a1, a2


BLOSUM62 = {}
amino = "ARNDCQEGHILKMFPSTWYV"
matrix =[[ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0],

 [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3],

 [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3],

 [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3],

 [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],

 [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2],

 [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2],

 [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3],

 [-2, 0, 1,-1,-3,-0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3],

 [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3],

 [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1],

 [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2],

 [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1],

 [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1],

 [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2],

 [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2],

 [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0],

 [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3],

 [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1],

 [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4]
]

#BLOSUM62 initialization
def create_blosum():
    for i in range(len(amino)):
        for j in range(len(amino)):
            temp = (amino[i],amino[j])
            BLOSUM62[temp] = matrix[i][j]
    return

def needleman_wunsch_affine(s1, s2, gap_open, gap_extend):
    n = len(s1)
    m = len(s2)
    d = [[0 for i in range(m+1)] for j in range(n+1)]       #Table d indicates the matches/mismatches
    a = [[0 for i in range(m+1)] for j in range(n+1)]       #Table a indicates the gaps extended in s1
    b = [[0 for i in range(m+1)] for j in range(n+1)]       #Table b indicates the gaps extended in s2
    t = [[0 for i in range(m+1)] for j in range(n+1)]       #Table t is utilized for tracing back the resulting alignment
    string_show = False

    #Initializing d
    d[0][0] = 0
    d[0][1] = gap_open + gap_extend
    d[1][0] = gap_open + gap_extend
    for i in range(2,len(d)):
        d[i][0] = d[i-1][0] + gap_extend
    for i in range(2,len(d[0])):
        d[0][i] = d[0][i-1] + gap_extend

    #Inititalizing a
    a[0][1] = gap_open + gap_extend
    for i in range(len(a)):
        a[i][0] = -10000000000
    for i in range(2,len(a[0])):
        a[0][i] = a[0][i-1] + gap_extend
    
    #Initializing b
    b[1][0] = gap_open + gap_extend
    for i in range(len(b[0])):
        b[0][i] = -10000000000
    for i in range(2,len(b)):
        b[i][0] = b[i-1][0] + gap_extend

    #Initializing t
    for i in range(len(t)):
        t[i][0] = 'b'
    for i in range(len(t[0])):
        t[0][i] = 'a'
    t[0][0] = '-'

    #Updating all the states/matrices
    for i in range(1,len(d)):
        for j in range(1,len(d[0])):
            b[i][j] = max(b[i-1][j] + gap_extend, d[i-1][j] + gap_extend + gap_open)
            a[i][j] = max(a[i][j-1] + gap_extend, d[i][j-1] + gap_extend + gap_open)
            scr = BLOSUM62[(s1[i-1],s2[j-1])]
            d[i][j] = max(d[i-1][j-1] + scr, a[i][j], b[i][j])
            if d[i][j] == d[i-1][j-1] + scr:
                t[i][j] = 'd'
            elif d[i][j] == a[i][j]:
                t[i][j] = 'a'
            elif d[i][j] == b[i][j]:
                t[i][j] = 'b'
        
    res1 = ""
    res2 = ""
    i = n
    j = m

    #Trace back to obtain resultant strings
    while(i > 0 or j > 0):
        if t[i][j] == 'd':
            res1 += s1[i-1]
            res2 += s2[j-1]
            i -= 1
            j -= 1
        elif t[i][j] == 'b':
            res1 += s1[i-1]
            res2 += '_'
            i -= 1
        elif t[i][j] == 'a':
            res1 += '_'
            res2 += s2[j-1]
            j -= 1

    #Reversing the resultant strings
    res1 = res1[::-1]
    res2 = res2[::-1]

    if(string_show):
        print("\nResultant strings: \n")
        print(res1)
        print(res2)

    #for i in range(len(d)):
    #   print(*d[i])
    #for i in range(len(a)):
    #   print(*a[i])
    #for i in range(len(b)):
    #   print(*b[i])
    return res1, res2, d[n][m]

def show_identity(s1,s2):
    iden = 0
    for i in range(len(s1)):
        if s1[i] == s2[i]:
            iden += 1
    #print("Identity: " + str(round(float(iden/len(s1)),3)))
    return float(iden/len(s1))

def show_similarity(s1,s2):
    simi = 0
    for i in range(len(s1)):
        if s1[i] == "_" or s2[i] == "_":
            continue
        if BLOSUM62[(s1[i],s2[i])] > 0:
            simi += 1
    #print("Similarity: " + str(round(float(simi/len(s1)),3)))
    return float(simi/len(s1))

create_blosum()
s1 = input()
s2 = input()

res1, res2, score = needleman_wunsch_affine(s1,s2,-10.1,-10.1)
#print(res1)
#print(res2)
#print(0.505*show_identity(res1,res2) + 0.005 * show_similarity(res1,res2) + 0.49 * score)
#exit()

#res1, res2, score = needleman_wunsch_affine(s1,s2,-20,-20)
iden = [[0 for _ in range(2000)] for lol in range(2000)]
sim = [[0 for _ in range(2000)] for lol in range(2000)]
scores = [[0 for _ in range(2000)] for lol in range(2000)]
max_score = -1000000000
min_score = 1000000000
max_sim = -1000000000
min_sim = 1000000000
max_iden = -1000000000
min_iden = 1000000000
counter1 = 0
counter2 = 0

for i1 in range(-20,0,1): 
    counter2 = 0
    i = (float)(i1/10)
    for j1 in range(i1,0,1):
        j = (float)(j1/10)
        res1, res2, score = needleman_wunsch_affine(s1,s2,i,j)
        sim[counter1][counter2] = show_similarity(res1,res2)
        iden[counter1][counter2] = show_identity(res1,res2)
        max_score = max(max_score,score)
        min_score = min(min_score,score)
        max_sim = max(max_sim,sim[counter1][counter2])
        min_sim = min(min_sim,sim[counter1][counter2])
        max_iden = max(max_iden,iden[counter1][counter2])
        min_iden = min(min_iden,iden[counter1][counter2])
        scores[counter1][counter2] = score
        counter2 += 1
    counter1 += 1

counter1 = 0
counter2 = 0
for i1 in range(-20,0,1):
    counter2 = 0
    i = (float)(i1/10)
    for j1 in range(-20,0,1):
        j = (float)(j1/10)
        if max_score == min_score:
            scores[counter1][counter2] = max_score
        else:
            scores[counter1][counter2] = (float)((scores[counter1][counter2] - min_score)/(max_score-min_score))
        if max_sim == min_sim:
            sim[counter1][counter2] = max_sim
        else:
            sim[counter1][counter2] = (float)((sim[counter1][counter2] - min_sim)/(max_sim - min_sim))
        if max_iden == min_iden:
            iden[counter1][counter2] = max_iden
        else:
            iden[counter1][counter2] = (float)((iden[counter1][counter2] - min_iden)/(max_iden - min_iden))
        counter2 += 1     
    counter1 +=  1

new_scores = [[0 for _ in range(202)] for lol in range(202)]
counter1 = 0
counter2 = 0

suma = 0
upper_bound = 1

for gamma1 in range(0,10,1):
    counter2 = 0
    gamma = (float)(gamma1/10)
    for beta1 in range(0,10-gamma1,1):
        beta = (float)(beta1/10)
        alpha = 1 - gamma - beta
        for i in range(0,20,1) : 
            for j in range(0,20,1):    
                new_scores[counter1][counter2] += (float)(alpha * iden[i][j] + beta * sim[i][j] + gamma * scores[i][j])           
                #print(alpha,beta,gamma,iden[i][j],sim[i][j],scores[i][j],new_scores[counter1][counter2])
        new_scores[counter1][counter2]/=(400)
        counter2 += 1
    counter1 += 1

#print(new_scores)
best = -10000000

for i in range(0,20):
    for j in range(0,i):
        if best < new_scores[i][j]:
            best = new_scores[i][j]
            final_gamma = (float)(i/20)
            final_beta = (float)(j/20)
            final_alpha = (1 - final_gamma - final_beta)

for i in range(20):
    for j in range(i):
        if best == new_scores[i][j]:
            final_gamma = (float)(i/20)
            final_beta = (float)(j/20)
            final_alpha = (1 - final_gamma - final_beta)
            #print(final_alpha,final_beta,final_gamma)

maxi=-100000000
best_gap_open = -20.0
best_gap_ext = -20.0

for i1 in range(-20,0,1): 
    counter2 = 0
    i = (float)(i1/10)
    for j1 in range(i1,0,1):
        j = (float)(j1/10)
        val = (final_alpha * iden[counter1][counter2] + final_beta * sim[counter1][counter2] + final_gamma  * scores[counter1][counter2])
        if val > maxi:
            maxi = val
            best_gap_open = i
            best_gap_ext = j
        counter2+=1
    counter1+=1

print(final_alpha,final_beta,final_gamma)

print(best_gap_open,best_gap_ext)

res1,res2,score = needleman_wunsch_affine(s1,s2,best_gap_open,best_gap_open)
print(res1)
print(res2)
print(score)
