seq='SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDTVYCPRHVICTAEDMLNPNYEDLLIRKSNHSFLVQAGNVQLRVIGHSMQNCLLRLKVDTSNPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNHTIKGSFLNGSCGSVGF'

lengthOfSequence=len(seq)


def printTheSequence(sequence):
    print(sequence)
def printLines(lengthOfSequence):
    print('|' * lengthOfSequence)
def printFinalString(finalString):
    print(finalString)
#AP contains parameters for alpha prediction
alphaPrediction = {
    'E': 1.53,
    'A': 1.45,
    'L': 1.34,
    'H': 1.24,
    'M': 1.2,
    'Q': 1.17,
    'W': 1.14,
    'V': 1.14,
    'F': 1.12,
    'K': 1.07,
    'I': 1.00,
    'D': 0.98,
    'T': 0.82,
    'R': 0.79,
    'S': 0.79,
    'C': 0.77,
    'N': 0.73,
    'Y': 0.61,
    'P': 0.59,
    'G': 0.53,
}
#otherAlphaPrediction contains reduced parameters for alpha prediction
otherAlphaPrediction = {
    'E': 1,
    'A': 1,
    'L': 1,
    'H': 1,
    'M': 1,
    'Q': 1,
    'W': 1,
    'V': 1,
    'F': 1,
    'K': 0.5,
    'I': 0.5,
    'D': 0,
    'T': 0,
    'R': 0,
    'S': 0,
    'C': 0,
    'N': -1,
    'Y': -1,
    'P': -1,
    'G': -1,
}

betaPrediction= {
    'M': 1.67,
    'V': 1.65,
    'I': 1.6,
    'C': 1.3,
    'Y': 1.29,
    'F': 1.28,
    'Q': 1.23,
    'L': 1.22,
    'T': 1.2,
    'W': 1.19,
    'A': 0.97,
    'R': 0.9,
    'G': 0.81,
    'D': 0.8,
    'K': 0.74,
    'S': 0.72,
    'H': 0.71,
    'N': 0.65,
    'P': 0.62,
    'E': 0.26,
}

otherBetaPrediction= {
    'M': 1,
    'V': 1,
    'I': 1,
    'C': 1,
    'Y': 1,
    'F': 1,
    'Q': 1,
    'L': 1,
    'T': 1,
    'W': 1,
    'A': 0.5,
    'R': 0,
    'G': 0,
    'D': 0,
    'K': -1,
    'S': -1,
    'H': -1,
    'N': -1,
    'P': -1,
    'E': -1,
}

#initialising list which will store at indexes if possible secondary structure or not
alp=[]
for i in range(lengthOfSequence) :
    alp.append(0)
bet = []
for i in range(lengthOfSequence):
    bet.append(0)



def nextAlpha(seq, ind):   #will return the next index from where the alpha series will start
  traversalLength   = lengthOfSequence - 5
  for i in range(ind,traversalLength):
    currentSequene=seq[(i+1) -1:(i+6)  - 1]
    # print(currentSequene)

    hf = 0
    hb = 0
    score = 0

    #number of helix former and breaker, and initialise score
    for j in currentSequene:
      if otherAlphaPrediction[j]  <= 0:
        hb+=1
      else:
        hf+=1
      score+= otherAlphaPrediction[j]
    if(score>=4 and hb<2 and hf>3):   #these are the condition which need to be satisfied if have to be alpha
      return i
def scoringCondition(score):
    return score < 4 

def break3(val , ind , i):
    return val < 4 and i ==  ind-4
def checkBreakCondition(value , index , i):
    return value < 4 and i == index + 6
def compareAplhaSelection(seq, ind, alp):    #this will take index and see before and after the index and will return the index uptil which the alpha will propogate
  id1, id2= ind,ind+6
  for i in range(ind-4, -1, -1):
    currentSequene=seq[((i+1) -1):((i+5) -1 )]
    score=0
    for j in currentSequene:
      score+= alphaPrediction[j]
    # print("inside 1st, score=",score)
    if(break3(score ,  ind , i)):
      break
    if(scoringCondition(score)):
      id1=i+1
      break
    if(score>4 and i==0):
      id1=0
    
    for i in range(id2, lengthOfSequence-3):
        tempItr  =i +1 
        currentSequene=seq[tempItr -1 :tempItr +3]
        score=0
        for j in currentSequene:
            score+= alphaPrediction[j]
            if(checkBreakCondition(score , ind  , i)):
                break
            if(score<4 or i==lengthOfSequence-3):
                id2=i
                break
  # print("ind1, ind 2 = ", id1, id2)     
    for i in range(id1, id2):     #here marking 1 on list the part where it can be alpha
        alp[i]=1
  return id2, alp

#this function will use above 2 function and send us the list with possiblity of alpha is marked as 1 in that location
def possibleAlpha(seq):
  alp=[]
  for i in range(lengthOfSequence):
      alp.append(0)

  i=0
  while i<lengthOfSequence:
    ind=nextAlpha(seq,i)
    i=ind
    if(ind==None):
      break
    for j in range(6):
      alp[(i+j)]=1
    ind2, alp=compareAplhaSelection(seq, ind, alp)
    i=ind2-1
    i+=1
  return alp


bet=[]
for i in range(lengthOfSequence):
    bet.append(0)

def incrementChecker(val , hb, hf):
    return val > 3 and hb < 2 and hf > 2
def nextBeta(seq, ind):   #will return the next index from where the alpha series will start
  for i in range(ind,lengthOfSequence-4):
    itrRange = i +1
    currentSequene=seq[itrRange -1 :itrRange +4]
    # print(currentSequene)
    hf = 0
    hb = 0
    score = 0 
    #number of helix former and breaker, and initialise score
    for j in currentSequene:
      if otherBetaPrediction[j]  > 0:
        hf+=1
      else:
        hb+=1
      score+= otherAlphaPrediction[j]
    if(incrementChecker(score, hb ,hf)):
      return i


def anotherBreakCondition(val , ind  , i):
    return val < 4 and i == ind -4

def break2(val  , ind , i):
    return val < 4 and i==ind + 5

def compareBetaSelection(seq, ind, bet):    #this will take index and see before and after 4 value returns the index uptil which the alpha will propogate
  id1 = ind
  id2 = ind+5

  for i in range(ind-4, -1, -1):
    currentSequene=seq[i:i+4]
    score=0
    for j in currentSequene:
      score+= betaPrediction[j]
    if(anotherBreakCondition(score , ind , i)):
      break
    if(scoringCondition(score)):
      id1=i+1
      break
    if(score> 4 and i==0):
      id1=0

  for i in range(id2, lengthOfSequence-3):
    currentSequene=seq[i:i+4]
    score=0
    for j in currentSequene:
      score+= betaPrediction[j]
    if(break2(score,  ind , i)):
      break
    if(otherBreakCondition(score , i , lengthOfSequence) ):
      id2=i
      break
  for i in range(id1, id2):
    bet[i]=1
  return id2, bet

def otherBreakCondition(val , i , ind):
  return val < 4 or i == ind -3

def findNextBeta(seq):
  bet=[0]*lengthOfSequence
  i=0
  while i<lengthOfSequence:
    ind=nextBeta(seq,i)
    i=ind
    if(ind==None):
      break
    for j in range(5):
      bet[(i+j)]=1
    ind2, bet=compareBetaSelection(seq, ind, bet)
    i=ind2-1
    i+=1
  return bet

alpha= possibleAlpha(seq)
beta=findNextBeta(seq)

print(alpha)
print(beta)

#This is not the final answer, here we are not considering the conflicting case(this has been handled later)
finallist=['T']*len(seq)

def returnAddValue(i):
    if(alpha[i] == 1):
        return 'H'
def returnOtherAddValue(i):
    if(beta[i] == 'H'):
        return 'S'
for i  in range(lengthOfSequence):
  if(alpha[i]==1):
    finallist[i]=returnAddValue(i)
  if(beta[i]==1):
    finallist[i]=returnOtherAddValue(i)

# print(finallist)

#Here we are finding the indexes which has been predicted both alpha and beta by our method
conflict=[]
for i in range(lengthOfSequence):
  if(alpha[i]==1 and beta[i]==1):
    conflict.append(i)

print(conflict)

def returnH():
  return 'H'
def returnS():
  return 'S'

conflictvalue= [[0]*len(conflict)]
indices=[]
for i in range((len(conflict)+1)-2):

  if(conflict[i]+1== conflict[i+1] and i!= len(conflict)-2):
    indices.append(conflict[i])
  else:
    indices.append(conflict[i])
    if(i == len(conflict)-2):
      indices.append(conflict[i+1])
    # print(indices)
    salpha, sbeta=0,0
    for j in indices:
      salpha+= alphaPrediction[seq[j]]
      sbeta+= betaPrediction[seq[j]]
    if(salpha>=sbeta):
      for k in indices:
        finallist[k]=returnH()
    else:
      for k in indices:
        finallist[k]=returnS()
    indices=[]

print("The finallist= ", finallist)

finalstring=''
for i in finallist:
  finalstring+= i


printTheSequence(seq)

printLines(lengthOfSequence)

printFinalString(finalstring)

found1 = "TTHHHHHHTTTTTTTSSSSSSSSSSTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTHSSSSSSSSSSSSSHTTTTTTTHHHHHHTTTTTTTTTTHSSSSSTTTTSSSSSSTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
found2 = "TTTT     HHHHHH EEEEEETTEEEEEEEETTEEEEEGGGG  HHHHH   HHHHHHH  GGG EEEETTEEE EEEEEEETTEEEEEE   TTTT        TTTEEEEEEEEETTEEEEEEEEEETTTT B    TTTTTTTEE "
for i, j in enumerate(found1):
  print(f"At index {i+1} -  {j} and {found2[i]}")