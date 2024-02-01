import random

# Class for (only 2D) binary Matrix with operator overloading for addition and multiply
class BinaryMatrix:
    def __init__(self, M):
        self.Nrow=len(M)
        self.Ncol=len(M[0])
        self.M=[[(M[j][i] %2) for i in range(self.Ncol)] for j in range(self.Nrow)]
        self.Mt=[[(M[j][i] %2) for j in range(self.Nrow)] for i in range(self.Ncol)]
           
    def __add__(self, x):
        if self.Nrow != x.Nrow or self.Ncol != x.Ncol:
           return False
        mat = [[ ((self.M[j][i] + x.M[j][i])% 2) for i in range(self.Ncol)] for j in range(self.Nrow)]
        return BinaryMatrix(mat)

    def __mul__(self, x):
        if self.Ncol != x.Nrow:
            return False
        mat = [[0 for i in range(x.Ncol)] for j in range(self.Nrow)]
        for i in range(self.Nrow):
            for j in range(x.Ncol):
                mat[i][j] = 0 
                for k in range(self.Ncol):
                    mat[i][j] += self.M[i][k] * (x.M[k][j])
                mat[i][j] = mat[i][j] % 2 
        return BinaryMatrix(mat)

def NibbleTo2DbitArray(nible):
    M=[]
    for i in range(4):
        bit=[(nible&0x1)]
        M.append(bit)
        nible=nible>>1
    return M

def AnalyseString(Istring, MS_LS, GM, HM, RM, ErrorFraction):
    G=BinaryMatrix(GM)
    H=BinaryMatrix(HM)
    R=BinaryMatrix(RM)

    yLSstring=""
    yMSstring=""
    Errstring=""
    errorSyndromeString=""
    DetectOkSting=""
    DataCorrectString=""

    TotalNbits=0
    NflippedBits=0
    Nnibbles=0
    NwrongDetects=0
    NwrongReceived=0
    for element in Istring:
        Nnibbles=Nnibbles+1
        asci=(ord(element)) >> MS_LS
       
#encode 
        xM=NibbleTo2DbitArray(asci)
        x=BinaryMatrix(xM)
        y=G*x
#print encoded data
        sum1=0
        for n in range(y.Nrow):
            sum1=sum1+y.M[n][0]*2**n
        yLSstring=yLSstring+f"{sum1:02x}"[1]
        yMSstring=yMSstring+f"{sum1:02x}"[0]   
        
#Introduce random Errors in data with parity
        ErrorVector=random.choices([1, 0], weights=[ErrorFraction, 1-ErrorFraction], k=len(y.M))
        NumberOfErrors=ErrorVector.count(1)
        NflippedBits=NflippedBits+NumberOfErrors
        TotalNbits=TotalNbits+len(y.M)
        Errstring=Errstring+str(NumberOfErrors)
        ErrorV=BinaryMatrix([ErrorVector])
        ErrorM=ErrorV.Mt
        Error=BinaryMatrix(ErrorM)
        yError =Error+y
#Error detect
        errorSyndrome=H*yError
        print (errorSyndrome.M)
        SyndromeValue=0
        for n in range(errorSyndrome.Nrow):
            SyndromeValue=SyndromeValue+errorSyndrome.M[n][0]*2**n
        print (SyndromeValue)
        errorSyndromeString=errorSyndromeString+str(SyndromeValue)
#errorcorrect for hamming [7,4]
        if SyndromeValue!=0 and H.Ncol == 7:
           toggleM=[[0] for j in range(H.Ncol)]
           toggleM[SyndromeValue-1][0] = 1
           toggle=BinaryMatrix(toggleM)
           yError = toggle+yError
           print (toggleM)
           print (ErrorM)
#errorcorrect for hamming [8,4]
        if SyndromeValue>7 and H.Ncol == 8:   #parity error (single bit errors)
           SyndromeValue = SyndromeValue-8
           print (ErrorM)
           toggleM=[[0] for j in range(H.Ncol)]
           toggleM[SyndromeValue-1][0] = 1
           toggle=BinaryMatrix(toggleM)
           yError = toggle+yError
           print (toggleM)
        elif SyndromeValue>1 and H.Ncol == 8:   #parity error (double (or more) bit errors)
           print (f"{ErrorM} double")
           toggleM=[[0] for j in range(H.Ncol)]
           toggleM[SyndromeValue-1][0] = 1
           toggle=BinaryMatrix(toggleM)
           yError = toggle+yError
           print (toggleM)
        
#Decode
        yDecode=R*yError
        print(x.M)
        print(yDecode.M)
#evaluate
        if ((yDecode.M!=xM and SyndromeValue==0) or (yDecode.M==xM and SyndromeValue!=0) ):
                DetectOkSting=DetectOkSting+'X'
                NwrongDetects=NwrongDetects+1
        else:
                DetectOkSting=DetectOkSting+'-'
        if (yDecode.M==xM):
                DataCorrectString=DataCorrectString+'-'
        else:
                DataCorrectString=DataCorrectString+'X'
                NwrongReceived=NwrongReceived+1

    print('yLS  ', yLSstring)
    print('yMS  ', yMSstring)
    print('#Err ', Errstring,' = ',f"{(100*NflippedBits/TotalNbits):.3g}",'%')
    print('ErDet', errorSyndromeString)
    print('DetOk', DetectOkSting,' = ',f"{(100*NwrongDetects/Nnibbles):.3g}",'%')
    print('DatOk', DataCorrectString,' = ',f"{(100*NwrongReceived/Nnibbles):.3g}",'%')
    return

GparM= [[1, 0, 0, 0],
[0, 1, 0, 0],
[0, 0, 1, 0],
[0, 0, 0, 1],
[1, 1, 1, 1]]

HparM = [[1, 1, 1, 1, 1]]

RparM= [[1, 0, 0, 0, 0],
[0, 1, 0, 0, 0],
[0, 0, 1, 0, 0],
[0, 0, 0, 1, 0]]

GHam74 = [[1, 0, 0, 0],
          [0, 1, 0, 0],
          [0, 0, 1, 0],
          [0, 0, 0, 1],
          [0, 1, 1, 1],
          [1, 0, 1, 1],
          [1, 1, 0, 1]]

HHam74  = [[1, 0, 1, 0, 1, 0, 1],
           [0, 1, 1, 0, 0, 1, 1],
           [0, 0, 0, 1, 1, 1, 1]] 

RHam74  =  [[1, 0, 0, 0, 0, 0, 0],
           [0, 1, 0, 0, 0, 0, 0],
           [0, 0, 1, 0, 0, 0, 0],
           [0, 0, 0, 1, 0, 0, 0]]

GHam84 = [[1, 1, 0, 1],
          [1, 0, 1, 1],
          [1, 0, 0, 0],
          [0, 1, 1, 1],
          [0, 1, 0, 0],
          [0, 0, 1, 0],
          [0, 0, 0, 1],
          [1, 1, 1, 0]]

HHam84  = [[1, 0, 1, 0, 1, 0, 1, 0],
           [0, 1, 1, 0, 0, 1, 1, 0],
           [0, 0, 0, 1, 1, 1, 1, 0],
           [1, 1, 1, 1, 1, 1, 1, 1]] 

RHam84  = [[0, 0, 1, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 1, 0, 0, 0],
           [0, 0, 0, 0, 0, 1, 0, 0],
           [0, 0, 0, 0, 0, 0, 1, 0]]

ErrorFraction = float(input('Give fraction of bit errors (number between 0 and 1)'))
Istring = input('Text :')
print('parity')
AnalyseString(Istring, 0, GparM, HparM, RparM, ErrorFraction)
AnalyseString(Istring,  4, GparM, HparM, RparM, ErrorFraction)
print('hamming7_4')
AnalyseString(Istring, 0, GHam74,HHam74 , RHam74, ErrorFraction)
AnalyseString(Istring,  4, GHam74, HHam74, RHam74, ErrorFraction)
print('hamming8_4')
AnalyseString(Istring, 0, GHam84, HHam84, RHam84, ErrorFraction)
AnalyseString(Istring,  4, GHam84, HHam84, RHam84, ErrorFraction)

