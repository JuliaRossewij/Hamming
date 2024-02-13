import random
from re import X
from timeit import default_timer

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

def ByteTo2DbitArray(byte):
    M=[]
    for i in range(8):
        bit=[(byte&0x1)]
        M.append(bit)
        byte=byte>>1
    return M


def bitArray2Value(y):
    sum1=0
    for n in range(y.Nrow):
        sum1=sum1+y.M[n][0]*2**n
    return sum1


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
    CorrectedDataString = ""

    TotalNbits=0
    NflippedBits=0
    Nnibbles=0
    NwrongDetects=0
    NwrongReceived=0
    start = default_timer()
    for element in Istring:
        Nnibbles=Nnibbles+1
        asci=(ord(element)) >> MS_LS
       
#encode 
        xM=NibbleTo2DbitArray(asci)
        x=BinaryMatrix(xM)
        y=G*x
#print encoded data
        sum1=bitArray2Value(y)
#        sum1=0
#        for n in range(y.Nrow):
#            sum1=sum1+y.M[n][0]*2**n
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
        SyndromeValue=bitArray2Value(errorSyndrome)
#        SyndromeValue=0
#        for n in range(errorSyndrome.Nrow):
#            SyndromeValue=SyndromeValue+errorSyndrome.M[n][0]*2**n
        errorSyndromeString=errorSyndromeString+hex(SyndromeValue)[2]
#errorcorrect for hamming [7,4]
        if SyndromeValue!=0 and H.Ncol == 7:
           toggleM=[[0] for j in range(H.Ncol)]
           toggleM[SyndromeValue-1][0] = 1
           toggle=BinaryMatrix(toggleM)
           yError = toggle+yError
#errorcorrect for hamming [8,4]
        if SyndromeValue>7 and H.Ncol == 8:   #parity error (single bit errors)
           SyndromeValue = SyndromeValue-8
           toggleM=[[0] for j in range(H.Ncol)]
           toggleM[SyndromeValue-1][0] = 1
           toggle=BinaryMatrix(toggleM)
           yError = toggle+yError
        elif SyndromeValue>1 and H.Ncol == 8:   #parity error (double (or more) bit errors)
           toggleM=[[0] for j in range(H.Ncol)]
           toggleM[SyndromeValue-1][0] = 1
           toggle=BinaryMatrix(toggleM)
           yError = toggle+yError                  
#Decode
        yDecode=R*yError
        CorrectedDataString=CorrectedDataString+hex(bitArray2Value(yDecode))[2]
#        print(x.M)
#        print(yDecode.M)
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
    duration = default_timer() - start
    print('yLS  ', yLSstring)
    print('yMS  ', yMSstring)
    print('#Err ', Errstring,' = ',f"{(100*NflippedBits/TotalNbits):.3g}",'%')
    print('ErSyn', errorSyndromeString)
    print('DetOk', DetectOkSting,' = ',f"{(100*NwrongDetects/Nnibbles):.3g}",'%')
    print('ycorr', CorrectedDataString)
    print('DatOk', DataCorrectString,' = ',f"{(100-100*NwrongReceived/Nnibbles):.3g}",'%')
    print('Time=', duration)
    return

#This function was added only to measure the time gain with lookup tables (LUTs) for hamming84. Is not tested for other cdong schemes
def AnalyseStringLUTHamming84(Istring, MS_LS, GM, HM, RM, ErrorFraction):
    G=BinaryMatrix(GM)
    H=BinaryMatrix(HM)
    R=BinaryMatrix(RM)

    Glut = []
    for nible in range(0, 16):
        xM=NibbleTo2DbitArray(nible)
        x=BinaryMatrix(xM)
        y=G*x
        Glut.append(y.M)

    Hlut = []
    Rlut = []
    for byte in range(0, 256):
       yM=ByteTo2DbitArray(byte)
       y=BinaryMatrix(yM)
       syndrome=H*y
       Hlut.append(bitArray2Value(syndrome))
       yreceived=R*y
       Rlut.append(yreceived.M)

    yLSstring=""
    yMSstring=""
    Errstring=""
    errorSyndromeString=""
    DetectOkSting=""
    DataCorrectString=""
    CorrectedDataString = ""

    TotalNbits=0
    NflippedBits=0
    Nnibbles=0
    NwrongDetects=0
    NwrongReceived=0
    start = default_timer()
    for element in Istring:
        Nnibbles=Nnibbles+1
        asci=(ord(element)) >> MS_LS
#encode 
        xM=NibbleTo2DbitArray(asci)
 #       x=BinaryMatrix(xM)
 #       y=G*x
 #       print(y.M)
 #       print(Glut[asci&0xF])
        y=BinaryMatrix(Glut[asci&0xF])
#print encoded data
        sum1=bitArray2Value(y)
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
#        errorSyndrome=H*yError
#        SyndromeValue=bitArray2Value(errorSyndrome)
#        print(SyndromeValue)
        SyndromeValue=Hlut[bitArray2Value(yError)]
#        print(Hlut[bitArray2Value(yError)])
        errorSyndromeString=errorSyndromeString+hex(SyndromeValue)[2]
#errorcorrect for hamming [7,4]
        if SyndromeValue!=0 and H.Ncol == 7:
           toggleM=[[0] for j in range(H.Ncol)]
           toggleM[SyndromeValue-1][0] = 1
           toggle=BinaryMatrix(toggleM)
           yError = toggle+yError
#errorcorrect for hamming [8,4]
        if SyndromeValue>7 and H.Ncol == 8:   #parity error (single bit errors)
           SyndromeValue = SyndromeValue-8
           toggleM=[[0] for j in range(H.Ncol)]
           toggleM[SyndromeValue-1][0] = 1
           toggle=BinaryMatrix(toggleM)
           yError = toggle+yError
        elif SyndromeValue>1 and H.Ncol == 8:   #parity error (double (or more) bit errors)
           toggleM=[[0] for j in range(H.Ncol)]
           toggleM[SyndromeValue-1][0] = 1
           toggle=BinaryMatrix(toggleM)
           yError = toggle+yError                  
#Decode
#        yDecode=R*yError
        yDecodeM=Rlut[bitArray2Value(yError)]
        yDecode=BinaryMatrix(yDecodeM)
        CorrectedDataString=CorrectedDataString+hex(bitArray2Value(yDecode))[2]
#        print(xM)
#        print(yDecode.M)
#        print(yDecodeM)
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
    duration = default_timer() - start
    print('yLS  ', yLSstring)
    print('yMS  ', yMSstring)
    print('#Err ', Errstring,' = ',f"{(100*NflippedBits/TotalNbits):.3g}",'%')
    print('ErSyn', errorSyndromeString)
    print('DetOk', DetectOkSting,' = ',f"{(100*NwrongDetects/Nnibbles):.3g}",'%')
    print('ycorr', CorrectedDataString)
    print('DatOk', DataCorrectString,' = ',f"{(100-100*NwrongReceived/Nnibbles):.3g}",'%')
    print('Time=', duration)
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
#AnalyseString(Istring,  4, GparM, HparM, RparM, ErrorFraction)
print('hamming7_4')
AnalyseString(Istring, 0, GHam74,HHam74 , RHam74, ErrorFraction)
#AnalyseString(Istring,  4, GHam74, HHam74, RHam74, ErrorFraction)
print('hamming8_4')
AnalyseString(Istring, 0, GHam84, HHam84, RHam84, ErrorFraction)
#AnalyseString(Istring,  4, GHam84, HHam84, RHam84, ErrorFraction)
print('hamming8_4 with lookup tables')
AnalyseStringLUTHamming84(Istring, 0, GHam84, HHam84, RHam84, ErrorFraction)