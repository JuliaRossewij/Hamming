import random
import string
import numpy
import matplotlib.pyplot as plt
from re import X
from timeit import default_timer

# Class for (only 2D) binary Matrix with operator overloading for addition and multiplying with another Matrix 
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

# This function turns a nibble (least significant or most significant bit set) in a list representation
def NibbleTo2DbitArray(nible):
    M=[]
    for i in range(4):
        bit=[(nible&0x1)]
        M.append(bit)
        nible=nible>>1
    return M

# This function turns a byte in a list representation
def ByteTo2DbitArray(byte):
    M=[]
    for i in range(8):
        bit=[(byte&0x1)]
        M.append(bit)
        byte=byte>>1
    return M

# This function takes a Binary Matrix object as input. 
# It can only be used when it is a vector object.
# and converts the vector into a Value and returns this.
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
    sumSingelOperations=0

    # Start timing total time
    # Select most significant or least significant nibble depending on value ML_LS
    start = default_timer()
    for element in Istring:
        Nnibbles=Nnibbles+1
        asci=(ord(element)) >> MS_LS
        
# encode by matrix multiplication with generator matrix
# time the process of encoding the matrix
        xM=NibbleTo2DbitArray(asci)
        x=BinaryMatrix(xM)
        startSingleOperation = default_timer()
        y=G*x
        sumSingelOperations = sumSingelOperations + default_timer()-startSingleOperation  
#print encoded data
        sum1=bitArray2Value(y)
        yLSstring=yLSstring+f"{sum1:02x}"[1]
        yMSstring=yMSstring+f"{sum1:02x}"[0]   
        
# Introduce random Errors in data with parity
# by adding the error vector to the original parity/Hamming code
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
# start timing the process of the multiplication of y with the parity check matrix
# Calculate the SyndromeValue which denotes the place of the inverted bit
        startSingleOperation = default_timer()
        errorSyndrome=H*yError
        sumSingelOperations = sumSingelOperations + default_timer()-startSingleOperation 
        SyndromeValue=bitArray2Value(errorSyndrome)
        errorSyndromeString=errorSyndromeString+hex(SyndromeValue)[2]
#errorcorrect for hamming [7,4]
# adds a vector with 1 at the place of the denoted error to the Hammingcode
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
#  Time decoding process
#  Decode with matrix multiplication of the decoding matrix R
#  Print the string with the corrected data
        startSingleOperation = default_timer()
        yDecode=R*yError
        sumSingelOperations = sumSingelOperations + default_timer()-startSingleOperation 
        CorrectedDataString=CorrectedDataString+hex(bitArray2Value(yDecode))[2]

#evaluate
# Evaluate performance of error detection and correction
# by comparing original Hammingcode and received hammingcode and SyndromeValue
# print the total time and the joined time of the matrix multiplications
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
    print('TotalTime=', duration, "the matrix operations = ", sumSingelOperations)
    return duration, sumSingelOperations

#This function was added only to measure the time gain with lookup tables (LUTs) for hamming84. Is not tested for other cdong schemes
def AnalyseStringLUTHamming84(Istring, MS_LS, GM, HM, RM, ErrorFraction):
    G=BinaryMatrix(GM)
    H=BinaryMatrix(HM)
    R=BinaryMatrix(RM)
    # G,H and R LUTS are generated
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
    sumSingelOperations = 0
    start = default_timer()
    for element in Istring:
        Nnibbles=Nnibbles+1
        asci=(ord(element)) >> MS_LS
#encode 
        startSingleOperation = default_timer()
        xM=NibbleTo2DbitArray(asci)
        sumSingelOperations = sumSingelOperations + default_timer()-startSingleOperation
        y=BinaryMatrix(Glut[asci&0xF])
#print encoded data
        sum1=bitArray2Value(y)
        yLSstring=yLSstring+f"{sum1:02x}"[1]
        yMSstring=yMSstring+f"{sum1:02x}"[0]          
#Introduce random Errors in data with parity. 
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
        startSingleOperation = default_timer()
        SyndromeValue=Hlut[bitArray2Value(yError)]
        sumSingelOperations = sumSingelOperations + default_timer()-startSingleOperation
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
        startSingleOperation = default_timer()
        yDecodeM=Rlut[bitArray2Value(yError)]
        sumSingelOperations = sumSingelOperations + default_timer()-startSingleOperation
        yDecode=BinaryMatrix(yDecodeM)
        CorrectedDataString=CorrectedDataString+hex(bitArray2Value(yDecode))[2]

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
    print('TotalTime=', duration, "the lut operations = ", sumSingelOperations)
    return duration, sumSingelOperations

# the defenitions of the G, H and R matrices for parity and Hamming
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

# Start of the main programm

ErrorFraction = float(input('Give fraction of bit errors (number between 0 and 1)'))
x_axis=[]
LUToverallAVG=[]
LUToverallSTD=[]
MatrixOverallAVG=[]
MatrixOverallSTD=[]
LUTonlyAVG=[]
LUTonlySTD=[]
MatrixOnlyAVG=[]
MatrixOnlySTD=[]
ratioLUTmatrixOverall=[]
ratioLUTmatrixOnly=[]
for n in range(1, 200, 10):
    Istring = "".join(random.choice(string.ascii_lowercase) for i in range(n))
    print(Istring)
    LUToverallList=[]
    MatrixOverallList=[]
    LUTonlyList=[]
    MatrixOnlyList=[]
    for x in range(100):
        MatrixOverall, MatrixOnly = AnalyseString(Istring, 0, GHam84, HHam84, RHam84, ErrorFraction)
        LUToverall, LUTonly =AnalyseStringLUTHamming84(Istring, 0, GHam84, HHam84, RHam84, ErrorFraction)
        LUToverallList.append(LUToverall)
        MatrixOverallList.append(MatrixOverall)    
        LUTonlyList.append(LUTonly)
        MatrixOnlyList.append(MatrixOnly)    
    x_axis.append(n)
    LUToverallAVG.append(numpy.average(LUToverallList))
    LUToverallSTD.append(numpy.std(    LUToverallList))
    MatrixOverallAVG.append(numpy.average(MatrixOverallList))
    MatrixOverallSTD.append(numpy.std(    MatrixOverallList))
    LUTonlyAVG.append(numpy.average(LUTonlyList))
    LUTonlySTD.append(numpy.std(    LUTonlyList))
    MatrixOnlyAVG.append(numpy.average(MatrixOnlyList))
    MatrixOnlySTD.append(numpy.std(    MatrixOnlyList))
    ratioLUTmatrixOverall.append(100*numpy.average(LUToverallList)/numpy.average(MatrixOverallList))
    ratioLUTmatrixOnly.append(   100*numpy.average(LUTonlyList   )/numpy.average(MatrixOnlyList   ))


plt.errorbar(x_axis, MatrixOverallAVG, yerr=MatrixOverallSTD, label="Matrix overall")
plt.errorbar(x_axis, MatrixOnlyAVG   , yerr=MatrixOnlySTD   , label="Matrix only")
plt.errorbar(x_axis, LUToverallAVG   , yerr=LUToverallSTD   , label="LUT overall")
plt.errorbar(x_axis, LUTonlyAVG      , yerr=LUTonlySTD      , label="LUT only")
plt.xlabel("StringSize (=number of nibbles sent)")
plt.ylabel("ExecutionTime (s)")
plt.legend(loc="upper left")    
plt.show()

plt.plot(x_axis, ratioLUTmatrixOverall, label="ratio overall Execution time LUT / matrix ")
plt.plot(x_axis, ratioLUTmatrixOnly   , label="ratio Execution time LUT / matrix (calculation only)")
plt.xlabel("StringSize (=number of nibbles sent)")
plt.ylabel("ratio (%)")
plt.legend(loc='center')
plt.ylim((0,80))
plt.grid(True)
plt.show()
