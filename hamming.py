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

GparM= [[1, 0, 0, 0],
[0, 1, 0, 0],
[0, 0, 1, 0],
[0, 0, 0, 1],
[1, 1, 1, 1]]
Gpar=BinaryMatrix(GparM)

HparM = [[1, 1, 1, 1, 1]]
Hpar=BinaryMatrix(HparM)

RparM= [[1, 0, 0, 0, 0],
[0, 1, 0, 0, 0],
[0, 0, 1, 0, 0],
[0, 0, 0, 1, 0]]
Rpar=BinaryMatrix(RparM)

LsHexstring=""
LsParstring=""
LsErrstring=""
MsHexstring=""
MsParstring=""
ParErrorDetectString=""
ParDetectOkSting=""
DataCorrectString=""

TotalNbits=0
NflippedBits=0
Nnibbles=0
NwrongDetects=0
NwrongReceived=0

ErrorFraction = float(input('Give fraction of bit errors (number between 0 and 1)'))
Istring = input('Text :')
for element in Istring:
    Nnibbles=Nnibbles+1
    asci=ord(element)
    asciHex=hex(asci)
    MsHexstring=MsHexstring+asciHex[2]
    LsHexstring=LsHexstring+asciHex[3]
#encode with parity
 #   print(hex(asci))
    xM=NibbleTo2DbitArray(asci)
 #   print(xM)
    x=BinaryMatrix(xM)
    y=Gpar*x
#print parity
    parHex=hex(y.M[4][0])                                 
    LsParstring=LsParstring+parHex[2]
#Introduce random Errors in data with parity
    ErrorVector=random.choices([1, 0], weights=[ErrorFraction, 1-ErrorFraction], k=len(y.M))
    NumberOfErrors=ErrorVector.count(1)
    NflippedBits=NflippedBits+NumberOfErrors
    TotalNbits=TotalNbits+len(y.M)
    LsErrstring=LsErrstring+str(NumberOfErrors)
    ErrorV=BinaryMatrix([ErrorVector])
    ErrorM=ErrorV.Mt
    Error=BinaryMatrix(ErrorM)
    yError =Error+y
#Error detect with parity
    ErrorDetect=Hpar*yError
    ParErrorDetectString=ParErrorDetectString+str(ErrorDetect.M[0][0])
#Decode with parity
    yDecode=Rpar*yError
    if ((yDecode.M!=xM and ErrorDetect.M[0][0]==0) or (yDecode.M==xM and ErrorDetect.M[0][0]!=0) ):
            ParDetectOkSting=ParDetectOkSting+'X'
            NwrongDetects=NwrongDetects+1
    else:
            ParDetectOkSting=ParDetectOkSting+'-'
    if (yDecode.M==xM):
            DataCorrectString=DataCorrectString+'-'
    else:
            DataCorrectString=DataCorrectString+'X'
            NwrongReceived=NwrongReceived+1
    
    m=NibbleTo2DbitArray(asci>>4)
    x=BinaryMatrix(m)
    y=Gpar*x
    parHex=hex(y.M[4][0])
    MsParstring=MsParstring+parHex[2]

print('LShex', LsHexstring)
print('LSpar', LsParstring)
print('#Err ', LsErrstring,' = ',f"{(100*NflippedBits/TotalNbits):.3g}",'%')
print('ErDet', ParErrorDetectString)
print('DetOk', ParDetectOkSting,' = ',f"{(100*NwrongDetects/Nnibbles):.3g}",'%')
print('DatOk', DataCorrectString,' = ',f"{(100*NwrongReceived/Nnibbles):.3g}",'%')
print('MShex', MsHexstring)
print('MSpar', MsParstring)



