# Hamming

Dit is de code van Julia Rossewij voor het Hamming practicum

Manual
The program encodes and decodes a asci tekst given by the user with 3 different coding schemes:
-	Parity
-	Hamming 74
-	Hamming 84
After encoding, random errors (bitflips) can be introduced with a probability given by the user. Execution time is also measured for each coding scheme. In addition hamming 84 is implemented with Look Up Tables (LUTs) instead of matrix multiplications.
Input
You can first enter the probability of bit errors occuring, this is a number greater or equal than 0 and smaller or equal than 1. Then enter the text you want to encode and encode. 
Output
-	yLS & yMS: Here the two most significant and the two least significant bits depending on the characters from your message will be shown.
-	#Err: #Err will show the number of errors were generated in each block.
-	ErSyn: This is the syndrome. For parity if 1, there is an error detected. For Haming, it indicates the position of the single bit error.
-	DetOk: If the error detection was a false positive (error detection while the input was equal to the output) or false negative (error detection was missed) than DetOK will show a X
-	Ycorr: Ycorr will show the received message
-	DatOk: If the decoded character is the different from the input character then DatOk displays a X.


