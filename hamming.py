s = input ()
string = ""


for i in s :
    binary = str(bin(ord(i)))[:5] + " " + str(bin(ord(i)))[5:] + " "
    binary = binary.replace("b","")
    string  = string + binary

print (string)
lijst = string.split(" ")
