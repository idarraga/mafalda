import sys

prefix=sys.argv[1]
print("Merging = %s(8 chips)"%prefix)

# open all 8 files 
fd=[]
for chip in range(1,9) :
    fn = prefix+"_"+str(chip)+".txt"
    fd.append( open(fn) )

# assembly file
foutput = open(prefix+"_assembly.txt","w")

# Chips are placed like this
# 1 2 3 4
# 8 7 6 5

# Do first the top part
for j in range(0,256):
    fullline=""
    for chip in range(0,4): # indexes 0,1,2,3 in the list of descriptors
        line=fd[chip].readline()
        line=line.rstrip('\n')            # get rid of the '\n'
        fullline = fullline + " " + line  # an add the needed space
    #print fullline
    foutput.write(fullline+'\n')          # print to file now with the last '\n' once per line.

# Bottom now
for j in range(0,256):
    fullline=""
    for chip in range(7,3,-1):  # watch the order ! 8,7,6,5 --> indexes 7,6,5,4 in the list of descriptors
        line=fd[chip].readline()
        line=line.rstrip('\n')            # get rid of the '\n'
        fullline = fullline + " " + line  # an add the needed space
    #print fullline
    foutput.write(fullline+'\n')          # print to file now with the last '\n' once per line.

#foutput.write("test")

print("[PYTHON] done.");

