#!/usr/bin/python

#Psuedocode
#1. User need to choose either to input the sequence or load the sequence from file.
#2. User need to input minimum palindrome length from user, if no input from the user, the default value is 4.
#3. If user choose to load the sequence from file, prompt user to choose file format.
#4. Function is created to read and extract sequence in FASTA and GenBank files.
#5. Compute the reverse complement for the sequence
#6. Compute short sequences from the main sequence.
#7. Extract all the palindromes from the short sequences
#8. Compute classification of the palindromes into normal palindromes and spacer palindromes.
#9.    Normal Palindrome:
#       if the palindrome is common with its reverse complement sequence and in even number,
#      Spacer Palindrome:
#       if the palindrome is not exactly common with its reverse complement, contains '_' and in odd/even number.
#7. Extract and output normal gene palindromes
#8. Extract and output non-repeating palindromes with an intervening spacer region.

# Import modules for CGI handling

import cgi, cgitb
cgitb.enable(display=0, logdir="/path/to/logdir")


#Create instance of FieldStorage
form = cgi.FieldStorage()


######## GET VALUES FROM FIELD ###########
if form.getvalue('format'):
	format = form.getvalue('format')
else:
	format = "no Format"


if form.getvalue('rawSeq1'):
	if form.getvalue('rawSeq1'):
		seq = str(form.getvalue('rawSeq1'))
		seq = seq.upper() #set to uppercase
		seq = seq.strip() #remove any whitespaces
		format = "Raw"
	else:
		format = "invalid"

elif format == "Raw":

	if form.getvalue('rawSeq2'):
		seq = str(form.getvalue('rawSeq2'))
		seq = seq.upper() #set to uppercase
		seq = seq.strip() #remove any whitespaces
	else:
		format = "invalid"


elif format == "FASTA":
	if form.getvalue('rawSeq2'):
		content = form.getvalue('rawSeq2')
	else:
		format = "invalid"

elif format == "GenBank":
	if form.getvalue('rawSeq2'):
		content = form.getvalue('rawSeq2')
	else:
		format = "invalid"

else:
	format = "invalid"

if form.getvalue('minPal1'):
	if form.getvalue('minPal1'):
		minLength = int(form.getvalue('minPal1'))
 
	else:
		minLength = 4 #default value of palindrome

elif form.getvalue('minPal2'):
	if form.getvalue('minPal2'):
		minLength = int(form.getvalue('minPal2'))
 
	else:
		minLength = 4 #default value of palindrome

###### END OF GET VALUES FROM FIELD ########



####### PYTHON CODE ########

#import library 
import re 

# read in FASTA & GENBANK
def readFASTA(text_content):
    wo_space = re.sub('[\s]', '', text_content)
    seq_gb = '[ACGTacgt]{5}.*' # regular expression
    seqline = re.search(seq_gb, wo_space)

    if seqline:
        seq = seqline.group(0)[:20]  # truncate sequence to 20 nucleotide bases
        valid = True
    else:
        seq = ""
        valid = False

    return seq.upper(), valid  # set return string to uppercase


# read in GENBANK
def readGB(text_content):
    wo_space = re.sub('[\s]', '', text_content)
    seq_gb = '(ORIGIN)(.*)(\/\/$)' # regular expression
    seqline = re.search(seq_gb, wo_space)

    if seqline:
        seq = seqline.group(2)[:20]  # truncate sequence to 20 nucleotide bases
        seq = re.sub('[\d]', '', seq)
        valid = True
    else:
        seq = ""
        valid = False

    return seq.upper(), valid  # set return string to uppercase

 
#Compute the reverse complement seq
def reverseComplement(seq):
    Bases = {'A':'T','C':'G','G':'C','T':'A','_':'_'} #dictionary for base pairs
    comp = ''.join([Bases[i] for i in seq]) #concatenate all complement bases
    rcomp = comp[::-1]
    
    return rcomp  #return complement & reverse complement seq


#Compute short seqs from the main seq
def shortSeq(seq,minLength): 
    seqLength = len(seq) #compute the sequence length
    shortseq = [] #empty list to store all short seqs
    
    #extract short seq and store it in a list
    for i in range(seqLength,minLength-1,-1): #loop from the reverse of sequence until minimum palindrome length
        for k in range(seqLength-i+1): #loop in the length of short sequence
            shortseq.append(seq[k:i+k]) #append the short sequence into list

    return shortseq #return all the short seqs extracted from the sequence

#compute for all palindromes
def allPalindromes(shortseqs): 
    revcompsseq = "" #empty string to store reverse complement of short seq
    allPalindromeMatches = [] #empty list to store all palindromes matches

    for sseq in shortseqs: #for every shortseq in short seqs 
        revcompsseq = reverseComplement(sseq) #compute the reverse complement for shortseq
        
        #check if the first index of both shortseq and its reverse comp matches 
        #skip to the next iteration of shortseq if not match
        if (sseq[0] != revcompsseq[0]): 
            continue 

        #if the first and second index of both shortseq and its reverse comp matches, append it to the all palindrome matches list   
        elif (sseq[0] == revcompsseq[0] and sseq[1] == revcompsseq[1]):

            #to check if the short seq is already present in the list
            if sseq not in allPalindromeMatches: 
                allPalindromeMatches.append(sseq) #if not, append it to list 
            else: 
                break #break if short seq already present in list                                    
              
    if len(allPalindromeMatches): #true if list is not empty
        return(allPalindromeMatches) #return list that contains all palindrome matches
    else: 
        pass #pass if list is empty 

#compute classification of the palindromes into normal palindromes and spcaer palindromes
def classifyPalindrome(palindromes):
    normalP = [] #empty list to store normal palindromes
    spacerP = [] #empty list to store spacer palindromes
    cntnP = 0    #initiate count for normal palindromes
    cntsP = 0    #initiate count for spacer palindromes


    for palindrome in palindromes: 
        #if the palindrome is equal to its reverse complement seq, length is in even number, '_' absent in palindrome and palindrome not in normalPalindrome list
        if palindrome == reverseComplement(palindrome) and (len(palindrome) % 2) == 0 and '_' not in palindrome and palindrome not in normalP:
            normalP.append(palindrome) #append palindrome to nP list
        
        else:
            spacerP.append(palindrome) #else append palindrome to sP list
    
    #output for normal Palindrome
    if len(normalP): #returns true if normalP present
        print("<b>No. of Non-Spacer Palindromes found: %d </b><br><br>" % (len(normalP)))
        print("<table>")
        print("<tr><th><b>No</b></th><th><b>Palindromes with Non-Spacer Region</b></th></tr>")  # Adds column headers
        for n in normalP: 
            cntnP += 1
            print("<tr><td><b>%s </b></td>" % (cntnP)) 
            print("<td><b> %s </b></td></tr>" % (str(n)))
        print("</table>")

    else: 
        print("<b>There is no non-spacer palindromes found.</b><br>")

    
    #output for spacer Palindrome 
    if len(spacerP): #returns true if spacerP present
        print "<b><br><br><br>No. of Palindromes with Spacer Region found: %d </b><br><br>" % (len(spacerP))
        print("<table>")
        print("<tr><th><b>No</b></th><th><b>Palindromes with Spacer Region</b></th></tr>")  # Adds column headers
        for n in spacerP: 
            cntsP += 1
            print("<tr><td><b>%s </b></td>" % (cntsP)) 
            print("<td><b> %s </b></td></tr>" % (str(n)))
        print("</table>")

    else:
        print ("<b>There is no palindromes with an intervening spacer region found.</b><br>")


def result(seq, minLength):
       
     #call functions to get values
     rcomp = reverseComplement(seq) #call function reverseComplement()
     shortseqs = shortSeq(seq,minLength) #call function shortSeq()
     palindromes = allPalindromes(shortseqs) #call function allPalindromes()

     # print result 
     print "<p><br><b>Format Selected: %s</b></p><br>" % (format) 
     print "<p><b>Your input sequence is: %s </b></p><br>" % (seq)
     print "<p><b>Number of Base Pair: %d</b></p><br>" % (len(seq))
     
     print "<p><b>Minimum Palindrome Length: %d</b></p><br><br><br>" % (minLength)

     classifyPalindrome(palindromes) #call function classifyPalindrome()

####### END OF PYTHON CODE ########


####### WEBPAGE ########
print"Content-type:text/html\r\n\r\n";
print"<!DOCTYPE html>"
print"<html lang='en'>"
print"<head>"
print"<meta charset='UTF-8'>"
print"<meta http-equiv='X-UA-Compatible' content='IE=edge'>"
print"<meta name='viewport' content='width=device-width, initial-scale=1.0'>"
print"<title>PalindromeFinder</title>"
print"<link rel=stylesheet href=/siv3009/yoongsim/project/style.css>"
print"<link rel='icon' href='/siv3009/yoongsim/project/images/logo.png'>"
print"</head>"
print"<body>"
print"<div class=header>"
print"<img src=/siv3009/yoongsim/project/images/logo.png height=40 width=60>"
print"<h1>PalindromeFinder</h1>"
print"</div>"
 
print"<div id='header_nav'>"
print"<ul id='header_nav_ul'>"
print"<li id = 'header_li' ><a href='/siv3009/yoongsim/project/main.html'>Home</a></li>"
print"<li id = 'header_li' ><a href='/siv3009/yoongsim/project/seq.html' >Sequence Analysis</a></li>"
print"<li id = 'header_li'><a href='/siv3009/yoongsim/project/team.html' >About Us</a></li>"
print"</ul>"
print"</div>"
print"<b><b>"

print"<div>"
print"<section id='result'>"
print"<div style = 'text-align: center'>"
print"<br><br>"
print"<h2>Result of Sequence Analysis</h2>"
print"</div>"
print"<br><br>"
print"<table style='margin-left: auto; margin-right: auto; width: 90%;' border='0'>"
print"<tr><td>"


print"<fieldset>"


if format == "Raw":
     result(seq,minLength)
    
elif format == "FASTA":
     seq,valid = readFASTA(content)
     if valid:
          result(seq,minLength)
    
     else:
          print "<h2 style=color:red;>ERROR: please recheck input file</h2>"

elif format == "GenBank":
     seq,valid = readGB(content)
     if valid:
          result(seq,minLength)
     else:
          print "<h2 style=color:red;>ERROR: please recheck input file</h2>"

else:
     print "<h2 style=color:red;>ERROR: please recheck input file</h2>"
    

print"</fieldset>"
print"</td></tr>"
print"</table>"
print"</section>"
print"</div>"


print"<div class=footer>"
print"<div class=columns>"
print"<div class = 'columns small1'>"
print"<img src=/siv3009/yoongsim/project/images/logo.png height=50 width=75></img>"
print"</div >"

print"<div class= 'columns small' >"
print"<a><font color=#EDDC00 />Services</a><br>"
print"<b>Data resources and tools</b><br>"
print"<b>Data submission</b><br>"
print"<b>Support and feedback</b><br>"
print"<b>Licensing</b>"
print"</div>"
		
print"<div class= 'columns small' >"
print"<a><font color=#1A7F45 />Research</a><br>"
print"<b>Publications</b><br>"
print"<b>Reasearch groups</b><br>"
print"<b>Postdocs and PhDs</b>"
print"</div>"

print"<div class= 'columns small' >"
print"<a><font color=#2C70C6 />Training</a><br>"
print"<b>Live Training</b><br>"
print"<b>On-demand training</b><br>"
print"<b>Support for trainer</b><br>"
print"<b>Contact organisers</b>"
print"</div>"
		
print"<div class= 'columns small' >"
print"<a><font color=#FF0033 />Industry</b><br>"
print"<b>Members Area</b><br>"
print"<b>Contact Industry Team</b>"
print"</div>"
		
print"<div class= 'columns small' >"
print"<a><font color=#9D1027 />About PalindromeFinder</a><br>"
print"<b>Contact us</b><br>"
print"<b>Events</b><br>"
print"<b>Jobs</b><br>"
print"<b>News</b><br>"
print"<b>People and groups</b>"
print"</div>"
print"</div>"
print"</div>"

print"</body>"
print"</html>"
