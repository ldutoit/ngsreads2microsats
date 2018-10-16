#!/usr/bin/env python

##WRITTEN BY LUDOVIC DUTOIT, OCTOBER 2018
#This file aims at cleaning the output of msatcommander by checking for overlapping sequences, joining files


####EXAMPLE USAGE
##simple:
#python2.7 clean_msatcomm_output.py genome.fasta

##extended:
#python2.7 clean_msatcomm_output.py  -microsat_file MICROSAT_FILE -primer_file PRIMER_FILE -output_folder OUTPUT_FOLDER genome.fasta




#usage: clean_msatcomm_output.py [-h] [-microsat_file MICROSAT_FILE]
#                                [-primer_file PRIMER_FILE]
#                                [-output_folder OUTPUT_FOLDER]
#                                genome#

#positional arguments:
#  genome                Fasta file used as input for msatcommander#

#optional arguments:
#  -h, --help            show this help message and exit
#  -microsat_file MICROSAT_FILE
#                        microsat file output of msatcommander
#                        default=msatcommander.microsatellites.csv
#  -primer_file PRIMER_FILE
#                        primer file output of msatcommander
#                        default=msatcommander.primers.csv
#  -output_folder OUTPUT_FOLDER
#                        output_folder


import os,itertools, argparse
from collections import Counter
import pandas as pd
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment 
import argparse



###PARAMETERS

parser = argparse.ArgumentParser() # add the parser
parser.add_argument("genome",help="Fasta file used as input for msatcommander") # add the parser
parser.add_argument("-microsat_file",help="microsat file output of msatcommander default=msatcommander.microsatellites.csv",type=str,default="msatcommander.microsatellites.csv") 
parser.add_argument("-primer_file",help="primer file output of msatcommander default=msatcommander.primers.csv",type=str,default="msatcommander.primers.csv") 
parser.add_argument("-output_folder",help="output_folder",type=str,default="primers") 
args = parser.parse_args()




###FUNCTIONS
def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])


##MAIN PROGRAM
# create folder if it doe not exits yet:

if os.path.exists(args.output_folder):
    raise Exception("The output folder specified already exists:\n \
     please specify another one or remove the currently existing one.")
else:
    os.mkdir(args.output_folder)

#DEAL WITH OUTPUT OF MSAT COMMANDER AND JOIN IT
#read files

print "Dealing with the output of msatcommander..."
microsat=pd.read_csv(args.microsat_file,sep=",")
primers=pd.read_csv(args.primer_file,sep=",",error_bad_lines=False)
print "Skipping two lines is normal"
assert primers.columns[2]=="msats_id"
primers =  primers.rename(columns={"msats_id":"id"})#

#merge the two files
joined =   pd.merge(primers, microsat, how='inner',on="id")
joined =  joined.rename(columns={"name_x":"name","records_id_x":"records_id"})#

# add the length of the columns
joined["length_motif"] =[len(a) for a in joined['motif']]
#remove duplicates
joined = joined[joined.duplicate == 0] 

##CREATE A TEMPORARY FILE WITH A SUBSET OF THE COLUMNS
joined_clean=joined[["name", "records_id","start","end" , "id","primer", "left_sequence", "right_sequence", "count", "motif", "length_motif"]].copy()
joined_clean.to_csv(args.output_folder+"/tempjoin.csv",index=False)



####PART 2, add to the formatted table sequence information

#read the genome

print "Appending sequences....\n...readind the genome"
input_file = open(args.genome)
my_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
output=open(args.output_folder+"/joined_microsat_withseq.csv","w")


#READ the microsat files and add the sequence info to the file
print "Adding the sequences"
dict_sequences={}
i=0
with open(args.output_folder+"/tempjoin.csv") as f:
    for line in f:
        i+=1
        if i==1:
            # order of the first four columns should be 'name', 'records_id', 'start', 'end', 'id', 'primer', 'left_sequence', 'right_sequence'
            assert line.split(",")[:8]==['name', 'records_id', 'start', 'end', 'id', 'primer', 'left_sequence', 'right_sequence']
            output.write(line.strip()+",lengthproduct,microsat_seq,primertoprimer_seq\n")
        else:
            msatid,start,end,seqname,leftseq,rightseq=line.split(",")[4],line.split(",")[2],line.split(",")[3],line.split(",")[0], line.split(",")[6], line.split(",")[7]
            microsat_seq = my_dict[seqname].seq[int(start):int(end)]
            ##check that the primers are only one in the sequence
            if my_dict[seqname].seq.count(leftseq) == my_dict[seqname].seq.count(reverse_complement(rightseq)) ==1:
                primertoprimer=my_dict[seqname].seq[my_dict[seqname].seq.find(leftseq):my_dict[seqname].seq.find(reverse_complement(rightseq))+len(rightseq)]
                length_product=len(primertoprimer)
                output.write(line.strip()+","+str(length_product)+","+str(microsat_seq)+","+str(primertoprimer)+"\n")
                dict_sequences["msatID_"+msatid+"_primertoprimer"]=str(primertoprimer)
                dict_sequences["msatID_"+msatid+"_leftprimer"]=str(leftseq)
                dict_sequences["msatID_"+msatid+"_rightprimer"]=str(rightseq)
            else:
                print "WARNING:the microsat id",line.split(",")[1],"is discared because one or both of the primers is/are found more than once in the contig"

output.close()


#3. Look for overlaps, if overlaps across different microsat_id
#produce alignments.log that contains all the alignments, the matching sequences ids and then remove all

print "...CHECKING  OVERLAP..."

output =open(args.output_folder+"/alignments.log","w")
scorethreshold =17
sequence_with_overlaps=set()
for key1,key2 in itertools.combinations(dict_sequences.keys(),2):
    if key1.split("_")[1]!=key2.split("_")[1]: # if they are different microsat only
        if (key1+key2).count("primertoprimer")<2:
            #print "makealign"
            alignments= pairwise2.align.localms(dict_sequences[key1],dict_sequences[key2],1,-1,-2,-1)
            alignments+= pairwise2.align.localms(dict_sequences[key1],reverse_complement(dict_sequences[key2]),1,-1,-2,-1)
            if any([alignment for alignment in alignments if alignment[2]>scorethreshold]):
                #print "WAAAAAA"
                sequence_with_overlaps.add(key1.split("_")[1])
                sequence_with_overlaps.add(key2.split("_")[1])
                print "OVERLAP FOUND, IGNORING:",key1,key2
                print(format_alignment(*[y for y in alignments if y[2]==max([x[2] for x in alignments])][0]))
                output.write("OVERLAP FOUND, IGNORING:"+key1+key2+"\n")
                output.write(format_alignment(*[y for y in alignments if y[2]==max([x[2] for x in alignments])][0])+"\n" )
output.write("\nDONE: Removed "+ str(len(sequence_with_overlaps))+"sequences:\n"+"\n".join([seq for seq in sequence_with_overlaps]))

output.close()
## 4 OUTPUTS

# REREAD MICROSAT_JOINED_SEQ and exclude overlaps
output=open(args.output_folder+"/joined_microsat_withseqNOOVERLAP.csv","w")
i=0
totalmsats = 0
with open(args.output_folder+"/joined_microsat_withseq.csv") as f:
    for line in f:
        i+=1
        if i==1: 
            assert line.split(",")[4] == "id"
            output.write(line)
        else:
            if not  line.split(",")[4] in sequence_with_overlaps:
                totalmsats +=1
                output.write(line)

output.close()

os.remove(args.output_folder+"/tempjoin.csv")

# sorting of microsat table
df =  pd.read_csv(args.output_folder+"/joined_microsat_withseqNOOVERLAP.csv")
df = df.sort_values(by=['length_motif',"count"],ascending=False)
df.to_csv(path_or_buf=args.output_folder+"/joined_microsat_withseqNOOVERLAP.csv", sep=',')

#make a fasta
output =open(args.output_folder+"joined_microsat_withseqNOOVERLAP.fasta","w")
for index, row in df.iterrows():
   output.write(">msatid_"+str(row["id"])+"_leftprimer\n"+row["left_sequence"]+"\n")
   output.write(">msatid_"+str(row["id"])+"_rightprimer\n"+row["right_sequence"]+"\n")
   output.write(">msatid_"+str(row["id"])+"_primertoprimer_seq\n"+row["primertoprimer_seq"]+"\n")

#FINAL MESSAGE
finalmessage = "SUMMARY\nAll alignments and overlaps reported in "+ args.output_folder+"/alignments.log\n\
All microsats before overlap filtering reported in "+ args.output_folder+"/joined_microsat_withseq.csv\n\
\nIDENTIFIED "+str(totalmsats)+" NON-OVERLAPPING MICROSATS\n\
Outputted as a fasta in "+ args.output_folder+"/joined_microsat_withseqNOOVERLAP.fasta\n\
and in "+ args.output_folder+"/joined_microsat_withseqNOOVERLAP.csv\nDONE"
output=open(args.output_folder+"/README.txt","w")
output.write(finalmessage)
output.close()
print "\n\n\n"+finalmessage
