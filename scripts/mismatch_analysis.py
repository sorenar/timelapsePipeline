
from optparse import OptionParser
import re



def getOptions():
    parser = OptionParser()
    parser.add_option("--input","-i", dest = "infile", help = "Input filtered sam file",
                      metavar = "FILE", type = "string", default = "")

    parser.add_option("--output_pre","-o", dest = "outfile", help = "Output prefix for sub samfiles",
                      metavar = "FILE", type = "string", default = "")

    parser.add_option("--SNPs","-s", dest = "SNP", help = "SNP file made from controls",
                      metavar = "FILE", type = "string", default = "")

    parser.add_option("--chr_name", "-c", dest = "chr_name",
                      help = "The chromosome ot be processed",default = "")

    parser.add_option("--chr_length", "-l", dest = "chr_len",
                      help = "The chromosome ot be processed",default = 0)

    (options, args) = parser.parse_args()
    return options


def main():
    options = getOptions()
    sam_file = options.infile
    out_prefix = options.outfile
    chr_name = options.chr_name
    SNP_file = options.SNP
    chr_len = options.chr_len
    
   # print("inputs are: ", sam_file,"\n",out_prefix,"\n",SNP_file,"\n" )

    header = []
    #chroms = set()
    transcripts_0TC = []
    transcripts_1TC = []
    transcripts_2TC = []
    transcripts_3TC = []
    transcripts_4TC = []
    transcripts_5pTC = []
    
    # making a list of T to C SNPs
    snp = []
    f = open(SNP_file,"r")
    for line in f:
       line = line.strip()
       snp.append("_".join(line.split("\t")[0:2]))
    f.close()
    
    #print("first three SNP loci: \n",snp[1:3])
    f_log = open(str(out_prefix + "_log.txt"),"a+")
    # line_counter = 0
    with open(sam_file, 'r') as f:
       for line in f:
    #     print("\n",p)
          line = line.strip()
          if line.startswith('@'):
             header.append(line)
             # line_counter+=1
          else:
             # line_counter+=1
             # print("line #:   ",line_counter)
             # chrom = line.split("\t")[2]
             # chroms.add(chrom)
             fields = line.split("\t")
             cigar = fields[5]
             seq = fields[9]
             qual = fields[10]
             name = fields[0]
             chr = fields[2]
             pos = fields[3]
             md_field = [m for m in fields[11:len(fields)] if m.startswith("MD")]
             if len(md_field) == 0:
                f_log.write("No MD field \n") 
                continue
             md = md_field[0].split(":")[2]
             md_ops = [x for x in re.split("[^a-zA-Z]*",md) if x]
             md_counts = [int(x) for x in re.split("[^0-9]*",md) if x]
             cigar_ops = [x for x in re.split("[^a-zA-Z]*",cigar) if x]
             cigar_counts = [int(x) for x in re.split("[^0-9]*",cigar) if x]        
             q_ind = 0
             r_ind = 0
             md_temp = md
             TCn=0
             SNPn=0
             f_log.write(str(fields) + "\n")
             # print("cigar = ",cigar,"md = ",md)
             for i in range(len(cigar_ops)):
                m = cigar_ops[i]
                n = cigar_counts[i]
                #f_log.write("n = " + str(n))
                if m == "M":
                   md_ind = 0
                   while md_ind < n:
                      token = re.search("^[0-9]*",md_temp).group(0)
                      if token != '':
                         # it's a match
                         md_temp = re.sub("^[0-9]*",'',md_temp)
                         if (md_ind+int(token)) > n:
                            md_temp = str(md_ind + int(token) - n) + md_temp
                            q_ind+=(n - md_ind)
                            r_ind+=(n - md_ind)
                            md_ind = n
                         else:              
                            q_ind+=int(token)
                            r_ind+=int(token)
                            md_ind+=int(token)
                         f_log.write("md_temp = "+ str(md_temp)+ "md_ind = "+ str(md_ind)+ "n = " + str(n)+"\n")       
                      else:
                         token = re.search("^[ATCGNnatcg]*",md_temp).group(0)
                         md_temp = re.sub("^[ATCGNnatcg]*",'',md_temp)
                         f_log.write("md_temp = "+ str(md_temp)+ "md_ind = "+ str(md_ind)+ "n = "+ str(n))
                         if len(token) > 0:
                            for j in range(len(token)):
                               # it's a mismatch
                               q_ind+=1
                               md_ind+=1
                               r_ind+=1
                               s = seq[q_ind-1]
                               r = token[j]
                               q = ord(qual[q_ind-1])-33
                               f_log.write("r = "+str(r)+"; s = "+str(s)+"; q = "+str(q)+"\n")
                               chr_pos = chr + "_" + str(int(pos)+r_ind-1)
                               if r == "T" and s == "C" and q >30 and not chr_pos in snp:
                                  TCn+=1
                               elif q > 30:
                                  SNPn+=1
                elif m == "D":
                   r_ind+=n
                   token = re.search("^\^[ATCGNnatcg]*",md_temp).group(0)
                   md_temp = md_temp[n+1:len(md_temp)]
                elif m == "S":
                   q_ind+=n
                elif m == "I":
                   q_ind+=n
                elif m == "N":
                   r_ind+=n
             if TCn == 0:
                transcripts_0TC.append(line)
             elif TCn == 1:
                transcripts_1TC.append(line)
             elif TCn == 2:
                transcripts_2TC.append(line)
             elif TCn == 3:
                transcripts_3TC.append(line)
             elif TCn == 4:
                transcripts_4TC.append(line)
             else:
                transcripts_5pTC.append(line)

    f.close()
    f_log.close()
    with open(str(out_prefix + "_py_0.sam"),"a+") as f:
         for i in range(len(header)):
             f.write(header[i]+"\n")
         for i in range(len(transcripts_0TC)):
             f.write(transcripts_0TC[i]+"\n")
    f.close()

    with open(str(out_prefix + "_py_1.sam"),"a+") as f:
         for i in range(len(header)):
             f.write(header[i]+"\n")
         for i in range(len(transcripts_1TC)):
             f.write(transcripts_1TC[i]+"\n")
    f.close()

    with open(str(out_prefix + "_py_2.sam"),"a+") as f:
         for i in range(len(header)):
             f.write(header[i]+"\n")
         for i in range(len(transcripts_2TC)):
             f.write(transcripts_2TC[i]+"\n")
    f.close()

    with open(str(out_prefix + "_py_3.sam"),"a+") as f:
         for i in range(len(header)):
             f.write(header[i]+"\n")
         for i in range(len(transcripts_3TC)):
             f.write(transcripts_3TC[i]+"\n")
    f.close()

    with open(str(out_prefix + "_py_4.sam"),"a+") as f:
         for i in range(len(header)):
             f.write(header[i]+"\n")
         for i in range(len(transcripts_4TC)):
             f.write(transcripts_4TC[i]+"\n")
    f.close()

    with open(str(out_prefix + "_py_5p.sam"),"a+") as f:
         for i in range(len(header)):
             f.write(header[i]+"\n")
         for i in range(len(transcripts_5pTC)):
             f.write(transcripts_5pTC[i]+"\n")
    f.close()
    


if __name__ == '__main__':

   main()




