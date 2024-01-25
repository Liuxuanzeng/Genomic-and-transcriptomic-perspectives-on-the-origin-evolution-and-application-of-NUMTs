# coding:utf-8
# by liuxuanzeng
# 2023-8-29
# isoseq数据与线粒体比对后寻找最长转录本
###输入5个

#第1个 minimap _sort 结果文件
#需要将NUMTblast文件_NUMT_flank_all.txt拷贝到当前目录 ,'_mt_flank_100bp.fa 和  _genome_numt.fa拷贝到当前目录
#2 Sp name
# 3 比对百分比flort
#4 运行两天序列相差长度
#5 输入一个 iso测序文件

import pyfastx
import sys
dilen = int(sys.argv[4])
isofa = pyfastx.Fasta(sys.argv[5])

inputfile = sys.argv[1]
inputfile2 = sys.argv[2] + '_NUMT_flank_all.txt'


per = float(sys.argv[3])
outname1 = sys.argv[2] + '_nosame_MTtrans.txt'

output_file1 = open(outname1, "w")


n = 0

with open(inputfile, 'r') as f:
    for line in f:
        lin = line.strip().split()
        seqlen = lin[1]
        alignlen = lin[9]
        mtstart = lin[7]
        mtend = lin[8]
        seqname = lin[0]
        if int(alignlen)/int(seqlen) > per:
            n = n + 1
            if n == 1:
                seqlen2 = seqlen
                alignlen2 = alignlen
                mtstart2 = mtstart
                mtend2 = mtstart
                lin2 = lin
            else:
                if abs(int(mtstart2)-int(mtstart)) > 5 or abs(int(mtend) - int(mtend2)) > 5:
                    for z in lin2:
                        output_file1.write(z + "\t")
                    output_file1.write("\n")
                    seqlen2 = seqlen
                    alignlen2 = alignlen
                    mtstart2 = mtstart
                    mtend2 = mtend
                    lin2 = lin
for z in lin2:
    output_file1.write(z + "\t")
output_file1.write("\n")
f.close()
output_file1.close()

dt = {}
dt2 = {}
with open(outname1, 'r') as f:
    for line in f:
        lin = line.strip().split()
        seqlen = lin[1]
        mtstart = lin[7]
        mtend = lin[8]
        seqname = lin[0]
        dt[seqname] = (int(seqlen), int(mtstart), int(mtend))
f.close()

outname2 = sys.argv[2] + '_NUMT_mapto_MTtrans.txt'
output_file2 = open(outname2, "w")
outname3 = sys.argv[2] + '_NUMT_mapto_MTtrans.fa'
output_file3 = open(outname3, "w")
#输出map上NUMT侧翼100bp
flankname = sys.argv[2] + '_mt_flank_100bp.fa'
flankfa = pyfastx.Fasta(flankname)
outname4 = sys.argv[2] + '_NUMT_mapto_flank100.fa'
output_file4 = open(outname4, "w")
dt3 = {}
# 输出map上的NUMT
NUMTfname = sys.argv[2] + '_genome_numt.fa'
NUMTfa = pyfastx.Fasta(NUMTfname)
outname5 = sys.argv[2] + '_NUMT_mapto_numt.fa'
output_file5 = open(outname5, "w")

with open(inputfile2, 'r') as f:
    for line in f:
        lin = line.strip().split()
        NUMTlen = lin[3]
        NUMTstart = lin[8]
        NUMTend = lin[9]
        NUMTname = lin[12]
        if ',' in NUMTstart:
            pass
        else:
            if int(NUMTend) < int(NUMTstart):
                a = NUMTend
                b = NUMTstart
                NUMTstart = a
                NUMTend = b
            for i in dt:
                transname = i
                seqlen = dt[i][0]
                mtstart = dt[i][1]
                mtend = dt[i][2]
                if abs(int(NUMTstart) - mtstart) < dilen and abs(int(NUMTend) - mtend) < dilen:
                    if mtstart - int(NUMTstart) > 6 or int(NUMTend) - mtend > 6:
                        pass
                    else:
                        outlist = [NUMTname, NUMTlen, NUMTstart, NUMTend, transname, str(seqlen), str(mtstart), str(mtend)]
                        for k in outlist:
                            output_file2.write(k + "\t")
                        output_file2.write("\n")
                        if NUMTname not in dt3:
                            dt3[NUMTname] = ['1']
                            Lname = NUMTname + '__L100bp'
                            Rname = NUMTname + '__R100bp'
                            NUname = NUMTname + '__G'
                            se = NUMTfa[NUname]
                            output_file5.write('>' + se.name + "\n" + se.seq + "\n")
                            s = flankfa[Lname]
                            output_file4.write('>' + s.name + "\n" + s.seq + "\n")
                            s = flankfa[Rname]
                            output_file4.write('>' + s.name + "\n" + s.seq + "\n")
                        if transname not in dt2:
                            dt2[transname] = ['1']
                            s = isofa[transname]
                            output_file3.write('>' + s.name + "\n" + s.seq + "\n")

f.close()
output_file2.close()
output_file3.close()
output_file4.close()
output_file5.close()
