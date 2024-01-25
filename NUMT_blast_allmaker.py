#coding:utf-8
# by liuxuanzeng
# 2023-3-8
# 处理blastn后的结果
###输入4个
#第一个基因组fa文件
#2    线粒体fa文件
#3 线粒体gff文件
#4 物种名

import pyfastx
import sys

genome = sys.argv[1]
mtfa = sys.argv[2]
mtgff = sys.argv[3]
specie_name = sys.argv[4]
#############################################################################
Mfa = pyfastx.Fasta(mtfa)

inputname0 = specie_name + '_sort_E6.out'
outputname0 = specie_name + '_NUMT_cicos.txt'
output_file0 = open(outputname0, "w")

n = 0
a = 2
Gstat1, Gend1, Mstat1, Mend1 = 0, 0, 0, 0
mtlen = Mfa.size
print ('MT_lenth ' + str(mtlen))
with open(inputname0, 'r') as f:
    for line in f:
        lin2 = line.strip().split()
        Gstat2, Gend2 = int(lin2[6]), int(lin2[7])
        Mstat2, Mend2 = int(lin2[8]), int(lin2[9])
        if (Gstat2 - Gend1) < 3 and (((mtlen - Mstat2) < 3 and Mend1  < 3) or ( (mtlen - Mend1) < 3 and Mstat2 < 3)):
            lin1[7] = str(Gend2)
            lin1[8] = str(Mstat1) + ',' + str(Mend1)
            lin1[9] = str(Mstat2) + ',' + str(Mend2)
            a = 0
        if Gstat2 >= Gstat1 and Gend2 <=Gend1:
            Gstat2, Gend2 = Gstat1, Gend1
        else:
            if n > 0 and a != 1:
                for i in lin1:
                    output_file0.write(i + "\t")
                output_file0.write("\n")
        lin1 = lin2
        Gstat1, Gend1 = Gstat2, Gend2
        Mstat1, Mend1 = Mstat2, Mend2
        n = n + 1
        a = a + 1
    if a != 1:
        for i in lin1:
            output_file0.write(i + "\t")
        output_file0.write("\n")
print('Read ' + str(n) + ' lines OK!')

f.close()
output_file0.close()




####step1###  输出 name_NUMT_flank_all.txt
##产生一个包含NUMT坐标，numt名字，numt侧翼坐标的文件
#默认在Blast.out文件后添加列，第7、8列为基因组中numt起始,第9、10列为mt基因组的起始，第13列开始添加为numt_name
#第14、15、16、17，为30bp wiondow的起始1和起始2
#第18、19、20、21，为100bp wiondow的起始1和起始2
#第22、23、24、25，为1000bp wiondow的起始1和起始2
inputname1 = outputname0



outputname1 =   specie_name + '_NUMT_flank_all.txt'
output_file1 = open(outputname1, "w")
flankwindow1 = 25
flankwindow2 = 100
flankwindow3 = 2000
num = 1
with open(inputname1, 'r') as f:
    for line in f:
        lin = line.strip().split()
        numt = "numt" + str(num) + "_MT" + lin[8] + "_" + lin[9]
        lin.append(numt)
        # 30
        lin.append(str(int(lin[6]) - int(flankwindow1)))
        lin.append(str(int(lin[6]) - 1))
        lin.append(str(int(lin[7]) + 1))
        lin.append(str(int(lin[7]) + int(flankwindow1)))
        # 100
        lin.append(str(int(lin[6]) - int(flankwindow2)))
        lin.append(str(int(lin[6]) - 1))
        lin.append(str(int(lin[7]) + 1))
        lin.append(str(int(lin[7]) + int(flankwindow2)))
        # 1000
        lin.append(str(int(lin[6]) - int(flankwindow3)))
        lin.append(str(int(lin[6]) - 1))
        lin.append(str(int(lin[7]) + 1))
        lin.append(str(int(lin[7]) + int(flankwindow3)))
        num = num + 1
        for i in lin:
            output_file1.write(i+"\t")
        output_file1.write("\n")
#最后有个空行
f.close()
output_file1.close()
print('step1 OK!')
####step2###  输出 name_genome_numt.fa , name_mt_numt.fa , name_flank_25bp.fa , name_flank_100bp.fa, name_flank_2000bp.fa ,name_numt_with25bp.fa

Gfa = pyfastx.Fasta(genome)

outputname2 = specie_name + '_genome_numt.fa'
outputname3 = specie_name + '_mt_numt.fa'
outputname4 = specie_name + '_mt_flank_25bp.fa'
outputname5 = specie_name + '_mt_flank_100bp.fa'
outputname6 = specie_name + '_mt_flank_2000bp.fa'
outputname7 = specie_name + '_numt_with25bp.fa'

output_file2 = open(outputname2, "w")
output_file3 = open(outputname3, "w")
output_file4 = open(outputname4, "w")
output_file5 = open(outputname5, "w")
output_file6 = open(outputname6, "w")
output_file7 = open(outputname7, "w")

with open(outputname1, 'r') as z:
    for line in z:
        lin = line.strip().split()
        ###############################################################################
        # 输出 name_genome_numt.fa
        chrname = lin[0]
        statpos = lin[6]
        endpos = lin[7]
        newseqname = '>'+ lin[12] + '__G'
        interval = (int(statpos), int(endpos))
        seq = Gfa.fetch(chrname, interval)      #srtand
        outseqname = [newseqname, chrname, statpos, endpos]
        for i in outseqname:
            output_file2.write(i + " ")
        output_file2.write("\n")
        output_file2.write(seq + "\n")
        #################################################################################
        # 输出 name_mt_numt.fa
        chrname = lin[1]
        if ',' in lin[8]:
            MT1 = lin[8].split(",")
            MT2 = lin[9].split(",")
            statpos = MT1[0]
            endpos = MT1[1]
            if int(statpos) > int(endpos):
                statpos = MT1[1]
                endpos = MT1[0]
            newseqname = '>' + lin[12] + '__M'
            interval = (int(statpos), int(endpos))
            seq = Mfa.fetch(chrname, interval)  # srtand
            outseqname = [newseqname, chrname, statpos, endpos]
            for i in outseqname:
                output_file3.write(i + " ")
            output_file3.write("\n" + seq + "\n")
        else:
            statpos = lin[8]
            endpos = lin[9]
            if int(statpos) > int(endpos):
                statpos = lin[9]
                endpos = lin[8]
            newseqname = '>' + lin[12] + '__M'
            interval = (int(statpos), int(endpos))
            seq = Mfa.fetch(chrname, interval)  # srtand
            outseqname = [newseqname, chrname, statpos, endpos]
            for i in outseqname:
                output_file3.write(i + " ")
            output_file3.write("\n" + seq + "\n")

        ########################################################################
        # 输出 name_flank_25bp.fa
        chrname = lin[0]
        # 输出左边25bp序列
        statpos = int(lin[13])
        endpos = int(lin[14])
        newseqname = '>' + lin[12] + '__L25bp'
        if statpos < 1:
            statpos = 1
            endpos = endpos + 2
        s1 = Gfa[chrname]
        lenth = len(s1)
        if  endpos > lenth:
            endpos = lenth
            statpos = endpos - 1
        interval = (statpos, endpos)
        seq = Gfa.fetch(chrname, interval)  # srtand
        outseqname = [newseqname, chrname, str(statpos), str(endpos)]
        for i in outseqname:
            output_file4.write(i + " ")
        output_file4.write("\n")
        output_file4.write(seq + "\n")
        # 输出右边25bp序列
        statpos = int(lin[15])
        endpos = int(lin[16])
        newseqname = '>' + lin[12] + '__R25bp'
        if statpos < 1:
            statpos = 1
            endpos = endpos + 2
        s1 = Gfa[chrname]
        lenth = len(s1)
        if endpos > lenth:
            endpos = lenth
            statpos = endpos - 1
        interval = (statpos, endpos)
        seq = Gfa.fetch(chrname, interval)  # srtand
        outseqname = [newseqname, chrname, str(statpos), str(endpos)]
        for i in outseqname:
            output_file4.write(i + " ")
        output_file4.write("\n")
        output_file4.write(seq + "\n")
        ########################################################################
        # 输出 name_flank_100bp.fa
        chrname = lin[0]
        # 输出左边100bp序列
        statpos = int(lin[17])
        endpos = int(lin[18])
        newseqname = '>' + lin[12] + '__L100bp'
        if statpos < 1:
            statpos = 1
            endpos = endpos + 2
        s1 = Gfa[chrname]
        lenth = len(s1)
        if endpos > lenth:
            endpos = lenth
            statpos = endpos - 1
        interval = (statpos, endpos)
        seq = Gfa.fetch(chrname, interval)  # srtand
        outseqname = [newseqname, chrname, str(statpos), str(endpos)]
        for i in outseqname:
            output_file5.write(i + " ")
        output_file5.write("\n")
        output_file5.write(seq + "\n")
        # 输出右边100bp序列
        statpos = int(lin[19])
        endpos = int(lin[20])
        newseqname = '>' + lin[12] + '__R100bp'
        if statpos < 1:
            statpos = 1
            endpos = endpos + 2
        s1 = Gfa[chrname]
        lenth = len(s1)
        if endpos > lenth:
            endpos = lenth
            statpos = endpos - 1
        interval = (statpos, endpos)
        seq = Gfa.fetch(chrname, interval)  # srtand
        outseqname = [newseqname, chrname, str(statpos), str(endpos)]
        for i in outseqname:
            output_file5.write(i + " ")
        output_file5.write("\n")
        output_file5.write(seq + "\n")
        ########################################################################
        # 输出 name_flank_2000bp.fa
        chrname = lin[0]
        # 输出左边2000bp序列
        statpos = int(lin[21])
        endpos = int(lin[22])
        newseqname = '>' + lin[12] + '__L2000bp'
        if statpos < 1:
            statpos = 1
            endpos = endpos + 2
        s1 = Gfa[chrname]
        lenth = len(s1)
        if endpos > lenth:
            endpos = lenth
            statpos = endpos - 1
        interval = (statpos, endpos)
        seq = Gfa.fetch(chrname, interval)  # srtand
        outseqname = [newseqname, chrname, str(statpos), str(endpos)]
        for i in outseqname:
            output_file6.write(i + " ")
        output_file6.write("\n")
        output_file6.write(seq + "\n")
        # 输出右边2000bp序列
        statpos = int(lin[23])
        endpos = int(lin[24])
        newseqname = '>' + lin[12] + '__R2000bp'
        if statpos < 1:
            statpos = 1
            endpos = endpos + 2
        s1 = Gfa[chrname]
        lenth = len(s1)
        if endpos > lenth:
            endpos = lenth
            statpos = endpos - 1
        interval = (statpos, endpos)
        seq = Gfa.fetch(chrname, interval)  # srtand
        outseqname = [newseqname, chrname, str(statpos), str(endpos)]
        for i in outseqname:
            output_file6.write(i + " ")
        output_file6.write("\n")
        output_file6.write(seq + "\n")
 ########################################################################
        # 输出 name_numt_with25bp.fa
        chrname = lin[0]
        statpos = int(lin[13])
        endpos = int(lin[16])
        newseqname = '>' + lin[12] + '__with25bp'
        if statpos < 1:
            statpos = 1
            endpos = endpos + 2
        s1 = Gfa[chrname]
        lenth = len(s1)
        if endpos > lenth:
            endpos = lenth
            statpos = endpos - 1
        interval = (statpos, endpos)
        seq = Gfa.fetch(chrname, interval)  # srtand
        outseqname = [newseqname, chrname, str(statpos), str(endpos)]
        for i in outseqname:
            output_file7.write(i + " ")
        output_file7.write("\n")
        output_file7.write(seq + "\n")
z.close()

output_file2.close()
output_file3.close()
output_file4.close()
output_file5.close()
output_file6.close()
output_file7.close()
print('fa extract OK!')
########################################################################
#做一些序列统计，输出 name_numt_statistics.txt
#包括第一部分GC含量的统计

outputname8 = specie_name + '_numt_statistics.txt'

output_file8 = open(outputname8, "w")

Gnumt = pyfastx.Fasta(outputname2)
F25 = pyfastx.Fasta(outputname4)
F100 = pyfastx.Fasta(outputname5)

output_file8.write('###' + outputname8 + '\n')


output_file8.write('###GC含量统计###'+ '\n')
output_file8.write (specie_name + '\t' + 'genome_GC_content' + '\t' + str(Gfa.gc_content) + '\n')
output_file8.write (specie_name + '\t' + 'genome_GC_skew' + '\t' + str(Gfa.gc_skew) + '\n')
output_file8.write (specie_name + '\t' + 'numt_GC_content' + '\t' + str(Gnumt.gc_content) + '\n')
output_file8.write (specie_name + '\t' + 'numt_GC_skew' + '\t' + str(Gnumt.gc_skew) + '\n')
output_file8.write (specie_name + '\t' + 'Flank25bp_GC_content' + '\t' + str(F25.gc_content) + '\n')
output_file8.write (specie_name + '\t' + 'Flank25bp_GC_skew' + '\t' + str(F25.gc_skew) + '\n')
output_file8.write (specie_name + '\t' + 'Flank100bp_GC_content' + '\t' + str(F100.gc_content) + '\n')
output_file8.write (specie_name + '\t' + 'Flank100bp_GC_skew' + '\t' + str(F100.gc_skew) + '\n')
output_file8.write('###  NUMT数量和总长度统计 ###'+ '\n')
output_file8.write (specie_name + '\t' + 'Numt_count' + '\t' + str(len(Gnumt)) + '\n')
output_file8.write (specie_name + '\t' + 'Numt_length' + '\t' + str(Gnumt.size) + '\n')
output_file8.write (specie_name + '\t' + 'MT_length' + '\t' + str(mtlen) + '\n')

###
###长度统计区间0-200bp ，200-400，400-600bp，800-1000，1000-1500，1500-2000，2000以上
regeion1 = 0
regeion2 = 0
regeion3 = 0
regeion4 = 0
regeion5 = 0
regeion6 = 0
regeion7 = 0
regeion8 =0
output_file8.write('###numt序列长度统计###'+ '\n')

for s in Gnumt:
    length = len(s)
    if length <= 200 :
        regeion1 =  regeion1 + 1
    elif length <= 400 :
        regeion2 = regeion2 + 1
    elif length <= 600 :
        regeion3 = regeion3 + 1
    elif length <= 800 :
        regeion4 = regeion4 + 1
    elif length <= 1000 :
        regeion5 = regeion5 + 1
    elif length <= 1500 :
        regeion6 = regeion6 + 1
    elif length <= 2000:
        regeion7 = regeion7 + 1
    else:
        regeion8 = regeion8 + 1
output_file8.write('name'+ '\t' + '<=200bp' + '\t'+ '201-400bp' + '\t'+ '401-600bp' + '\t'+ '601-800bp'  + '\t'+ '801-1000bp' + '\t'+ '1001-1500bp' + '\t'+ '1501-2000bp'+ '\t'+ '>2000bp' + '\n')
output_file8.write(specie_name+ '\t' + str(regeion1) + '\t' + str(regeion2)  + '\t' + str(regeion3)  + '\t'+ str(regeion4) + '\t'+ str(regeion5) + '\t'+ str(regeion6) + '\t'+ str(regeion7) + '\t'+ str(regeion8) +  '\n')


output_file8.write('###numt编码蛋白统计###'+ '\n')

dt = {}
num = {}
with open(mtgff, 'r') as h:
    for line in h:
        lin = line.strip().split()
        if len(lin) > 5 and lin[2] == 'gene' and ('gene=' in lin[8]):
                    list2 = lin[8].split("gene=")
                    if ';' in list2[1]:
                        namelist = list2[1].split(";")
                        genename = namelist[0]
                        startpos = int(lin[3])
                        endpos = int(lin[4])
                        dt[genename] = [startpos, endpos]
                    else :
                        genename = list2[1]
                        startpos = int(lin[3])
                        endpos = int(lin[4])
                        dt[genename] = [startpos,endpos]
h.close()

with open(outputname1, 'r') as z:
    for line in z:
        lin = line.strip().split()
        if ',' in lin[8]:
            MT1 = lin[8].split(",")
            MT2 = lin[9].split(",")
            stat = int(MT1[0])
            end = int(MT1[1])
            if star > end:
                star = int(MT1[1])
                end = int(MT1[0])
            for i in dt:
                genename = i
                genestart = dt[i][0]
                geneend = dt[i][1]
                if (star < int(genestart) and (end - int(genestart)) > 100) or (star >= int(genestart) and ( int(geneend) - star) > 100):
                    if genename in num:
                        num[genename] = num[genename] + 1
                    else:
                        num[genename] = 1
            ###
            stat = int(MT2[0])
            end = int(MT2[1])
            if star > end:
                star = int(MT2[1])
                end = int(MT2[0])
            for i in dt:
                genename = i
                genestart = dt[i][0]
                geneend = dt[i][1]
                if (star < int(genestart) and (end - int(genestart)) > 100) or (star >= int(genestart) and ( int(geneend) - star) > 100):
                    if genename in num:
                        num[genename] = num[genename] + 1
                    else:
                        num[genename] = 1
        else:
            star = int(lin[8])
            end = int(lin[9])
            if star > end:
                star = int(lin[9])
                end = int(lin[8])
            for i in dt:
                genename = i
                genestart = dt[i][0]
                geneend = dt[i][1]
                if (star < int(genestart) and (end - int(genestart)) > 100) or (star >= int(genestart) and ( int(geneend) - star) > 100):
                    if genename in num:
                        num[genename] = num[genename] + 1
                    else:
                        num[genename] = 1
for n in num:
    output_file8.write(n + '\t' + str(num[n]) + '\n' )

z.close()
output_file8.close()
print('Finish!')

