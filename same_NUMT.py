#coding:utf-8
# by liuxuanzeng
# 2023-4-1
# 处理self-blastn后的结果,挑选物种内部相同的numt，判断插入后复制
# 处理slef-blastn后的结果,将相同的NUMT找到chr link
###输入1个
#1 物种名

import copy
# import pyfastx
import sys


specie_name = sys.argv[1]

n = ''
outputname0 = specie_name + '_same_numt.out'
output_file0 = open(outputname0, "w")

with open('self25_blast.out', 'r') as z:
    for line in z:
        lin = line.strip().split()
        nameL = lin[0]
        nameR = lin[1]
        if nameL == n:
            if seqlen - int(lin[3]) < 30 and int(lin[6]) < 15 and (int(lin[8]) < 15 or int(lin[9]) < 15):
                name2list = lin[1].split("MT")
                name2list2 = name2list[1].split('_')
                start = name2list2[0]
                end = name2list2[1]
                if ',' in start:
                    Llist = start.split(",")
                    start1 = Llist[0]
                    end1 = Llist[1]
                    Rlist = end.split(",")
                    start2 = Rlist[0]
                    end2 = Rlist[1]
                    destance = abs(int(start1) - int(end1)) + abs(int(start2) - int(end2))
                else:
                    destance = abs(int(start) - int(end))
                if abs(destance1 - destance) < 10 and float(lin[2]) >= 90:
                    if lin[0] != lin[1]:
                        for i in lin:
                            output_file0.write(i + "\t")
                        output_file0.write("\n")
        else:
            n = nameL
            seqlen = int(lin[3])
            name2list = lin[0].split("MT")
            name2list2 = name2list[1].split('_')
            start = name2list2[0]
            end = name2list2[1]
            if ',' in start:
                Llist = start.split(",")
                start1 = Llist[0]
                end1 = Llist[1]
                Rlist = end.split(",")
                start2 = Rlist[0]
                end2 = Rlist[1]
                destance1 = abs(int(start1) - int(end1)) + abs(int(start2) - int(end2))
            else:
                destance1 = abs(int(start) - int(end))
            if lin[0] != lin[1]:
                for i in lin:
                    output_file0.write(i + "\t")
                output_file0.write("\n")
z.close()
output_file0.close()

dt1 = {1: ['a', 'b']}
dt2 = {1: ['a', 'b']}
t = 1
with open(outputname0, 'r') as j:
    for line in j:
        lin = line.strip().split()
        nameL = lin[0]
        nameR = lin[1]
        n = 0
        for key in dt1:
            value = dt1[key]
            if nameL in value:
                if nameR in value:
                    n = 1
                else:
                    a = dt1[key]
                    b = a + [nameR]
                    dt2[key] = b
                    n = 1
            else:
                if nameR in value:
                    a = dt1[key]
                    b = a + [nameL]
                    dt2[key] = b
                    n = 1
        if n == 0:
            dt2[t] = [nameL, nameR]
            t = t + 1
        dt1 = copy.deepcopy(dt2)
j.close()

outputname1 = specie_name + '_same_numt_statistics.txt'

output_file1 = open(outputname1, "w")
for key in dt1:
    value = dt1[key]
    s = len(value)
    output_file1.write(str(s) + "\t")
    for i in value:
        output_file1.write(i + "\t")
    output_file1.write("\n")
output_file1.close()


inputblast = '../' + specie_name + '_NUMT_flank_all.txt'
inputsameblast = specie_name + '_same_numt_statistics.txt'
outputname3 = specie_name + '_same_NUMT_link.txt'

output_file3 = open(outputname3, "w")


dt = {}

with open(inputblast, 'r') as h:
    for line in h:
        lin = line.strip().split()
        chrname = lin[0]
        NUMT = lin[12]
        startpos = int(lin[6])
        endpos = int(lin[7])
        if 'seq' in chrname:
            turenamelist = chrname.split("_seq")
            turechr = turenamelist[0]
            addlist = turenamelist[1].split("_")
            add = int(addlist[0]) - 1
            startpos = startpos + add
            endpos = endpos + add
            chrname = turechr
        dt[NUMT] = (chrname, startpos, endpos)
h.close()

with open(inputsameblast, 'r') as z:
    for line in z:
        lin = line.strip().split()
        n = len(lin) - 1
        a = 2
        while a <= n:
            numt1 = lin[a-1]
            list1 = numt1.split("__w")
            numt1 = list1[0]
            numt2 = lin[a]
            list2= numt2.split("__w")
            numt2 = list2[0]
            chr1 = dt[numt1][0]
            start1 = dt[numt1][1]
            end1 = dt[numt1][2]
            chr2 = dt[numt2][0]
            start2 = dt[numt2][1]
            end2 = dt[numt2][2]
            a = a + 1
            output_file3.write(chr1 + '\t' + str(start1) + '\t' + str(end1) + '\t' + chr2 + '\t' + str(start2) + '\t' + str(end2) + '\n')
z.close()
output_file3.close()
print('Finish!')
