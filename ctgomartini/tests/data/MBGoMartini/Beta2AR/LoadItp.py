# 20230422
import os, re
import pandas

def SplitLines(lines,pattern='\[ (.+) \]'):
    startlist=[]
    endlist=[]
    for i,line in enumerate(lines):
        if len(re.findall(pattern, line))!=0:
            startlist.append(i)
            endlist.append(i-1)
    endlist.append(i)
    endlist=endlist[1:]

    Data_Dict={}
    for pair in zip(startlist,endlist):
        key=re.findall(pattern, lines[pair[0]])[0]

        value=lines[pair[0]+1:pair[1]+1]
        value=ProcessLines(value)
        if key in Data_Dict:
            Data_Dict[key]+=value
        else:
            Data_Dict[key]=value
    return Data_Dict


def ReadItp(file):
    with open(file,'r') as f:
        lines=f.readlines()
    newlines=[]
    for line in lines:
        line=line[:line.find(';')] if line.find(';') != -1 else line
        sline=line.strip()
        if sline=='': continue

        newlines.append(sline)

    Data_Dict=SplitLines(newlines,pattern='\[ (.+) \]')
    return Data_Dict
def ChangeTypes(alist):
    newlist=[]
    for item in alist:
        if re.match('\-?\d+\.\d',item):
            item=float(item)
        elif re.match('^\-?\d+$',item):
            item=int(item)
        else:
            item=item
        newlist.append(item)
    return newlist
def ProcessLines(lines):
    newlines=[]
    for line in lines:
        sline=line.strip().split()
        sline=ChangeTypes(sline)
        newlines.append(sline)
    return newlines


def WriteItp(itp_dict, outfile):
    lines=['\n']
    for key, values in itp_dict.items():
        lines+=[f'[ {key} ]\n']
        for value in values:
            line=' '.join(map(str, value))+'\n'
            lines+=[line]

        lines+=['\n']

    with open(outfile, 'w') as fp:
        fp.writelines(lines)
    return lines


if __name__ == "__main__":
    file='martini_v3.0.0.itp'
    Data_Dict=ReadItp(file)
