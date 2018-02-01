#-*- coding: utf-8 -*-  
import sys
import random
import numpy as np

def countCG(strs):
    strs = strs.upper()
    return float((strs.count("C")+strs.count("G")))/(len(strs))

def count(strs):
    strs = strs.upper()
    return float(strs.count("A"))/len(strs),float(strs.count("C"))/len(strs),float(strs.count("G"))/len(strs),float(strs.count("T"))/len(strs)

def getString(fileStr):
    file = open(fileStr, 'r')
    gen_seq = ""
    lines = file.readlines()
    for line in lines:
        line = line.strip()
        gen_seq += line

    gen_seq = gen_seq.upper()
    return gen_seq

def get_reverse_str(str):
    str = str.upper()
    str_new=""
    for i in xrange(len(str)):
        if(str[i]=="T"):
            str_new+="A"
        elif(str[i]=="A"):
            str_new+="T"
        elif(str[i]=="G"):
            str_new+="C"
        elif(str[i]=="C"):
            str_new+="G"
        else:
            str_new+=str[i]
    return str_new

def getSubSeq(str, pos, K):
    n = len(str)
    l = pos - K
    r = pos + K + 1
    if l > r or l < 0 or r > n - 1:
        return 0

    elif "N" in str[l:r]:
        return 0

    return str[l:r]

def DNA2Sentence(dna, K):

    sentence = ""
    length = len(dna)

    for i in xrange(length - K + 1):
        sentence += dna[i: i + K] + " "

    #delete extra space
    sentence = sentence[0 : len(sentence) - 1]
    return sentence
