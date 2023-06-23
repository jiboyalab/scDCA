import pandas as pd
import csv
import numpy as np
import argparse
from scipy.stats import gmean
def ReadMyCsv(SaveList, fileName):
    csv_reader = csv.reader(open(fileName))
    for row in csv_reader:  
        SaveList.append(row)
    return
def ReadMyTsv(SaveList, fileName):
    csv_reader = csv.reader(open(fileName),delimiter="\t")
    for row in csv_reader:  
        SaveList.append(row)
    return
def StorFile(data, fileName):
    with open(fileName, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)
    return
if __name__ =="__main__":
    parser = argparse.ArgumentParser(usage="it's usage tip.", description="Obtain the intersection of LR pairs output by 4 cellular communication tools, which are required to be found by at least 2 tools and have expression in single cell data.")
    parser.add_argument("--lr_cellphonedb", required=True, help="The results of LR pairs output by cellphonedb")
    parser.add_argument("--lr_cellchat", required=True, help="The results of LR pairs output by cellchat")
    parser.add_argument("--lr_nichenet", required=True, help="The results of LR pairs output by nichenet")
    parser.add_argument("--lr_icellnet", required=True, help="The results of LR pairs output by icellnet")
    parser.add_argument("--count", required=True, help="Count matrix / normalized count matrix path")
    parser.add_argument("--output", required=True, help="Output file")
    args = parser.parse_args()
    
    file1=args.lr_cellchat
    file2=args.lr_cellphonedb
    file3=args.lr_icellnet
    file4=args.lr_nichenet
    data1=[]
    ReadMyCsv(data1,file1)
    data2=[]
    ReadMyCsv(data2,file2)
    data3=[]
    ReadMyCsv(data3,file3)
    data4=[]
    ReadMyCsv(data4,file4)

    data=[]
    for i in range(1,len(data1)):
        data.append(data1[i][0])
    for i in range(1,len(data2)):
        data.append(data2[i][0])
    for i in range(1,len(data3)):
        data.append(data3[i][0])
    for i in range(1,len(data4)):
        data.append(data4[i][0])
    #de-duplication
    data_set=list(set(data))
    final=[]
    pair=[]
    pair.append("ligand")
    pair.append("receptor")
    final.append(pair)
    #Determine if the lr gene is in the expression
    expressiondata=pd.read_table(args.count,sep="\t",header=0,index_col=0)
    expressiondatalist=expressiondata.index.values.tolist()
    for i in range(len(data_set)):
        liglist=data_set[i].split("_")[0].split("+")
        reclist=data_set[i].split("_")[1].split("+")
        if set(liglist).issubset(set(expressiondatalist)) and set(reclist).issubset(set(expressiondatalist)):
            liglist_expresslist=[]
            for j in range(len(liglist)):
                liglist_expresslist.append(expressiondata.loc[liglist[j]].tolist())
            reclist_expresslist=[]
            for j in range(len(reclist)):
                reclist_expresslist.append(expressiondata.loc[reclist[j]].tolist())
            if np.all(gmean(liglist_expresslist,axis = 0) ==0) or np.all(gmean(reclist_expresslist,axis = 0) ==0):
                continue
            else:
                pair=[]
                pair.append(data_set[i].split("_")[0])
                pair.append(data_set[i].split("_")[1])
                final.append(pair)
    

    #Exists in at least two tools
    final2=[]
    pair=[]
    pair.append("ligand")
    pair.append("receptor")
    final2.append(pair)
    for i in range(1,len(final)):
        tep=[final[i][0]+"_"+final[i][1]]
        if tep in data1:
            label_data1=1
        else:
            label_data1=0
        if tep in data2:
            label_data2=1
        else:
            label_data2=0
        if tep in data3:
            label_data3=1
        else:
            label_data3=0
        if tep in data4:
            label_data4=1
        else:
            label_data4=0
        
        if label_data1+label_data2+label_data3+label_data4 >=2:
            pair=[]
            pair.append(final[i][0])
            pair.append(final[i][1])
            final2.append(pair)



    StorFile(final2,args.output)