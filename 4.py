# -*- coding: utf-8 -*-
# @Author: Liu Yang
# @Date:   2023-03-10 20:30:48
# @Last Modified by:   Liu Yang
# @Last Modified time: 2023-03-14 14:46:44


###################################################################################################################################3
#处理Reviewed数据





import os
os.chdir('D://abishe//Fifth//Reviewed')
import pandas as pd

#####################################################################################################
#处理突变格式
data1=pd.read_csv('Reviewed_shaixuan3.tsv',header=0,sep='\t')
mut=[]
for i in range(len(data1)):
	a1=str(data1['Feature range(s)'][i]).split('-')[0]
	a2=str(data1['Original sequence'][i])
	a3=str(data1['Resulting sequence'][i])
	a4=a2+a1+a3
	mut.append(a4)
mut_data=pd.DataFrame({'Mutations':mut})
b=data1['type']
c=data1['Affected protein AC']
d=data1['Interaction AC']
e=data1['Partners']
f=data1['Partners_1']
g=data1['Feature range(s)']
h=data1['Original sequence']
data2=pd.concat([d,c,mut_data,e,f,b,g,h],axis=1)
data2.to_csv('Reviewed.tsv',index=False,sep='\t')#39451


###############################################################################################################################################
#获取序列特征

#Reviewed
a=pd.read_csv('Reviewed.tsv',header=0,sep='\t')
#读取fasta格式的uniprot现有且审查蛋白数据,并进行处理
f=open('uniprot_sprot.fasta')
seq={}
for line in f:
	if line.startswith('>'):
		name=line.replace('\n','').split()[0]
		seq[name]=''
	else:
		seq[name]+=line.replace('\n','').strip()
P_ID=[]
for i in seq.keys():
	start=i.index('|')
	end=i.index('|',i.index('|')+1)
	P_ID.append(i[start+1:end])
P_seq=[]
for j in seq.values():
	P_seq.append(j)
data_uniprot={'Protein_ID':P_ID,'Protein_seq':P_seq}
Uniprot=pd.DataFrame(data_uniprot)
f.close()

seqs1=[]
seqs2=[]
for i in range(len(a)):
	temp1=eval(a['Partners_1'][i])[0][10:]
	temp2=eval(a['Partners_1'][i])[1][10:]
	seq1=list(Uniprot[Uniprot['Protein_ID']==temp1]['Protein_seq'])[0]
	seq2=list(Uniprot[Uniprot['Protein_ID']==temp2]['Protein_seq'])[0]
	seqs1.append(seq1)
	seqs2.append(seq2)
seq1_dat=pd.DataFrame({'A sequence':seqs1})
seq2_dat=pd.DataFrame({'B sequence':seqs2})
data=pd.concat([a,seq1_dat,seq2_dat],axis=1)
data.to_csv('Reviewed_seq.tsv',index=False,sep='\t')




###############################################
#检查突变位点与原序列一致性
import re
data=pd.read_csv('Reviewed_seq.tsv',header=0,sep='\t')


index_cuo=[]
for i in range(len(data)):
	a1=int(re.findall(r'\d+',data['Mutations'][i])[0])
	a2=data['Mutations'][i][0]
	a3=data['Affected protein AC'][i]
	a4=eval(data['Partners_1'][i])
	a_index=a4.index(a3)
	if a_index==0:
		temp1=data['A sequence'][i]
	else:
		temp1=data['B sequence'][i]

	if temp1[a1-1] !=a2:
		index_cuo.append(i)

data_1=data.loc[index_cuo]
data1=data.drop(index_cuo)
data_1.to_csv('Reviewed_cuowu.tsv',index=False,sep='\t')
data1.to_csv('Reviewed_zhengque.tsv',index=False,sep='\t')



#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#去除冲突数据

#Reviewed
import os
os.chdir('D://abishe//Fifth//Reviewed')
import pandas as pd

data=pd.read_csv('Reviewed_zhengque.tsv',header=0,sep='\t')
data1=list(data.groupby(['Affected protein AC','Mutations','Interaction AC']))
index2=[]
for j in range(len(data1)):
	if len(data1[j][1])==2:
		index2.append(j)

data2=pd.DataFrame()
for i in index2:
	x=data1[i][1]
	data2=pd.concat([data2,x])

qu_index=list(data2.index)


##########################
data1=list(data.groupby(['Affected protein AC','Mutations','Interaction AC','type']))
index2=[]
for j in range(len(data1)):
	if len(data1[j][1])>1:
		index2.append(j)

data2=pd.DataFrame()
for i in index2:
	x=data1[i][1]
	data2=pd.concat([data2,x])

qu_index1=list(data2.index)

##############
index1=list(set(qu_index)-set(qu_index1))
data_chuli=data.loc[index1]
data_result1=data.drop(index1)

data_result1.to_csv('Reviewed_last1.tsv',index=False,sep='\t')
data_chuli.to_csv('Reviewed_data_chuli.tsv',index=False,sep='\t')



a=pd.read_csv('Reviewed_data_chuli.tsv',header=0,sep='\t')
a1=list(a.groupby(['Interaction AC','Affected protein AC','Mutations']))
chongtu_type1=['Decreasing','Increasing']
chongtu_type2=['Decreasing','No effect']
chongtu_type3=['Increasing','No effect']
chongtu_type4=['Disrupting','No effect']
chongtu_type5=['Disrupting','Increasing']
index1=[]
for i in range(len(a1)):
	for j in range(len(a1[i][1])):
		if list(a1[i][1]['type'])[0] != list(a1[i][1]['type'])[1]:
			if list(a1[i][1]['type'])[0] in chongtu_type1 and list(a1[i][1]['type'])[1] in chongtu_type1:
				index1.append(i)
			elif list(a1[i][1]['type'])[0] in chongtu_type2 and list(a1[i][1]['type'])[1] in chongtu_type2:
				index1.append(i)
			elif list(a1[i][1]['type'])[0] in chongtu_type3 and list(a1[i][1]['type'])[1] in chongtu_type3:
				index1.append(i)
			elif list(a1[i][1]['type'])[0] in chongtu_type4 and list(a1[i][1]['type'])[1] in chongtu_type4:
				index1.append(i)
			elif list(a1[i][1]['type'])[0] in chongtu_type5 and list(a1[i][1]['type'])[1] in chongtu_type5:
				index1.append(i)

data_chongtu=pd.DataFrame()
for p in index1:
	a2=a1[p][1]
	data_chongtu=pd.concat([data_chongtu,a2])
data_chongtu1=data_chongtu.drop_duplicates(subset=['Interaction AC','Affected protein AC','Mutations','type'],keep='first')
data_chongtu1.to_csv('data_chongtu.tsv',index=False,sep='\t')

index_chongtu=data_chongtu1.index
data2=a.drop(index_chongtu)
data1=pd.read_csv('Reviewed_last1.tsv',header=0,sep='\t')
data3=pd.concat([data1,data2])
data3.to_csv('Reviewed_last_result.tsv',index=False,sep='\t')



'''
shuju_3=pd.read_csv('Reviewed_last_result.tsv',header=0,sep='\t')
print(len(shuju_3))
print(len(shuju_3.drop_duplicates(subset=['Affected protein AC','Mutations'],keep='first')))
print(len(shuju_3.drop_duplicates(subset=['Affected protein AC'],keep='first')))
print(len(shuju_3.drop_duplicates(subset=['Interaction AC'],keep='first')))
print(len(shuju_3.drop_duplicates(subset=['Partners_1'],keep='first')))
'''



