# -*- coding: utf-8 -*-
# @Author: Liu Yang
# @Date:   2023-03-09 20:42:29
# @Last Modified by:   Liu Yang
# @Last Modified time: 2023-03-14 22:06:37

#########################################################################################

import os
os.chdir('D:/abishe/Fifth')
import pandas as pd

#获取单核苷酸错义突变
#问题1：check一下，是否有silence突变
data=pd.read_csv('SNV_mutations.tsv',header=0,sep='\t')
index1=[]
for i in range(len(data)):
	if data['Original sequence'][i]==data['Resulting sequence'][i]:
		index1.append(i)#无

Aa=['G','A','L','I','V','P','F','M','W','S','Q','T','C','N','Y','D','E','K','R','H']
index2=[]
for j in range(len(data)):
	if data['Original sequence'][j] not in Aa or data['Resulting sequence'][j] not in Aa:
		index2.append(j)
data1=data.drop(index2)
data1.to_csv('SNV.tsv',index=False,sep='\t')

'''
print(len(shuju_3))
print(len(shuju_3.drop_duplicates(subset=['Affected protein AC','Feature range(s)','Original sequence','Resulting sequence'],keep='first')))
print(len(shuju_3.drop_duplicates(subset=['Affected protein AC'],keep='first')))
print(len(shuju_3.drop_duplicates(subset=['Interaction AC'],keep='first')))
print(len(shuju_3.drop_duplicates(subset=['Partners_1'],keep='first')))
'''
#只保留Affected protein AC是uniprotkb的
shuju=pd.read_csv('SNV.tsv',header=0,sep='\t')
quchu_type=[]
for i in range(len(shuju)):
    x=str(shuju['Affected protein AC'][i])
    if x[:9]!='uniprotkb':
        quchu_type.append(i)
shuju1=shuju.drop(quchu_type)
shuju_1=shuju.loc[quchu_type]
shuju1.to_csv('shuju.tsv',index=False,sep='\t')


#


####################################################################################################################
#数据分成两类，map突变与相互作用
#将相互作用参与者提取出来
data=pd.read_csv('shuju.tsv',header=0,sep='\t')
data.index=range(len(data))
all_ID=[]
for p in range(len(data)):
	x=data['Interaction participants'][p].split('|')
	ID=[]
	partner_dat=pd.DataFrame()
	for q in range(len(x)):
		y=x[q][:x[q].index('(')]
		ID.append(y)
	all_ID.append(ID)
part_dic={'Partners':all_ID}
part_dat=pd.DataFrame(part_dic)
data_ID=pd.concat([data,part_dat],axis=1)
data_ID.to_csv('data_ID.tsv',index=False,sep='\t')

###################################################################################################################################################
#验证每个Interaction ID对应的Interaction participants一样
data_ID=pd.read_csv('data_ID.tsv',header=0,sep='\t')
AC_ID=data_ID['Interaction AC']
partners=data_ID['Partners']
x=len(list(set(list(AC_ID))))
z=len(data_ID.drop_duplicates(subset=['Interaction AC','Partners'],keep='first'))


data=pd.read_csv('data_ID.tsv',header=0,sep='\t')
a_index=[]
b_index=[]
for i in range(len(data)):
	partner_a=set(eval(data['Partners'][i]))
	if len(partner_a)<3:
		a_index.append(i)
	else:
		b_index.append(i)
two_data=data.loc[a_index]
duo_data=data.loc[b_index]
two_data.to_csv('two_data.tsv',index=False,sep='\t')#59081
duo_data.to_csv('duo_data.tsv',index=False,sep='\t')#4313




#从intact1_jieguo种提取突变中的相互作用
intact1_jieguo=pd.read_csv('intact_3.txt',header=0,sep='\t')
data1=pd.read_csv('two_data.tsv',header=0,sep='\t')

intact1_AC=intact1_jieguo['Interaction AC']
data1_AC=data1['Interaction AC']
data1_result=data1[data1['Interaction AC'].isin(intact1_AC)]
intact1_last=intact1_jieguo[intact1_jieguo['Interaction AC'].isin(data1_AC)]

'''
intact1=list(intact1_AC)
index1=[]
for i in range(len(data1)):
	if data1['Interaction AC'][i] not in intact1:
		index1.append(i)
'''


intact1_last.index=range(len(intact1_last))
#将相互作用Interactor A和Interactor B合在一个列表Partners_1中
Partners=[]
for i in range(len(intact1_last)):
	a=[intact1_last['#ID(s) interactor A'][i],intact1_last['ID(s) interactor B'][i]]
	Partners.append(a)
Partners_data=pd.DataFrame({'Partners_1':Partners})
intact_data=pd.concat([intact1_last,Partners_data],axis=1)
intact_quchong=intact_data.drop_duplicates(subset=['Interaction AC'],keep='first')

data1_last=pd.merge(data1_result,intact_quchong,how='inner',on='Interaction AC')
data1_last.index=range(len(data1_last))
data1_last1=data1_last.drop('#ID(s) interactor A',axis=1)
data1_last2=data1_last1.drop('ID(s) interactor B',axis=1)
data1_last2.to_csv('data_last.tsv',index=False,sep='\t')


##########################################################################################################
#第一次筛选

#删除不是uniprotkb ID的，以及Affected protein AC不在Partners_1里的
data_last=pd.read_csv('data_last.tsv',header=0,sep='\t')
quchu_index_1=[]
quchu_index_2=[]
for i in range(len(data_last)):
	a=data_last['Affected protein AC'][i]
	b=eval(data_last['Partners_1'][i])
	for j in range(len(b)):
		if a not in b:
			quchu_index_1.append(i)
		elif b[j][:9] !='uniprotkb':
			quchu_index_2.append(i)
data_quchuindex1=data_last.loc[quchu_index_1]
data_quchuindex2=data_last.loc[quchu_index_2]

index1=list(set(quchu_index_1+quchu_index_2))
data_1=data_last.drop(index1)
data_quchuindex1.to_csv('data_quchuindex1.tsv',index=False,sep='\t')
data_quchuindex2.to_csv('data_quchuindex2.tsv',index=False,sep='\t')
data_1.to_csv('data_1.tsv',index=False,sep='\t')









#确定互作蛋白是否被uniprot数据库移除
data=pd.read_csv('data_1.tsv',header=0,sep='\t')
pro=[]
for i in range(len(data)):
	for j in range(len(eval(data['Partners_1'][i]))):
		a=eval(data['Partners_1'][i])[j]
		pro.append(a)
pro_1=list(set(pro))
pro_dat=pd.DataFrame({'Protein ID':pro_1})
for p in range(len(pro_dat)):
	temp=pro_dat['Protein ID'][p][10:]
	pro_dat['Protein ID'][p]=temp
pro_dat.to_csv('pro_uniprot.txt',index=False,sep='\t')#将涉及的全部蛋白质导出

#导入uniprot数据库工具ID mapping中搜索
a=pd.read_csv('uniprot-compressed_true_download_true_fields_accession_2Creviewed_2C-2023.03.10-05.27.54.54.tsv',header=0,sep='\t')
a1=a['From']
a2=a['Entry']
index0=[]
#确定uniprotkb id有没有被替换
for i in range(len(a)):
	if '-' in a1[i]:
		temp_a1=a1[i][:a1[i].index('-')]
	else:
		temp_a1=a1[i]

	if temp_a1 != a2[i]:
		index0.append(i)
a_1=a.drop(index0)
a1=list(a_1['From'])

index1=[]
for q in range(len(data)):
	for r in range(len(eval(data['Partners_1'][q]))):
		temp1=eval(data['Partners_1'][q])[r][10:]
		if temp1 not in a1:
			index1.append(q)
			break
data_1=data.loc[index1]
data1=data.drop(index1)
data1.to_csv('baoliu.tsv',index=False,sep='\t')#54649
data_1.to_csv('yichu.tsv',index=False,sep='\t')#1435



########################################################################




#第二次筛选
def baoliu1(shuju):
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

	index_shaixuan1=[]
	for i in range(len(shuju)):
		for j in range(len(eval(shuju['Partners_1'][i]))):
			if eval(shuju['Partners_1'][i])[j][10:] not in P_ID:
				index_shaixuan1.append(i)
	Unreviewed_result=shuju.loc[index_shaixuan1]
	result=shuju.drop(index_shaixuan1)
	result.to_csv('Reviewed_result.tsv',index=False,sep='\t')
	Unreviewed_result.to_csv('Unreviewed_result.tsv',index=False,sep='\t')


a=pd.read_csv('Reviewed_result.tsv',header=0,sep='\t')#42586
b=pd.read_csv('Unreviewed_result.tsv',header=0,sep='\t')#14101



#第三次筛选
def quchu_type(shuju):
	Decreasing=['mutation decreasing(MI:0119)','mutation decreasing strength(MI:1133)','mutation decreasing rate(MI:1130)']
	Increasing=['mutation increasing(MI:0382)','mutation increasing strength(MI:1132)','mutation increasing rate(MI:1131)']
	No_effect=['mutation with no effect(MI:2226)']
	Disrupting=['mutation disrupting(MI:0573)','mutation disrupting strength(MI:1128)','mutation disrupting rate(MI:1129)']
	result=pd.DataFrame()
	qu_index=[]
	for i in range(len(shuju)):
		if shuju['Feature type'][i] in Decreasing:
			type=['Decreasing']
			data_type=shuju[i:i+1]
			data_type['type']=type
			result=pd.concat([result,data_type])
		elif shuju['Feature type'][i] in Increasing:
			type=['Increasing']
			data_type=shuju[i:i+1]
			data_type['type']=type
			result=pd.concat([result,data_type])
		elif shuju['Feature type'][i] in No_effect:
			type=['No effect']
			data_type=shuju[i:i+1]
			data_type['type']=type
			result=pd.concat([result,data_type])
		elif shuju['Feature type'][i] in Disrupting:
			type=['Disrupting']
			data_type=shuju[i:i+1]
			data_type['type']=type
			result=pd.concat([result,data_type])
		else:
			qu_index.append(i)
	return qu_index
	result.to_csv('Reviewed_shaixuan3.tsv',index=False,sep='\t')

#Reviewed蛋白质数量为39451
#Unreviewed蛋白质数量为12813


##################################################################################################################
#查看突变影响蛋白质，有没有isoform，需匹配到canonical上

a=pd.read_csv('Reviewed_shaixuan3.tsv',header=0,sep='\t')#39451
b=pd.read_csv('Unreviewed_shaixuan3.tsv',header=0,sep='\t')#12813

def pipei(shuju):
	f=open('uniprot_sprot_varsplic.fasta')
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
	isoform_dic={'isoform_ID':P_ID,'isoform_seq':P_seq}
	all_isoform=pd.DataFrame(isoform_dic)
	f.close()

	iso_index=[]
	P_ID=list(set(P_ID))
	for p in range(len(shuju)):
		iso=shuju['Affected protein AC'][p]
		if iso in P_ID:
			iso_index.append(p)
	iso_dat=pd.DataFrame({'iso_index':iso_index})
	iso_dat.to_csv('iso_dat.tsv',index=False,sep='\t')
	data_iso=shuju.loc[iso_index]
	data_iso.to_csv('Unreviewed_iso.tsv',index=False,sep='\t')

#Reviewed都没有isoform
#Unreviewed有isoform，



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


'''
data1=pd.read_csv('Unreviewed_shaixuan3.tsv',header=0,sep='\t')
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
data2=pd.concat([d,c,mut_data,e,f,b],axis=1)
data2.to_csv('Unreviewed.tsv',index=False,sep='\t')
'''


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

data=pd.read_csv('Reviewed_seq.tsv',header=0,sep='\t')


index_cuo=[]
for i in range(len(data)):
	a1=int(str(data['Feature range(s)'][i]).split('-')[0])
	a2=data['Original sequence'][i]
	a3=data['Affected protein AC'][i]
	a4=eval(data['Partners_1'][i])
	a_index=a4.index(a3)
	if a_index==0:
		temp1=data['A sequence'][i]
	else:
		temp1=data['B sequence'][i]

	if temp1[a1-1] !=a2:
		index_cuo.append(i)





















###########################
#Unreviewed

###########################################################
#zc数据
a=pd.read_csv('Unreviewed_zc_data.tsv',header=0,sep='\t')
Un_pro=[]
for i in range(len(a)):
	a1=eval(a['Partners_1'][i])[0][10:]
	a2=eval(a['Partners_1'][i])[1][10:]
	Un_pro.append(a1)
	Un_pro.append(a2)
Un_pro1=list(set(Un_pro))
Un_pro_dat=pd.DataFrame({'protein ID':Un_pro1})
Un_pro_dat.to_csv('Unreviewed_zc_id.txt',index=False,sep='\t')

#添加序列
b=pd.read_csv('Unreviewed_zc_data.tsv',header=0,sep='\t')
#读取fasta格式的uniprot现有且审查蛋白数据,并进行处理
f=open('uniprot-compressed_true_download_true_format_fasta-2023.03.10-09.05.18.25.fasta')
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
	if '-' in temp1:
		temp1_1=eval(a['Partners_1'][i])[0][10:eval(a['Partners_1'][i])[0].index('-')]
	else:
		temp1_1=eval(a['Partners_1'][i])[0][10:]
	if '-' in temp2:
		temp2_1=eval(a['Partners_1'][i])[1][10:eval(a['Partners_1'][i])[1].index('-')]
	else:
		temp2_1=eval(a['Partners_1'][i])[1][10:]

	seq1=list(Uniprot[Uniprot['Protein_ID']==temp1_1]['Protein_seq'])[0]
	seq2=list(Uniprot[Uniprot['Protein_ID']==temp2_1]['Protein_seq'])[0]
	seqs1.append(seq1)
	seqs2.append(seq2)
seq1_dat=pd.DataFrame({'A sequence':seqs1})
seq2_dat=pd.DataFrame({'B sequence':seqs2})
data=pd.concat([a,seq1_dat,seq2_dat],axis=1)
data.to_csv('Unevriewed_zc_seq.tsv',index=False,sep='\t')




