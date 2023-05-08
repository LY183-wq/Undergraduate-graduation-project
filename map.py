

###################################################################################################################################################
###################################################################################################################################################
#方案1：
#获取相互作用中蛋白b的结构信息，最后与a的pdb信息取交集
#获取a的信息
import os
os.chdir('C:\\Users\\jiang\\Desktop')
import pandas as pd

import requests
import json
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import PDB

a=pd.read_csv('Reviewed1.tsv',header=0,sep='\t')
a1=list(a['Affected protein AC'])
IDs=list(set(a1))

def get_message(x):
	for k in list(x.keys()):
		x[k]=[x[k]]
	x_dat=pd.DataFrame(x)
	return x_dat

no_pdb=[]
yes_pdb=[]
result1=pd.DataFrame()
for i in IDs:
	url='https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures/'+i
	file_i=requests.get(url,verify=False)
	r=json.loads(file_i.content.decode())
	if r=={}:
		no_pdb.append(i)
	else:
		uni_id=list(r.keys())[0]
		pdb_message=list(r.values())[0]
		pdb_nums=len(pdb_message)
		uni_ids=[uni_id]*pdb_nums
		uni_dat=pd.DataFrame({'Uniprot_id':uni_ids})
		j_dats=pd.DataFrame()
		for j in range(pdb_nums):
			j_dat=get_message(pdb_message[j])
			j_dats=pd.concat([j_dats,j_dat])
		j_dats.index=range(len(j_dats))
		result_i=pd.concat([uni_dat,j_dats],axis=1)
		result1=pd.concat([result1,result_i])

result1.to_csv('result1.tsv',index=False,sep='\t')

#获取b的信息
import os
os.chdir('C:\\Users\\jiang\\Desktop')
import pandas as pd

import requests
import json
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import PDB

a=pd.read_csv('Reviewed1.tsv',header=0,sep='\t')
a1=list(a['b id'])
IDs=list(set(a1))

def get_message(x):
	for k in list(x.keys()):
		x[k]=[x[k]]
	x_dat=pd.DataFrame(x)
	return x_dat

no_pdb=[]
yes_pdb=[]
result1=pd.DataFrame()
for i in IDs:
	url='https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures/'+i
	file_i=requests.get(url,verify=False)
	r=json.loads(file_i.content.decode())
	if r=={}:
		no_pdb.append(i)
	else:
		uni_id=list(r.keys())[0]
		pdb_message=list(r.values())[0]
		pdb_nums=len(pdb_message)
		uni_ids=[uni_id]*pdb_nums
		uni_dat=pd.DataFrame({'Uniprot_id':uni_ids})
		j_dats=pd.DataFrame()
		for j in range(pdb_nums):
			j_dat=get_message(pdb_message[j])
			j_dats=pd.concat([j_dats,j_dat])
		j_dats.index=range(len(j_dats))
		result_i=pd.concat([uni_dat,j_dats],axis=1)
		result1=pd.concat([result1,result_i])

result1.to_csv('result_b.tsv',index=False,sep='\t')


#取pdb交集
import os
os.chdir('C:\\Users\\jiang\\Desktop\\get_pdb\\Reviewed_method1')
import pandas as pd
a=pd.read_csv('result2.tsv',header=0,sep='\t')
b=pd.read_csv('result2_b.tsv',header=0,sep='\t')

a1=list(a.groupby('Uniprot_id'))
b1=list(b.groupby('Uniprot_id'))

pro_a_dic=dict()
for i in range(len(a1)):
	pro_id=a1[i][0]
	pdb_group_i=list(a1[i][1]['pdb_id'])
	pro_a_dic[pro_id]=pdb_group_i

pro_b_dic=dict()
for j in range(len(b1)):
	pro_id_b=b1[j][0]
	pdb_group_j=list(b1[j][1]['pdb_id'])
	pro_b_dic[pro_id_b]=pdb_group_j

c=pd.read_csv('Reviewed1.tsv',header=0,sep='\t')
c_a=list(set(list(c['Affected protein AC'])))
no_pdb_a=[x for x in c_a if x not in list(a['Uniprot_id'])]
c_b=list(set(list(c['b id'])))
no_pdb_b=[y for y in c_b if y not in list(b['Uniprot_id'])]

pdb_maps=[]
for m in range(len(c)):
	temp_a=c['Affected protein AC'][m]
	temp_b=c['b id'][m]
	if temp_a in no_pdb_a or temp_b in no_pdb_b:
		pdb_maps.append('0')
	else:
		pdb_a=pro_a_dic[temp_a]
		pdb_b=pro_b_dic[temp_b]
		pdb_both=[x for x in pdb_a if x in pdb_b]
		pdb_maps.append(pdb_both)
pdb_dat=pd.DataFrame({'pdb_maps':pdb_maps})
c1=pd.concat([c,pdb_dat],axis=1)
c1_false1=c1[c1['pdb_maps']=='0']
c1_f1_index=list(c1_false1.index)
c2=c1.drop(c1_f1_index)
c2.index=range(len(c2))
c2_f2_index=[]
for p in range(len(c2)):
	if len(c2['pdb_maps'][p])==0:
		c2_f2_index.append(p)
c3=c2.drop(c2_f2_index)
c3.index=range(len(c3))
#去除重复pdb
all_pdb=[]
for q in range(len(c3)):
	temp_pdb=list(set(c3['pdb_maps'][q]))
	all_pdb=all_pdb+temp_pdb
all_pdb1=list(set(all_pdb))


import os
import requests
import json
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import PDB
def get_message(x):
	pdb=list(x.keys())[0]
	pdb_dat=pd.DataFrame({'pdb_id':[pdb]})
	m=list(list(x.values())[0].values())[0]
	dat1=pd.DataFrame()
	index2=0
	for i in list(m.keys()):
		b1=str(list(m.keys()).index(i))
		m_i=m[i]
		map1=m_i['mappings']
		if len(map1)>1:
			for p in range(len(map1)):
				p1=str(p)
				index1=str(index2)
				m_i_dic={'Uniprot_id_'+index1:i,'identifier_'+index1:[m_i['identifier']],'entity_id_'+index1:[m_i['mappings'][p]['entity_id']],'chain_id_'+index1:[m_i['mappings'][p]['chain_id']],'unp_end_'+index1:[m_i['mappings'][p]['unp_end']],
				'unp_start_'+index1:[m_i['mappings'][p]['unp_start']],'pdb_start_'+index1:[m_i['mappings'][p]['pdb_start']],'pdb_end_'+index1:[m_i['mappings'][p]['pdb_end']],
				'is_canonical_'+index1:[m_i['mappings'][p]['is_canonical']]}
				m_i_dat=pd.DataFrame(m_i_dic)
				dat1=pd.concat([dat1,m_i_dat],axis=1)
				index2=index2+1
		else:
			index1=str(index2)
			m_i_dic={'Uniprot_id_'+index1:i,'identifier_'+index1:[m_i['identifier']],'entity_id_'+index1:[m_i['mappings'][0]['entity_id']],'chain_id_'+index1:[m_i['mappings'][0]['chain_id']],'unp_end_'+index1:[m_i['mappings'][0]['unp_end']],
			'unp_start_'+index1:[m_i['mappings'][0]['unp_start']],'pdb_start_'+index1:[m_i['mappings'][0]['pdb_start']],'pdb_end_'+index1:[m_i['mappings'][0]['pdb_end']],
			'is_canonical_'+index1:[m_i['mappings'][0]['is_canonical']]}
			m_i_dat=pd.DataFrame(m_i_dic)
			dat1=pd.concat([dat1,m_i_dat],axis=1)
			index2=index2+1
	dat2=pd.concat([pdb_dat,dat1],axis=1)
	return dat2

no_pro=[]
duo=[]
data1=pd.DataFrame()
for j in all_pdb1:
	url='https://www.ebi.ac.uk/pdbe/api/mappings/all_isoforms/'+j
	file_j=requests.get(url,verify=False)
	r=json.loads(file_j.content.decode())
	if r=={}:
		no_pro.append(j)
	else:
		b=get_message(r)
		data1=pd.concat([data1,b])
#data1.to_csv('Reviewed_map_415_1.tsv',index=False,sep='\t',na_rep='_NAN_')
#xxx=pd.read_csv('Reviewed_map_415_1.tsv',header=0,sep='\t',keep_default_na=False, na_values=['_NAN_'])
#因为有链的id为'NA'，pandas会自动将其读取成NaN
#因此用na_rep参数将真正的NaN先保存为'_NAN_'，读取的时候再退回来


###################################################################################################################################################
###################################################################################################################################################
#Unreviewed
#按上面流程，运行best structures和SIFTS mapping两个接口，获取结果文件Unreviewed_map.tsv
import os
os.chdir('/data/liuyang/bishe_data/Reviewed_method1')
import pandas as pd

#取交集
a=pd.read_csv('Unreviewed_best_A.tsv',header=0,sep='\t')
b=pd.read_csv('Unreviewed_best_B.tsv',header=0,sep='\t')
a1=list(a.groupby('Uniprot_id'))
b1=list(b.groupby('Uniprot_id'))
pro_a_dic=dict()
for i in range(len(a1)):
	pro_id=a1[i][0]
	pdb_group_i=list(a1[i][1]['pdb_id'])
	pro_a_dic[pro_id]=pdb_group_i
pro_b_dic=dict()
for j in range(len(b1)):
	pro_id_b=b1[j][0]
	pdb_group_j=list(b1[j][1]['pdb_id'])
	pro_b_dic[pro_id_b]=pdb_group_j
c=pd.read_csv('Unreviewed_right.tsv',header=0,sep='\t')
c_a=list(set(list(c['Affected protein AC'])))
no_pdb_a=[x for x in c_a if x not in list(a['Uniprot_id'])]
c_b=list(set(list(c['b id'])))
no_pdb_b=[y for y in c_b if y not in list(b['Uniprot_id'])]
pdb_maps=[]
for m in range(len(c)):
	temp_a=c['Affected protein AC'][m]
	temp_b=c['b id'][m]
	if temp_a in no_pdb_a or temp_b in no_pdb_b:
		pdb_maps.append('0')
	else:
		pdb_a=pro_a_dic[temp_a]
		pdb_b=pro_b_dic[temp_b]
		pdb_both=[x for x in pdb_a if x in pdb_b]
		pdb_maps.append(pdb_both)
pdb_dat=pd.DataFrame({'pdb_maps':pdb_maps})
c1=pd.concat([c,pdb_dat],axis=1)
c1_false1=c1[c1['pdb_maps']=='0']
c1_f1_index=list(c1_false1.index)
c2=c1.drop(c1_f1_index)
c2.index=range(len(c2))
c2_f2_index=[]
for p in range(len(c2)):
	if len(c2['pdb_maps'][p])==0:
		c2_f2_index.append(p)
c3=c2.drop(c2_f2_index)
c3.index=range(len(c3))    #c3即为有交集的数据，即存在复合物pdb的突变数据
#去除重复pdb
pdb_maps1=[]
for k in range(len(c3)):
	temp_c=list(set(c3['pdb_maps'][k]))
	pdb_maps1.append(temp_c)
pdb_dat1=pd.DataFrame({'pdb_id':pdb_maps1})
c4=pd.concat([c3,pdb_dat1],axis=1)
#将pdb数据提取出来
c5=c4.explode('pdb_id')

m=pd.read_csv('Unreviewed_map.tsv',header=0,sep='\t')
#pdb数据整合到突变数据上
c6=pd.merge(c5,m,how='left',on='pdb_id')

#删除不是复合物的数据，输入pandas数据（通过Uniprot_id_1有无蛋白质id来进行判断）
def quchu_danti(x):
	index_qu=[]
	for i in range(len(x)):
		if pd.isnull(x['Uniprot_id_1'][i]):
			index_qu.append(i)
	x1=x.drop(index_qu)
	x1.index=range(len(x1))
	return x1

#########################################################################
#判断每个复合物上：Affected protein链的数量，b protein链的数量，总共链的数量
def panduan_A(x):
	column1=list(x.columns.values)
	all_num=int(len([y for y in column1 if y[-1].isdigit()])/9)    #和单体蛋白质相关的信息列，最后一个字符都为数字
	A_chains=[]
	for i in range(len(x)):
		A_chain_i=[]
		for j in range(all_num):
			uniprot_j=x['Uniprot_id_'+str(j)][i]
			if uniprot_j==x['Affected protein AC'][i]:
				chain_j=x['chain_id_'+str(j)][i]
				A_chain_i.append(chain_j)
				A_chain_i_1=list(set(A_chain_i))
		A_chains.append(A_chain_i_1)
	A_chains_dat=pd.DataFrame({'A chains':A_chains})
	x1=pd.concat([x,A_chains_dat],axis=1)
	return x1

#判断b蛋白质链的数量
def panduan_B(x):
	column1=list(x.columns.values)
	all_num=int(len([y for y in column1 if y[-1].isdigit()])/9)
	B_chains=[]
	for i in range(len(x)):
		B_chain_i=[]
		for j in range(all_num):
			uniprot_j=x['Uniprot_id_'+str(j)][i]
			if uniprot_j==x['b id'][i]:
				chain_j=x['chain_id_'+str(j)][i]
				B_chain_i.append(chain_j)
				B_chain_i_1=list(set(B_chain_i))
		B_chains.append(B_chain_i_1)
	B_chains_dat=pd.DataFrame({'B chains':B_chains})
	x1=pd.concat([x,B_chains_dat],axis=1)
	return x1

#判断总的蛋白链数
def panduan_all(x):
	column1=list(x.columns.values)
	all_num=int(len([y for y in column1 if y[-1].isdigit()])/9)
	all_chains=[]
	for i in range(len(x)):
		all_chain_i=[]
		for j in range(all_num):
			uniprot_j=x['Uniprot_id_'+str(j)][i]
			chain_j=x['chain_id_'+str(j)][i]
			if pd.isnull(uniprot_j):
				break
			else:
				all_chain_i.append(chain_j)
		all_chain_i_1=list(set(all_chain_i))
		all_chains.append(all_chain_i_1)
	all_chains_dat=pd.DataFrame({'all chains':all_chains})
	x1=pd.concat([x,all_chains_dat],axis=1)
	return x1

c7=quchu_danti(c6)
c8=panduan_A(c7)
c9=panduan_B(c8)
c10=panduan_all(c9)


		

#添加A和B的覆盖率、pdb分辨率，实验方法
a1=pd.read_csv('Unreviewed_best_A.tsv',header=0,sep='\t',keep_default_na=False, na_values=['_NAN_'])
b1=pd.read_csv('Unreviewed_best_B.tsv',header=0,sep='\t',keep_default_na=False, na_values=['_NAN_'])
a2=a1.drop_duplicates(subset=['Uniprot_id','pdb_id','chain_id'],keep='first')
a2.index=range(len(a2))
a=a2.drop(['end','start','unp_end','unp_start','tax_id'],axis=1)
b2=b1.drop_duplicates(subset=['Uniprot_id','pdb_id','chain_id'],keep='first')
b2.index=range(len(b2))
b=b2.drop(['end','start','unp_end','unp_start','tax_id'],axis=1)
chains_a=list(c10['A chains'])
chains_b=list(c10['B chains'])
A_cover_group=[]
B_cover_group=[]
for i in range(len(c10)):
	temp_pro=c10['Affected protein AC'][i]
	temp_b_pro=c10['b id'][i]
	temp_pdb=c10['pdb_id'][i]
	A_coverage=[]
	B_coverage=[]
	for j in chains_a[i]:
		a_temp=a.loc[(a['Uniprot_id']==temp_pro) & (a['pdb_id']==temp_pdb) & (a['chain_id']==j)]
		a_temp.index=range(len(a_temp))
		a_coverage_j=a_temp['coverage'][0]
		A_coverage.append(a_coverage_j)
	for k in chains_b[i]:
		b_temp=b.loc[(b['Uniprot_id']==temp_b_pro) & (b['pdb_id']==temp_pdb) & (b['chain_id']==k)]
		b_temp.index=range(len(b_temp))
		b_coverage_k=b_temp['coverage'][0]
		B_coverage.append(b_coverage_k)
	A_cover_group.append(A_coverage)
	B_cover_group.append(B_coverage)
A_cover_dat=pd.DataFrame({'A coverage':A_cover_group})
B_cover_dat=pd.DataFrame({'B coverage':B_cover_group})
c11=pd.concat([c10,A_cover_dat,B_cover_dat],axis=1)

a_pdb=a.drop_duplicates(subset=['pdb_id'],keep='first')
a_pdb.index=range(len(a_pdb))
a_resolution=pd.concat([a_pdb['pdb_id'],a_pdb['resolution'],a_pdb['experimental_method']],axis=1)
c12=pd.merge(c11,a_resolution,how='left',on='pdb_id')
c12.to_csv('Unreviewed_map_result.tsv',index=False,sep='\t')

#Reviewed_map_418.tsv和Unreviewed_map_417.tsv

#要不要对结构进行筛选
'''
#结构筛选原则，以Reviewed数据为例
#对于只有一个结构的数据，保留该数据
#对于有多个结构的数据，通过覆盖率和分辨率两个指标来筛选出最佳结构
import os
os.chdir('/data/liuyang/bishe_data/Reviewed_method1')
import pandas as pd
a=pd.read_csv('Reviewd_map_result.tsv',header=0,sep='\t')
#去除重复数据
a1=a.drop_duplicates(subset=['Affected protein AC','Mutations','b id','pdb_id'],keep='first')

'''




#首先将一定能够用的数据提取出来,，以Reviewed_map_418.tsv为例
import os
os.chdir('/data/liuyang/bishe_data/Reviewed_method1/3')
import pandas as pd
a1=pd.read_csv('Reviewed_map_418.tsv',header=0,sep='\t')
a=a1.drop_duplicates(subset=['Affected protein AC','Mutations','b id','pdb_id'],keep='first')
a.index=range(len(a))

def panduan_all_A(x):
	column1=list(x.columns.values)
	all_num=int(len([y for y in column1 if y[-1].isdigit()])/9)    #和单体蛋白质相关的信息列，最后一个字符都为数字
	A_chains=[]
	for i in range(len(x)):
		A_chain_i=[]
		for j in range(all_num):
			uniprot_j=x['Uniprot_id_'+str(j)][i]
			if uniprot_j==x['Affected protein AC'][i]:
				chain_j=x['chain_id_'+str(j)][i]
				A_chain_i.append(chain_j)
				A_chain_i_1=A_chain_i
		A_chains.append(A_chain_i_1)
	A_chains_dat=pd.DataFrame({'A all chains':A_chains})
	x1=pd.concat([x,A_chains_dat],axis=1)
	return x1
a_all=panduan_all_A(a)

#判断突变位点在不在该复合物对应链上，在的话，直接作为突变链
def mut_in_chains(x):
	import re
	x.index=range(len(x))
	column1=list(x.columns.values)
	all_num=int(len([y for y in column1 if y[-1].isdigit()])/9)
	mut_chains=[]
	all_a_chains=x['A all chains'][0]
	loc1=float(re.findall(r'\d+',x['Mutations'][0])[0])
	id=[]
	for j in all_a_chains:
		for p in range(all_num):
			if x['chain_id_'+str(p)][0]==j and loc1 >= float(x['unp_start_'+str(p)][0]) and loc1 <= float(x['unp_end_'+str(p)][0]):
				mut_chains.append(j)
				id.append(p)
				break
		else:
			continue
		break
	mut_dat=pd.DataFrame({'mut_chain':mut_chains})
	id_dat=pd.DataFrame({'id':id})
	x1=pd.concat([x,mut_dat],axis=1)
	return x1

data1=pd.DataFrame()
for i in range(len(a_all)):
	temp_a=a_all[i:i+1]
	temp_data=mut_in_chains(temp_a)
	data1=pd.concat([data1,temp_data])


#验证结果，判断data2和index是不是一样长度
import math
data1.index=range(len(data1))
index1=[]
for i in range(len(data1)):
	if math.isnan(data1['id'][i]):
		index1.append(i)
data2=data1.drop(index1)
data2.index=range(len(data2))
index2=[]
for i in range(len(data2)):
	loc1=float(re.findall(r'\d+',data2['Mutations'][i])[0])
	if float(data2['unp_start_'+str(int(data2['id'][i]))]) <= loc1 and float(data2['unp_end_'+str(int(data2['id'][i]))]) >= loc1:
		index2.append(i)


a_all=data2
data3=pd.concat([a_all['Affected protein AC'],a_all['b id'],a_all['type'],a_all['Mutations'],
	a_all['pdb_id'],a_all['resolution'],a_all['A chains'],a_all['B chains'],a_all['all chains'],
	a_all['A coverage'],a_all['B coverage'],a_all['A all chains'],a_all['mut_chain']],axis=1)
data3.to_csv('Reviewed_before_select.tsv',index=False,sep='\t')

#Reivewed中共有48532

#将两套数据移入路径PremPS2下
#还要去掉partner1和partner2为同一条链的
import os
os.chdir('/data/liuyang/bishe_data/PremPS2')
import pandas as pd
a=pd.read_csv('Reviewed_before_select.tsv',header=0,sep='\t')
index1=[]
for i in range(len(a)):
    if eval(a['A chains'][i])==eval(a['B chains'][i]) and len(eval(a['A chains'][i]))==1:
        index1.append(i)
a1=a.drop(index1)
a1.to_csv('Reviewed_before_select_1.tsv',index=False,sep='\t')


#每个突变数据，选出一个最佳结构，先通过A和B覆盖率，再通过分辨率来进行筛选
import os
os.chdir('/data/liuyang/bishe_data/PremPS2')   #
import pandas as pd
a=pd.read_csv('Unreviewed_before_select_1.tsv',header=0,sep='\t')

#判断b链的覆盖率是否都一样##########
index_a=[]
for i in range(len(a)):
	if len(set(eval(a['B coverage'][i])))>1:
		index_a.append(i)
#经过判断Unreviewed中，b链覆盖率都一样
#Reviewed中，有313个b链覆盖率不一样

#函数best_pdb，在每个突变数据中，按照A、B覆盖率和分辨率的顺序筛选出一个最佳结构
def best_pdb(x):
	x1=list(x.groupby(['Affected protein AC','b id','type','Mutations']))
	data_last=pd.DataFrame()
	for i in range(len(x1)):
		x1_i=x1[i][1]
		x1_i.index=range(len(x1_i))
		a_list=[]
		b_list=[]
		for j in range(len(x1_i)):
			chain_a_j=x1_i['mut_chain'][j]
			index_a_j=eval(x1_i['A chains'][j]).index(chain_a_j)
			cover_a_j=eval(x1_i['A coverage'][j])[index_a_j]
			cover_b_j=max(eval(x1_i['B coverage'][j]))
			a_list.append(cover_a_j)
			b_list.append(cover_b_j)
		a_dat=pd.DataFrame({'A cover':a_list})
		b_dat=pd.DataFrame({'B cover':b_list})
		x1_all=pd.concat([x1_i,a_dat,b_dat],axis=1)
		x2=x1_all.sort_values(['A cover','B cover','resolution'],ascending=[False,False,True])
		x2.index=range(len(x2))
		pdb_i=x2[0:1]
		data_last=pd.concat([data_last,pdb_i])
	return data_last
result=best_pdb(a)
result.to_csv('Unreviewed_selected.tsv',index=False,sep='\t')
#Unreviewed剩余366个
#Reviewed剩余10971个

#以Unreviewed_selected.tsv为例，将数据移入新路径/data/liuyang/bishe_data/PremPS1
import re
import os
os.chdir('/data/liuyang/bishe_data/PremPS2')
import pandas as pd
a=pd.read_csv('Unreviewed_selected.tsv',header=0,sep='\t')
a1=pd.concat([a['Affected protein AC'],a['pdb_id'],a['A chains'],a['Mutations'],a['mut_chain']],axis=1)
chain_id=[]
locations=[]
for i in range(len(a1)):
	chain_id.append(a1['mut_chain'][i])
	loc1=int(re.findall(r'\d+',a1['Mutations'][i])[0])
	locations.append(loc1)
chain_dat=pd.DataFrame({'chain_id':chain_id})
loc_dat=pd.DataFrame({'unp_sites':locations})
a2=a1.drop('A chains',axis=1)
a3=pd.concat([a2,chain_dat],axis=1)
a4=a3.drop('Mutations',axis=1)
a5=pd.concat([a4,loc_dat],axis=1)
a6=list(a5.groupby(['Affected protein AC','pdb_id','chain_id']))
unp_sites_1=[]
pro=[]
pdb=[]
chains=[]
for j in range(len(a6)):
	mut_sites=list(a6[j][1]['unp_sites'])
	mut_set1=list(set(mut_sites))
	mut_set=sorted(mut_set1)
	result = ','.join(str(x) for x in mut_set)
	unp_sites_1.append(result)
	pro.append(a6[j][0][0])
	pdb.append(a6[j][0][1])
	chains.append(a6[j][0][2])
a7=a5.drop('unp_sites',axis=1)
unp_data_1=pd.DataFrame({'unp_sites':unp_sites_1})
pro_dat=pd.DataFrame({'unp_id':pro})
pdb_dat=pd.DataFrame({'pdb_id':pdb})
chain_dats=pd.DataFrame({'chain_id':chains})
a8=pd.concat([pro_dat,pdb_dat,chain_dats,unp_data_1],axis=1)
a8.to_csv('Unreviewed_before.tsv',index=False,sep='\t')


#学姐的pdb_profiling工具，将uniprot序列位点，转换为pdb位点
#以Reviewed数据为例
import re
import os
os.chdir('/data/liuyang/bishe_data/PremPS2')
import pandas as pd
a=pd.read_csv('Unreviewed_selected.tsv',header=0,sep='\t')
b=pd.read_csv('unreviewed_unp_230420.txt',header=0,sep='\t')
pro_pdb=dict()
for i in range(len(b)):
	temp_pro=b['uni_list'][i].split(',')
	temp_pdb=b['pdb_list'][i].split(',')
	for j in range(len(temp_pro)):
		pro_pdb[temp_pro[j]]=temp_pdb[j]
pdb_sites=[]
for m in range(len(a)):
	loc_m=int(re.findall(r'\d+',a['Mutations'][m])[0])
	temp_mutation=a['Affected protein AC'][m]+'_'+str(loc_m)
	if pro_pdb[temp_mutation] != 'NONE':
		pdb_site=str(pro_pdb[temp_mutation][pro_pdb[temp_mutation].index('_',pro_pdb[temp_mutation].index('_')+1)+1:])
		pdb_sites.append(pdb_site)
	else:
		pdb_sites.append('NONE')
pdb_dat=pd.DataFrame({'pdb_site':pdb_sites})
a1=pd.concat([a,pdb_dat],axis=1)
a1.to_csv('Unreviewed_after.tsv',index=False,sep='\t')

#Unreviewed:
#Reviewed:11095
##########################################################################################################################################################
##处理数据为PremPS运行格式
import os
os.chdir('/data/liuyang/bishe_data/PremPS2')
import pandas as pd

def panduan(x):
	from collections import Counter
	x1 = Counter(x)
	y=dict(x1)
	z=''
	for i in sorted(list(y.keys())):
		for j in range(y[i]):
			z=z+'.'+i+'_'+str(j+1)
	z1=z[1:]
	return z1

a1=pd.read_csv('Unreviewed_after.tsv',header=0,sep='\t')
#去除pdb_site为NONE的数据
a_qu=list(a1[a1['pdb_site']=='NONE'].index)
a=a1.drop(a_qu)
a.index=range(len(a))
pdb_ids=[]
chains_id=[]
MutChains=[]
Mutations=[]
result_ids=[]
pdb_chains=[]
for i in range(len(a)):
	pdb_id=a['pdb_id'][i].upper()+'.pdb'
	Chain=eval(a['A chains'][i])
	Chains=panduan(Chain)
	temp_chain=eval(a['B chains'][i])
	temp_Chains=panduan(temp_chain)
	all_Chain=Chains+'.'+temp_Chains
	MutChain=a['mut_chain'][i]+'_1'
	Mutation=a['Mutations'][i][0]+str(a['pdb_site'][i])+a['Mutations'][i][-1]
	pdb_chain=a['pdb_id'][i].lower()+'_'+a['mut_chain'][i].lower()

	pdb_ids.append(pdb_id)
	chains_id.append(all_Chain)
	MutChains.append(MutChain)
	Mutations.append(Mutation)
	pdb_chains.append(pdb_chain)

pdb_dat=pd.DataFrame({'PDBfile':pdb_ids})
chains_dat=pd.DataFrame({'Chains':chains_id})
MutChains_dat=pd.DataFrame({'MutChain':MutChains})
Mutations_dat=pd.DataFrame({'Mutation_PDB':Mutations})
pdb_chains_dat=pd.DataFrame({'pdb_chains':pdb_chains})
a1=pd.concat([pdb_dat,chains_dat,MutChains_dat,Mutations_dat,pdb_chains_dat],axis=1)

#由于同聚体的原因，一样的链会写两次，因此进行去重
x2_list=[]
for i in range(len(a1)):
	x=a1['Chains'][i].split('.')
	x1=sorted(list(set(x)))
	x2='.'.join(x1)
	x2_list.append(x2)
x2_dat=pd.DataFrame({'Chains':x2_list})
a1_1=a1.drop('Chains',axis=1)
a1_2=pd.concat([a1_1,x2_dat],axis=1)
a1=a1_2

a2=list(a1.groupby(['PDBfile','MutChain']))
isPI=[]
result_data=pd.DataFrame()
for j in range(len(a2)):
	temp_m=['m'+str(j+1)]*len(a2[j][1])
	temp_dat=pd.DataFrame({'isPI':temp_m})
	temp_a=a2[j][1]
	temp_a.index=range(len(temp_a))
	temp_data=pd.concat([temp_a,temp_dat],axis=1)
	result_data=pd.concat([result_data,temp_data])

result_data1=result_data
result_data1.index=range(len(result_data1))
temp_re=[x+1 for x in range(len(result_data1))]
result_dat=pd.DataFrame({'Result_Id':temp_re})
result_data2=result_data1
result_data3=pd.concat([result_data2,result_dat],axis=1)


#Unreviewed:
#Reviewed:
#Unreviewed数据中的m61的符合物6zji没有大型pdb结构
#删除其中没有pdb结构的数据
import requests
import json
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import PDB
import requests
import os
ok_pdb=[]
no_pdb=[]
pdbfiles=list(set(list(result_data3['PDBfile'])))
for i in pdbfiles:
	url='https://files.rcsb.org/download/'+i
	file_i=requests.get(url,verify=False)
	if file_i.status_code == 200:
		ok_pdb.append(i)
	else:
		no_pdb.append(i)

no_data=result_data3[result_data3['PDBfile'].isin(no_pdb)]
no_index=no_data.index
result_data4=result_data3.drop(no_index)
result_data4.to_csv('Unreviewed_v1_single_1.tsv',index=False,sep='\t')



#运行PremPS
from itertools import groupby
from collections import defaultdict,deque,Counter
from pandas import DataFrame

##设置工作路径
workdir = '/data/liuyang/PremPS/'
pathdata = workdir+'datasets/'

def produce_input():
    usedcols = ['PDBfile','Chains','MutChain','Mutation_PDB','Result_Id','isPI']

    f = pd.read_csv(pathdata+'Unreviewed_v1_single_1.tsv', sep='\t')
    f = f.loc[:,usedcols].drop_duplicates(['PDBfile','Chains','MutChain','Mutation_PDB'])
    for _,group in f.groupby(['isPI']):
        jobid = ''.join(list(set(group['isPI'])))
        x=''.join(list(set(group['PDBfile'])))
        pdbfile = x[:x.index('.')].upper()+'.pdb'
#         with open(workdir+jobid+'.sh','w') as fw:
#             fw.write('python premps.py -i %s > %s/%s.log 2>&1\npython premps_cyt.py -i %s > %s/%s.log1 2>&1' % (jobid,jobid,jobid,jobid,jobid,jobid))
        os.system('mkdir %s' % (workdir+jobid))
        os.system('wget https://files.rcsb.org/download/%s -O %s' % (pdbfile,workdir+jobid+'/'+pdbfile))
        group.loc[:,usedcols].to_csv(workdir+jobid+'/'+jobid+'.input',sep='\t',index=0)
    # alljobids,xargs运行需要调用
    with open(workdir+'jobids','w') as fw:
        fw.write('\n'.join(set(f['isPI'])))



#通过pdb结构，判断有几条链
import os
os.chdir('/data/liuyang/PremPS/datasets')
import pandas as pd
f = pd.read_csv('Reviewed_v1_single_1.tsv',header=0, sep='\t')
dict1=dict()
for j in range(len(f)):
    temp_1=f['PDBfile'][j][:f['PDBfile'][j].index('.')].lower()
    temp_2=f['isPI'][j]
    dict1[temp_1]=temp_2
a1=pd.read_csv('/data/liuyang/Reviewed_after.tsv',header=0,sep='\t')
a_qu=list(a1[a1['pdb_site']=='NONE'].index)
a2=a1.drop(a_qu)
a2.index=range(len(a2))
pdb1=list(dict1.keys())
a=a2[a2['pdb_id'].isin(pdb1)]
a.index=range(len(a))
chain_pdb=[]
workdir='/data/liuyang/Reviewed_PremPS/'
for i in range(len(a)):
    m_i=dict1[a['pdb_id'][i]]
    pdb_file=workdir+m_i+'/'+a['pdb_id'][i].upper()+'.pdb'
    if os.path.exists(pdb_file):
        from Bio.PDB import PDBParser
        parser = PDBParser()
        structure = parser.get_structure(a['pdb_id'][i].upper(),pdb_file)
        chain_ids = []
        for model in structure:
            for chain in model:
                chain_id = chain.get_id()
                if chain_id not in chain_ids:
                    chain_ids.append(chain_id)
        chain_pdb.append(chain_ids)
    else:
        chain_pdb.append('NONE')
chain_dat=pd.DataFrame({'pdb_chains':chain_pdb})
ax=pd.concat([a,chain_dat],axis=1)
ax.to_csv('Reviewed_after_PDBchains.tsv',index=False,sep='\t')



#读取PremPS计算结果，最终结果map到Unreviewed_selected.tsv上
import re
import os
os.chdir('/data/liuyang')
import pandas as pd
a=pd.read_csv('Unreviewed_selected.tsv',header=0,sep='\t')
b=pd.read_csv('Unreviewed_after_PDBchains.tsv',header=0,sep='\t')
b1=b[b['pdb_chains']!='NONE']
b2=b1.drop_duplicates(subset=['pdb_id'],keep='first')
b2.index=range(len(b2))
b3=pd.concat([b2['pdb_id'],b2['pdb_chains']],axis=1)
b_pdb=list(set(list(b1['pdb_id'])))
a1=a[a['pdb_id'].isin(b_pdb)]
a1.index=range(len(a1))
a2=pd.merge(a1,b3,how='left',on='pdb_id')

d=pd.read_csv('unreviewed_unp_230422.txt',header=0,sep='\t')
pro_pdb=dict()
for i in range(len(d)):
	temp_pro=d['uni_list'][i].split(',')
	temp_pdb=d['pdb_list'][i].split(',')
	for j in range(len(temp_pro)):
		pro_pdb[temp_pro[j]]=temp_pdb[j]

pdb_sites=[]
for m in range(len(a2)):
	loc_m=int(re.findall(r'\d+',a2['Mutations'][m])[0])
	temp_mutation=a2['Affected protein AC'][m]+'_'+str(loc_m)
	if pro_pdb[temp_mutation] != 'NONE':
		pdb_site=str(pro_pdb[temp_mutation][pro_pdb[temp_mutation].index('_',pro_pdb[temp_mutation].index('_')+1)+1:])
		pdb_sites.append(pdb_site)
	else:
		pdb_sites.append('NONE')
pdb_dat=pd.DataFrame({'pdb_site':pdb_sites})
a_m=pd.concat([a2,pdb_dat],axis=1)
a_n=a_m[a_m['pdb_site'] != 'NONE']
a_n.index=range(len(a_n))

c=pd.read_csv('Unreviewed_v1_single_1.tsv',header=0,sep='\t')
pdb_ids=[]
chains=[]
for i in range(len(c)):
	pdb_id=c['PDBfile'][i][:c['PDBfile'][i].index('.')].lower()
	pdb_ids.append(pdb_id)
	mut_chain=c['MutChain'][i][:c['MutChain'][i].index('_')]
	chains.append(mut_chain)
pdb_dat=pd.DataFrame({'pdb_id':pdb_ids})
chains_dat=pd.DataFrame({'mut_chain':chains})
c1=pd.concat([pdb_dat,chains_dat,c['isPI']],axis=1)
c2=c1.drop_duplicates(subset=['pdb_id','isPI'],keep='first')
c2.index=range(len(c2))
c2_pdb=list(set(list(c2['pdb_id'])))
a3=a_n[a_n['pdb_id'].isin(c2_pdb)]
a3.index=range(len(a3))
a4=pd.merge(a3,c2,how='left',on=['pdb_id','mut_chain'])
a4.index=range(len(a4))
#Reviewed:9508
#Unreviewed:356
		
#统计蛋白A的链长,x为数据（即上面a4），y为路径（即PremPS结果所在路径）
def get_seq(x,y):
	A_length=[]
	for i in range(len(x)):
		mutchain_i=x['mut_chain'][i]
		PI_i=x['isPI'][i]
		path_i=y+PI_i+'out'
		all_i=os.listdir(path_i)
		file1=x['pdb_id'][i].upper()+'_'+mutchain_i+'_1.seq'
		file_i=path_i+'/'+file1
		if os.path.exists(file_i):
			with open(file_i) as f:
				file_lines=f.readlines()
				len_all=0
				for j in file_lines:
					if j[0] != '>':
						len_j=len(j)
						len_all=len_all+len_j
			A_length.append(len_all)
		else:
			A_length.append('NO seq')
	A_dat=pd.DataFrame({'A length':A_length})
	x1=pd.concat([x,A_dat],axis=1)
	return x1


#获取PremPS计算结果,x为数据（即上面a4），y为路径（即PremPS结果所在路径）
def get_PremPS(x,y):
	PremPS_all=[]
	location_all=[]
	for i in range(len(x)):
		isPI_i=x['isPI'][i]
		temp_mutation=x['Mutations'][i]
		site_i=x['pdb_site'][i]
		mutation_i=temp_mutation[0]+site_i+temp_mutation[-1]
		path_i=y+isPI_i+'/'
		file_i=path_i+isPI_i+'.sunddg'
		if os.path.exists(file_i):
			all_i=pd.read_csv(file_i,header=0,sep='\t')
			all_mut=list(all_i['Mutation_PDB'])
			if mutation_i in all_mut:
				data_i=all_i[all_i['Mutation_PDB']==mutation_i]
				data_i.index=range(len(data_i))
				PremPS_i=data_i['PremPS'][0]
				location_i=data_i['Location'][0]
				PremPS_all.append(PremPS_i)
				location_all.append(location_i)
			else:
				PremPS_all.append('NO')
				location_all.append('NO')
		else:
			PremPS_all.append('NONE PDB')
			location_all.append('NONE PDB')
	PremPS_dat=pd.DataFrame({'PremPS':PremPS_all})
	location_dat=pd.DataFrame({'Location':location_all})
	x1=pd.concat([x,PremPS_dat,location_dat],axis=1)
	return x1

path4='/data/liuyang/Reviewed_PremPS/'
a5=get_seq(a4,path4)
a5_qu=a5[a5['A length']=='NO seq']
qu_index=list(a5_qu.index)
a6=a5.drop(qu_index)
a6.index=range(len(a6))

path5='/data/liuyang/Reviewed_PremPS/'
a7=get_PremPS(a6,path5)

a7.to_csv('Reviewed_PremPS_before_result.tsv',index=False,sep='\t')

#保留可用结果
import os
os.chdir('/data/liuyang')
import pandas as pd
a=pd.read_csv('Unreviewed_PremPS_before_result.tsv',header=0,sep='\t')
a1=a[a['PremPS']=='NO']
a2=a[a['PremPS']=='NONE PDB']
qu_index=list(a1.index)+list(a2.index)
a_result=a.drop(qu_index)
a_result.index=range(len(a_result))
a_result.to_csv('Unreviewed_PremPS_result.tsv',index=False,sep='\t')

#Reviewed:8986
#Unreviewed:337
#S4191:4132


#分析结果
import os
os.chdir('/data/liuyang')
import pandas as pd
a=pd.read_csv('Reviewed_PremPS_result.tsv',header=0,sep='\t')
b=pd.read_csv('Unreviewed_PremPS_result.tsv',header=0,sep='\t')
c=pd.read_csv('S4191_PremPS_result.tsv',header=0,sep='\t')
a1=pd.concat([a['type'],a['PremPS']],axis=1)
b1=pd.concat([b['type'],b['PremPS']],axis=1)
data=pd.concat([a1,b1])

pp=[]
for i in range(len(c)):
    if c['DDGexp'][i]>0:
        pp.append('destabilizing')
    else:
        pp.append('stabilizing')
pp1=pd.DataFrame({'type':pp})
c1=pd.concat([pp1,c['PremPS']],axis=1)

data.index=range(len(data))
data1=data[data['type']=='No effect']
dindex=list(data1.index)
data2=data.drop(dindex)
data1=pd.concat([data2,c1])
data1.index=range(len(data1))

t1=['Disrupting','Decreasing','destabilizing']
t2=['Increasing','stabilizing']
for j in range(len(data1)):
    if data1['type'][j] in t1:
        data1['type'][j]='destabilizing'
    elif data1['type'][j] in t2:
        data1['type'][j]='stabilizing'


import numpy as np
import matplotlib.pyplot as plt

data_m=data1[data1['type']=='destabilizing']
data_n=data1[data1['type']=='stabilizing']
a2=np.array(list(data_m['PremPS']))
c2=np.array(list(data_n['PremPS']))

import numpy as np
# 绘制箱线图
fig, ax = plt.subplots()
ax.boxplot([a2,c2],showfliers=True, showmeans=True)
# 设置图表标题和轴标签
ax.set_title('Boxplot of IntAct mutations datasets',fontsize=14)
ax.set_xlabel('Mutation type',fontsize=14)
ax.set_ylabel('PremPS',fontsize=14)
ax.set_xticklabels(['destabilizing', 'stabilizing'],fontsize=14)
# 显示图表
plt.show()

from scipy.stats import ttest_ind
t_statistic, p_value = ttest_ind(a2, c2)
p_value

#绘制密度分布图
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='white',font_scale=1.5)
g = sns.kdeplot(data_m['PremPS'],shade=False,color='red')
g = sns.kdeplot(data_n['PremPS'],shade=False,color='blue')
plt.legend(['destabilizing','stabilizing'],fontsize=15,bbox_to_anchor=(1.4, 1),loc='upper right')
plt.title('Density plot of 2 mutation type PremPS ΔΔG',fontsize=15)
threshold = 0.5
plt.axvline(x=threshold, color='r', linestyle='--')
# 显示图形
plt.show()

data_m.index=range(len(data_m))
data_n.index=range(len(data_n))
pred1=[]
for i in range(len(data_m)):
    if data_m['PremPS'][i]>0:
        pred1.append(1)
    else:
        pred1.append(0)
pred2=[]
for j in range(len(data_n)):
    if data_n['PremPS'][j]<0:
        pred2.append(1)
    else:
        pred2.append(0)
pred1_dat=pd.DataFrame({'predict':pred1})
pred2_dat=pd.DataFrame({'predict':pred2})
data_m1=pd.concat([data_m,pred1_dat],axis=1)
data_n1=pd.concat([data_n,pred2_dat],axis=1)

data1=pd.concat([data_m1,data_n1])
data1.index=range(len(data1))




#计算其他数值
#判断降低和打破的预测准确率
#判断增加的预测准确率
#判断真阳、假阳、真阴、假阴
down=['destabilizing']
up=['stabilizing']
TP_list=[]
FP_list=[]
TN_list=[]
FN_list=[]
for i in range(len(data1)):
	if data1['type'][i] in down and data1['predict'][i]==1:
		TP_list.append(i)
	elif data1['type'][i] in down and data1['predict'][i]==0:
		FP_list.append(i)
	elif data1['type'][i] in up and data1['predict'][i]==1:
		TN_list.append(i)
	elif data1['type'][i] in up and data1['predict'][i]==0:
		FN_list.append(i)

TP=len(TP_list)
FP=len(FP_list)
TN=len(TN_list)
FN=len(FN_list)





import math
accuracy=(TP+TN)/(TP+TN+FP+FN)
specificity=TN/(TN+FP)
recall=TP/(TP+FN)
f1=2*TP/(2*TP+FP+FN)
MCC = (TP * TN - FP * FN) / math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
#计算OR值，看降低的突变是ddg（PremPS）是否倾向于大于0
#看增加的突变是ddg（PremPS）是否倾向于小于0
#使用R进行计算
#构建2*2的列联表，使用公式计算OR值（大，说明降低蛋白质稳定性的突变的ddg更倾向于大于0）
l1=len(m2[m2['predict']==1])
l2=len(m2[m2['predict']==0])
l3=len(data3_1[data3_1['predict']==1])
l4=len(data3_1[data3_1['predict']==0])
OR1=(l1 * l4)/(l2 * l3)

#specificity和MCC马修斯相关系数小

#计算kappa系数
from sklearn.metrics import cohen_kappa_score
y_true=[]
for i in range(len(data_m)):
	if data_m['type'][i] in down:
		y_true.append(1)
	else:
		y_true.append(0)
y_pred=[]
for i in range(len(data_m)):
	if data_m['PremPS'][i]>0:
		y_pred.append(1)
	else:
		y_pred.append(0)
# y_true表示真实标签，y_pred表示预测标签
kappa = cohen_kappa_score(y_true, y_pred)









#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################

#处理成mutbind2形式
import os
os.chdir('/data/liuyang/bishe_data/Mutabind2')
import pandas as pd

def panduan(x):
	from collections import Counter
	x1 = Counter(x)
	y=dict(x1)
	z=''
	for i in sorted(list(y.keys())):
		for j in range(y[i]):
			z=z+'.'+i+'_'+str(j+1)
	z1=z[1:]
	return z1

a1=pd.read_csv('Unreviewed_after.tsv',header=0,sep='\t')
#去除pdb_site为NONE的数据
a_qu=list(a1[a1['pdb_site']=='NONE'].index)
a=a1.drop(a_qu)
a.index=range(len(a))

pdb_ids=[]
chains_id=[]
MutChains=[]
Mutations=[]
result_ids=[]
pdb_chains=[]
A_Chains=[]
B_Chains=[]
for i in range(len(a)):
	pdb_id=a['pdb_id'][i].upper()+'.pdb'
	A_Chain=eval(a['A chains'][i])
	A_Chain_1=panduan(A_Chain)
	B_Chain=eval(a['B chains'][i])
	B_Chain_1=panduan(B_Chain)
	
	MutChain=a['mut_chain'][i]+'_1'
	Mutation=a['Mutations'][i][0]+str(a['pdb_site'][i])+a['Mutations'][i][-1]

	pdb_ids.append(pdb_id)
	A_Chains.append(A_Chain_1)
	B_Chains.append(B_Chain_1)
	MutChains.append(MutChain)
	Mutations.append(Mutation)

pdb_dat=pd.DataFrame({'PDBfile':pdb_ids})
A_dat=pd.DataFrame({'Partner1':A_Chains})
B_dat=pd.DataFrame({'Partner2':B_Chains})
MutChains_dat=pd.DataFrame({'MutChain':MutChains})
Mutations_dat=pd.DataFrame({'Mutation_PDB':Mutations})
a1=pd.concat([pdb_dat,A_dat,B_dat,MutChains_dat,Mutations_dat],axis=1)

#删除其中没有pdb结构的数据
x=pd.read_csv('Unreviewed_v1_single_1.tsv',header=0,sep='\t')
x_pro=list(set(list(x['PDBfile'])))
a_x=a1[a1['PDBfile'].isin(x_pro)]
a_x.index=range(len(a_x))

#同聚体不能用来计算，因此删除同聚体的数据
index_i=[]
for i in range(len(a_x)):
	if a_x['Partner1'][i]==a_x['Partner2'][i]:
		index_i.append(i)
a_y=a_x.drop(index_i)
a_y.index=range(len(a_y))


a2=list(a_y.groupby(['PDBfile','MutChain']))
isPI=[]
result_data=pd.DataFrame()
for j in range(len(a2)):
	temp_m=['ly_n'+str(j+1)]*len(a2[j][1])
	temp_dat=pd.DataFrame({'isPI':temp_m})
	temp_a=a2[j][1]
	temp_a.index=range(len(temp_a))
	temp_data=pd.concat([temp_a,temp_dat],axis=1)
	result_data=pd.concat([result_data,temp_data])

result_data1=result_data
result_data1.index=range(len(result_data1))
temp_re=[x+1 for x in range(len(result_data1))]
result_dat=pd.DataFrame({'Result_Id':temp_re})
result_data2=result_data1
result_data3=pd.concat([result_data2,result_dat],axis=1)

#处理同聚体形式的数据，将突变链作为Partner1，其余链作为Partner2
for i in range(len(result_data3)):
	if result_data3['Partner1'][i]==result_data3['Partner2'][i]:
		result_data3['Partner1'][i]=result_data3['MutChain'][i]
		temp_B=result_data3['Partner2'][i].split('.')
		temp_B.remove(result_data3['MutChain'][i])
		temp_B1='.'.join(temp_B)
		result_data3['Partner2'][i]=temp_B1
	else:
		continue
result_data3.to_csv('Unreviewed_mutabind2.tsv',index=False,sep='\t')







from itertools import groupby
from collections import defaultdict,deque,Counter
from pandas import DataFrame

##设置工作路径
workdir = '/data/yangyan/latest_work'
pathdata = '/data/yangyan/latest_work'

def produce_input():
    usedcols = ['PDBfile','Partner1','Partner2','MutChain','Mutation_PDB','Result_Id','isPI']
    f = pd.read_csv('Unreviewed_mutabind2.tsv', sep='\t')
    f = f.loc[:,usedcols].drop_duplicates(['PDBfile','Partner1','Partner2','MutChain','Mutation_PDB'])
    for _,group in f.groupby(['isPI']):
        jobid = ''.join(list(set(group['isPI'])))
        x=''.join(list(set(group['PDBfile'])))
        pdbfile = x[:x.index('.')].upper()+'.pdb'
#         with open(workdir+jobid+'.sh','w') as fw:
#             fw.write('python premps.py -i %s > %s/%s.log 2>&1\npython premps_cyt.py -i %s > %s/%s.log1 2>&1' % (jobid,jobid,jobid,jobid,jobid,jobid))
        os.system('mkdir %s' % (workdir+jobid))
        os.system('wget https://files.rcsb.org/download/%s -O %s' % (pdbfile,workdir+jobid+'/'+pdbfile))
        group.loc[:,usedcols].to_csv(workdir+jobid+'/'+jobid+'.input',sep='\t',index=0)
    # alljobids,xargs运行需要调用
    with open(workdir+'jobids','w') as fw:
        fw.write('\n'.join(set(f['isPI'])))


nohup xargs -a jobids -P 20 -I {} /home/cheny/anaconda2/bin/python PremPS.py -i {} &
/usr/bin/python MutaBindS.py -i 9shd > 9shd/debug.txt 2>&1
/usr/bin/python MutaBindS.py -i ly_n1 > ly_n1/debug.txt 2>&1


###################################################################################################################################################################
import os
os.chdir('/data1/cheny/ly')
import pandas as pd
path1='/data1/cheny/ly'
data1 = [f for f in os.listdir(path1) if os.path.isdir(os.path.join(path1, f))]
m=[x for x in data1 if x[:3]=='ly_' and x[-3:] != 'out']
n=[y for y in data1 if y[-3:]=='out']
n1=[x[:-3] for x in n]
p=[z for z in m if z not in n1]
p1=list(set(p))

for i in range(len(p1)):
    x='/usr/bin/python MutaBindS.py -i '+ p1[i] + ' > ' +p1[i]+'/debug.txt 2>&1'
    x1=pd.DataFrame({'code':[x]})
    x1.to_csv(p1[i]+'.sh',header=False,index=False)



#读取mutabind2计算结果
import os
os.chdir('/data/liuyang/bishe_data/Mutabind2/Unreviewed_result')
import pandas as pd
a=pd.read_csv('/data/liuyang/bishe_data/Mutabind2/Unreviewed_mutabind2.tsv',header=0,sep='\t')
#只有异聚体
ddgs=[]
for i in range(len(a)):
	mut_i=a['Mutation_PDB'][i]
	path_i='/data/liuyang/bishe_data/Mutabind2/Unreviewed_result/'+a['isPI'][i]
	file_i=os.listdir(path_i)
	if os.path.exists(path_i+'/'+a['isPI'][i]+'.sunddg'):
		datax=pd.read_csv(path_i+'/'+a['isPI'][i]+'.sunddg',header=0,sep='\t')
		mutx=list(datax['Mutation_PDB'])
		if mut_i in mutx:
			x_i=datax[datax['Mutation_PDB']==mut_i]
			x_i.index=range(len(x_i))
			ddg_i=x_i['SunDDG'][0]
			ddgs.append(ddg_i)
		else:
			ddgs.append('NONE MUT')#储存突变位点没有匹配的
	else:
		ddgs.append('NONE PDB')#储存PDB结构不能运行的

ddg_dat=pd.DataFrame({'Mutabind2':ddgs})
a1=pd.concat([a,ddg_dat],axis=1)

pdbs=[]
chains=[]
mut=[]
for i in range(len(a1)):
	pdb_i=a1['PDBfile'][i][:a1['PDBfile'][i].index('.')].lower()
	pdbs.append(pdb_i)
	chain_i=a1['MutChain'][i][:a1['MutChain'][i].index('_')]
	chains.append(chain_i)
	mut_i=a1['Mutation_PDB'][i]
	mut.append(mut_i)
pdbs_dat=pd.DataFrame({'pdb_id':pdbs})
chains_dat=pd.DataFrame({'mut_chain':chains})
muts_dat=pd.DataFrame({'Mutations':mut})
a2=pd.concat([a1,pdbs_dat,chains_dat,muts_dat],axis=1)

b=pd.read_csv('Unreviewed_after_PDBchains.tsv',header=0,sep='\t')
muts=[]
for j in range(len(b)):
	mut_j=b['Mutations'][j][0]+str(b['pdb_site'][j])+b['Mutations'][j][-1]
	muts.append(mut_j)
mut_dat=pd.DataFrame({'Mutations':muts})
b1=pd.concat([b['pdb_id'],b['type'],b['mut_chain'],mut_dat],axis=1)

c=pd.merge(a2,b1,how='left',on=['pdb_id','mut_chain','Mutations'])
c.to_csv('Unreviewed_mutbind2_result.tsv',index=False,sep='\t')


##########################################################
#最终对mutabind2的预测效果进行评估
import os
os.chdir('/data/liuyang/bishe_data/Mutabind2')
import pandas as pd
a=pd.read_csv('Reviewed_mutabind2_result.tsv',header=0,sep='\t')
b=pd.read_csv('Unreviewed_mutabind2_result.tsv',header=0,sep='\t')
c=pd.read_csv('S4191.txt',header=0,sep='\t')
#去除a和b中mutabind2为NONE的数据
a1=a[a['Mutabind2']=='NONE PDB']
a2=a[a['Mutabind2']=='NONE MUT']
index1=list(a1.index)+list(a2.index)
ax=a.drop(index1)
ax.index=range(len(ax))
b1=b[b['Mutabind2']=='NONE PDB']
b2=b[b['Mutabind2']=='NONE MUT']
index2=list(b1.index)+list(b2.index)
bx=b.drop(index2)
bx.index=range(len(bx))
ax1=pd.concat([ax['type'],ax['Mutabind2']],axis=1)
bx1=pd.concat([bx['type'],bx['Mutabind2']],axis=1)
x=pd.concat([ax1,bx1])

t1=['Disrupting','Decreasing','destabilizing']
t2=['Increasing','stabilizing']
x_no=x[x['type']=='No effect']
no_index=list(x_no.index)
x1=x.drop(no_index)
x1.index=range(len(x1))
for j in range(len(x1)):
    if x1['type'][j] in t1:
        x1['type'][j]='destabilizing'
    elif x1['type'][j] in t2:
        x1['type'][j]='stabilizing'

c1=pd.concat([c['DDGexp'],c['MutaBind2']],axis=1)
types=[]
for i in range(len(c1)):
	if c1['DDGexp'][i]>0:
		type_i='destabilizing'
		types.append(type_i)
	else:
		type_i='stabilizing'
		types.append(type_i)
c_types=pd.DataFrame({'type':types})
c2=pd.concat([c_types,c1['MutaBind2']],axis=1)
c3=c2.rename(columns={'MutaBind2':'Mutabind2'})

#总数据
data=pd.concat([x1,c3])
data.index=range(len(data))
data['Mutabind2'] = data['Mutabind2'].astype(float)  #将所有Mutabind2字符串转换成float型

data.to_csv('Mutabind2_all.tsv',index=False,sep='\t')



#绘制箱线图
import numpy as np
import matplotlib.pyplot as plt
data_m=data[data['type']=='destabilizing']
data_n=data[data['type']=='stabilizing']
m1=list(data_m['Mutabind2'])
n1=list(data_n['Mutabind2'])
m2=[float(x) for x in m1]
n2=[float(y) for y in n1]
a2=np.array(m2)
c2=np.array(n2)
import numpy as np
# 绘制箱线图
fig, ax = plt.subplots()
ax.boxplot([a2,c2],showfliers=True, showmeans=True)
# 设置图表标题和轴标签
ax.set_title('Boxplot of the Mutabind2 ΔΔG of mutations datasets',fontsize=14)
ax.set_xlabel('Mutation type',fontsize=14)
ax.set_ylabel('Mutabind2 ΔΔG',fontsize=14)
ax.set_xticklabels(['destabilizing', 'stabilizing'],fontsize=14)
# 显示图表
plt.show()

#t检验
from scipy.stats import ttest_ind
t_statistic, p_value = ttest_ind(a2, c2)
p_value

#绘制密度分布图
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='white',font_scale=1.5)
g = sns.kdeplot(data_m['Mutabind2'],shade=False,color='red')
g = sns.kdeplot(data_n['Mutabind2'],shade=False,color='blue')
plt.legend(['destabilizing','stabilizing'],fontsize=15,bbox_to_anchor=(1.4, 1),loc='upper right')
plt.title('Density plot of 2 mutation type Mutabind2 ΔΔG',fontsize=15)
threshold = 0.5
plt.axvline(x=threshold, color='r', linestyle='--')
# 显示图形
plt.show()

data_m.index=range(len(data_m))
data_n.index=range(len(data_n))
pred1=[]
for i in range(len(data_m)):
    if data_m['Mutabind2'][i]>0:
        pred1.append(1)
    else:
        pred1.append(0)
pred2=[]
for j in range(len(data_n)):
    if data_n['Mutabind2'][j]<0:
        pred2.append(1)
    else:
        pred2.append(0)

pred1_dat=pd.DataFrame({'predict':pred1})
pred2_dat=pd.DataFrame({'predict':pred1})
data_m1=pd.concat([data_m,pred1_dat],axis=1)
data_n1=pd.concat([data_n,pred2_dat],axis=1)

data1=pd.concat([data_m1,data_n1])
#计算其他数值
#判断降低和打破的预测准确率
#判断增加的预测准确率
#判断真阳、假阳、真阴、假阴
down=['destabilizing']
up=['stabilizing']
TP_list=[]
FP_list=[]
TN_list=[]
FN_list=[]
for i in range(len(data1)):
	if data1['type'][i] in down and data1['predict'][i]==1:
		TP_list.append(i)
	elif data1['type'][i] in down and data1['predict'][i]==0:
		FP_list.append(i)
	elif data1['type'][i] in up and data1['predict'][i]==1:
		TN_list.append(i)
	elif data1['type'][i] in up and data1['predict'][i]==0:
		FN_list.append(i)

TP=len(TP_list)
FP=len(FP_list)
TN=len(TN_list)
FN=len(FN_list)


import math
accuracy=(TP+TN)/(TP+TN+FP+FN)
specificity=TN/(TN+FP)
recall=TP/(TP+FN)
f1=2*TP/(2*TP+FP+FN)
MCC = (TP * TN - FP * FN) / math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
#计算OR值，看降低的突变是ddg（PremPS）是否倾向于大于0
#看增加的突变是ddg（PremPS）是否倾向于小于0
#使用R进行计算
#构建2*2的列联表，使用公式计算OR值（大，说明降低蛋白质稳定性的突变的ddg更倾向于大于0）

l1=len(data_m1[data_m1['predict']==1])
l2=len(data_m1[data_m1['predict']==0])
l3=len(data_n1[data_n1['predict']==1])
l4=len(data_n1[data_n1['predict']==0])
OR1=(l1 * l4)/(l2 * l3)

#specificity和MCC马修斯相关系数小

#计算kappa系数
from sklearn.metrics import cohen_kappa_score
y_true=[]
for i in range(len(data_m)):
	if data_m['type'][i] in down:
		y_true.append(1)
	else:
		y_true.append(0)
y_pred=[]
for i in range(len(data_m)):
	if data_m['Mutabind2'][i]>0:
		y_pred.append(1)
	else:
		y_pred.append(0)
# y_true表示真实标签，y_pred表示预测标签
kappa = cohen_kappa_score(y_true, y_pred)

###################################################################################################################################################################
###################################################################################################################################################################
###################################################################################################################################################################
###################################################################################################################################################################
###################################################################################################################################################################
###################################################################################################################################################################
'''
index_all=[]
for p in range(len(a_all)):
	for q in a_all['A all chains'][p]:
		if a_all['A all chains'][p].count(q)>1:
			index_all.append(p)
index_all_1=list(set(index_all))
#a1:受影响蛋白A，链分成两段以上的
a1=a_all.loc[index_all_1]
a1.index=range(len(a1))
#a2：受影响蛋白A，链未分成两段的
#首先找出所有符合情况的链（包含突变位点）
a2=a_all.drop(index_all_1)
a2.index=range(len(a2))
column1=list(a2.columns.values)
all_num=int(len([y for y in column1 if y[-1].isdigit()])/9)
may_mut_chain=[]
for m in range(len(a2)):
	may_A=[]
	loc1=int(re.findall(r'\d+',a2['Mutations'][m])[0])
	if len(set(eval(a2['A coverage'][m]))) != 1:
		for n in eval(a['A chains'][m]):
			for p in range(all_num):
				if a2['chain_id_'+str(p)][m]==n and loc1>=a2['unp_start_'+str(p)][m] and loc1<=a2['unp_end'+str(p)][m]:
					may_A.append(a2['chain_id_'+str(p)][m])
					may_A1=list(set(may_A))
	else:
		may_A=eval(a2['A chains'][m])
		may_A1=list(set(may_A))
	may_mut_chain.append(may_A1)
may_dat=pd.DataFrame({'may_mutchain':may_mut_chain})
a2_1=pd.concat([a2,may_dat],axis=1)
#如果蛋白A只有一条链，直接去预测
index1=[]
for i in range(len(a2_1)):
	if len(a2_1['may_mutchain'][i])==1:
		index1.append(i)
a2_single=a2_1.loc[index1]
a2_single.index=range(len(a2_single))
a2_single.to_csv('Reviewed_PremPS_single.tsv',index=False,sep='\t',na_rep='_NAN_')



#data1.to_csv('Reviewed_map_415_1.tsv',index=False,sep='\t',na_rep='_NAN_')
#xxx=pd.read_csv('Reviewed_map_415_1.tsv',header=0,sep='\t',keep_default_na=False, na_values=['_NAN_'])
#因为有链的id为'NA'，pandas会自动将其读取成NaN
#因此用na_rep参数将真正的NaN先保存为'_NAN_'，读取的时候再退回来






index1=[]
for i in range(len(a)):
	temp_chains=eval(a['A chains'][i])
	if len(temp_chains)==1:
		index1.append(i)
ok_data=a.loc[index1]
ok_data.index=range(len(ok_data))
column1=list(ok_data.columns.values)
all_num=int(len([y for y in column1 if y[-1].isdigit()])/9)
index2=[]
for m in range(len(ok_data)):
	start=[]
	end=[]
	for n in range(all_num):
		if ok_data['Uniprot_id_'+str(n)][m]==ok_data['Affected protein AC'][m]:
			start.append(int(float(ok_data['unp_start_'+str(n)][m])))
			end.append(int(float(ok_data['unp_end_'+str(n)][m])))
	loc1=int(re.findall(r'\d+',ok_data['Mutations'][m])[0])
	if max(start)==min(end)+1 and loc1>=min(start) and loc1<=max(end):
		index2.append(m)
	elif max(start)!=min(end)+1:
		if loc1>min(start) and loc1<min(end):
			index2.append(m)
		elif loc1>max(start) and loc1<max(end):
			index2.append(m)
ok_data1=ok_data.loc[index2]
ok_data1.to_csv('Reviewed_PremPS.tsv',index=False,sep='\t')




#以Reviewed_PremPS_single.tsv为例，将数据移入新路径/data/liuyang/bishe_data/PremPS
import re
import os
os.chdir('/data/liuyang/bishe_data/PremPS1')
import pandas as pd
a=pd.read_csv('Reviewed_PremPS_single.tsv',header=0,sep='\t',keep_default_na=False, na_values=['_NAN_'])
a1=pd.concat([a['Affected protein AC'],a['pdb_id'],a['A chains'],a['Mutations']],axis=1)
chain_id=[]
locations=[]
for i in range(len(a1)):
	chain_id.append(eval(a1['A chains'][i])[0])
	loc1=int(re.findall(r'\d+',a1['Mutations'][i])[0])
	locations.append(loc1)
chain_dat=pd.DataFrame({'chain_id':chain_id})
loc_dat=pd.DataFrame({'unp_sites':locations})
a2=a1.drop('A chains',axis=1)
a3=pd.concat([a2,chain_dat],axis=1)
a4=a3.drop('Mutations',axis=1)
a5=pd.concat([a4,loc_dat],axis=1)
a6=list(a5.groupby(['Affected protein AC','pdb_id','chain_id']))
unp_sites_1=[]
pro=[]
pdb=[]
chains=[]
for j in range(len(a6)):
	mut_sites=list(a6[j][1]['unp_sites'])
	mut_set1=list(set(mut_sites))
	mut_set=sorted(mut_set1)
	result = ','.join(str(x) for x in mut_set)
	unp_sites_1.append(result)
	pro.append(a6[j][0][0])
	pdb.append(a6[j][0][1])
	chains.append(a6[j][0][2])
a7=a5.drop('unp_sites',axis=1)
unp_data_1=pd.DataFrame({'unp_sites':unp_sites_1})
pro_dat=pd.DataFrame({'unp_id':pro})
pdb_dat=pd.DataFrame({'pdb_id':pdb})
chain_dats=pd.DataFrame({'chain_id':chains})
a8=pd.concat([pro_dat,pdb_dat,chain_dats,unp_data_1],axis=1)
a8.to_csv('Reviewed_unp_single.tsv',index=False,sep='\t')

#学姐的pdb_profiling工具，将uniprot序列位点，转换为pdb位点
#以Reviewed数据为例
import re
import os
os.chdir('/data/liuyang/bishe_data/Reviewed_method1/2')
import pandas as pd
a=pd.read_csv('Reviewed_PremPS_single.tsv',header=0,sep='\t')
b=pd.read_csv('reviewed_unp_first_230418.txt',header=0,sep='\t')
pro_pdb=dict()
for i in range(len(b)):
	temp_pro=b['uni_list'][i].split(',')
	temp_pdb=b['pdb_list'][i].split(',')
	for j in range(len(temp_pro)):
		pro_pdb[temp_pro[j]]=temp_pdb[j]
pdb_sites=[]
for m in range(len(a)):
	loc_m=int(re.findall(r'\d+',a['Mutations'][m])[0])
	temp_mutation=a['Affected protein AC'][m]+'_'+str(loc_m)
	if pro_pdb[temp_mutation] != 'NONE':
		pdb_site=str(pro_pdb[temp_mutation][pro_pdb[temp_mutation].index('_',pro_pdb[temp_mutation].index('_')+1)+1:])
		pdb_sites.append(pdb_site)
	else:
		pdb_sites.append('NONE')
pdb_dat=pd.DataFrame({'pdb_site':pdb_sites})
a1=pd.concat([a,pdb_dat],axis=1)
a1.to_csv('Reviewed_1_single.tsv',index=False,sep='\t')

#注意pdb中有的位点，可能会包含字母，如28G，因此使用字符串格式


#处理数据为PremPS运行格式
import os
os.chdir('/data/liuyang/bishe_data/Reviewed_method1/2')
import pandas as pd

def panduan(x,z):
	from collections import Counter
	x1 = Counter(x)
	y=dict(x1)
	for i in sorted(list(y.keys())):
		for j in range(y[i]):
			z=z+'.'+i+'_'+str(j+1)
	return z

a1=pd.read_csv('Reviewed_1_single.tsv',header=0,sep='\t')
#去除pdb_site为NONE的数据
a_qu=list(a1[a1['pdb_site']=='NONE'].index)
a=a1.drop(a_qu)
a.index=range(len(a))
pdb_ids=[]
chains_id=[]
MutChains=[]
Mutations=[]
result_ids=[]
pdb_chains=[]
for i in range(len(a)):
	pdb_id=a['pdb_id'][i].upper()+'.pdb'
	Chains=eval(a['A chains'][i])[0]+'_1'
	temp_chain=eval(a['B chains'][i])
	temp_Chains=panduan(temp_chain,Chains)
	MutChain=eval(a['A chains'][i])[0]+'_1'
	Mutation=a['Mutations'][i][0]+str(a['pdb_site'][i])+a['Mutations'][i][-1]
	pdb_chain=a['pdb_id'][i]+'_'+eval(a['A chains'][i])[0].lower()
	pdb_ids.append(pdb_id)
	chains_id.append(temp_Chains)
	MutChains.append(MutChain)
	Mutations.append(Mutation)
	pdb_chains.append(pdb_chain)
pdb_dat=pd.DataFrame({'PDBfile':pdb_ids})
chains_dat=pd.DataFrame({'Chains':chains_id})
MutChains_dat=pd.DataFrame({'MutChain':MutChains})
Mutations_dat=pd.DataFrame({'Mutation_PDB':Mutations})
pdb_chains_dat=pd.DataFrame({'pdb_chains':pdb_chains})
a1=pd.concat([pdb_dat,chains_dat,MutChains_dat,Mutations_dat,pdb_chains_dat],axis=1)
a2=list(a1.groupby('PDBfile'))
isPI=[]
result_data=pd.DataFrame()
for j in range(len(a2)):
	temp_m=['m'+str(j+1)]*len(a2[j][1])
	temp_dat=pd.DataFrame({'isPI':temp_m})
	temp_a=a2[j][1]
	temp_a.index=range(len(temp_a))
	temp_data=pd.concat([temp_a,temp_dat],axis=1)
	result_data=pd.concat([result_data,temp_data])

result_data1=result_data.drop_duplicates(subset=['PDBfile','MutChain','Mutation_PDB','pdb_chains'],keep='first')
result_data1.index=range(len(result_data1))
temp_re=[x for x in range(len(result_data1))]
result_dat=pd.DataFrame({'Result_ID':temp_re})
result_data2=result_data1
result_data3=pd.concat([result_data2,result_dat],axis=1)
#去除同聚体造成的重复
temp_bs=[]
for p in range(len(result_data3)):
	temp_b=result_data3['Chains'][p].split('.')
	temp_b1=list(set(temp_b))
	temp_b2='.'.join(temp_b1)
	temp_bs.append(temp_b2)
temp_bs_data=pd.DataFrame({'Chains':temp_bs})
result_data4=result_data3.drop('Chains',axis=1)
result_data5=pd.concat([result_data4,temp_bs_data],axis=1)

result_data5.to_csv('Reviewed_v1_single.tsv',index=False,sep='\t')
#接着使用PremPS进行计算
'''

###########################################################################################################################################################
###########################################################################################################################################################
##########################################################################################################################################
#################################################################################################
########################################################################################################################################################
#################################################################################################################################
############################################################################################################################








###################################################################################################################################################
###################################################################################################################################################
'''
#方案2:
#将相互作用所有相互匹配的对，展开
import os
os.chdir('/data/liuyang/bishe_data/get_pdb/Reviewed/map')
import pandas as pd
import math
a=pd.read_csv('Reviewed.tsv',header=0,sep='\t')
b=pd.read_csv('pdb_ligand.tsv',header=0,sep='\t')
pro_all=[]
for i in range(len(a)):
	pro_a=a['Affected protein AC'][i]
	pro_b=a['b id'][i]
	pro_list=sorted([pro_a,pro_b])
	pro_str=pro_list[0]+'_'+pro_list[1]
	pro_all.append(pro_str)
pro_dat=pd.DataFrame({'pro_str':pro_all})
a1=pd.concat([a,pro_dat],axis=1)

pro_lists=[]
for m in range(len(b)):
	numbers=(b.shape[1]-3)/9
	pro_list=[]
	for n in range(int(numbers)):
		pro_n=b['Uniprot_id_'+str(n)][m]
		chain_n=b['chain_id_'+str(n)][m]
		if pd.isnull(pro_n) or pd.isnull(chain_n):
			break
		else:
			pro_n1=str(pro_n)+'_'+str(chain_n)
			pro_list.append(pro_n1)
	pro_list1=list(set(pro_list))
	pro_lists.append(pro_list1)
pro_dat1=pd.DataFrame({'pro_list':pro_lists})
b1=pd.concat([b,pro_dat1],axis=1)

#将不是复合物的数据去掉
b1_mono=b1[pd.isnull(b1['Uniprot_id_1'])]
mono_index=list(b1_mono.index)
b2=b1.drop(mono_index)
b2.index=range(len(b2))
#将蛋白列表pro_list拆分成一个个组合
for p in range(len(b2)):
	temp_list=sorted(b2['pro_list'][p])
	if len(temp_list)==1:
		pro_p=temp_list[0][:temp_list[0].index('_')]+'_'+temp_list[0][:temp_list[0].index('_')]
	else:
		for q in range(len(temp_list)):
'''















