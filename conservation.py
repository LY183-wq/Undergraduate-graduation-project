# -*- coding: utf-8 -*-
# @Author: Liu Yang
# @Date:   2023-03-16 15:36:02
# @Last Modified by:   Liu Yang
# @Last Modified time: 2023-03-17 15:08:10


#########################################################################################################################################
#########################################################################################################################################
#########################################################################################################################################



###########################################################################################################################################
#Reviewed
import os
os.chdir('/data/liuyang/bishe_data/alphafold/Reviewed')
import pandas as pd
import re
a=pd.read_csv('Reviewed_completed.tsv',header=0,sep='\t')
#按Affected protein进行分组
a1=list(a.groupby('Affected protein AC'))
dat1=pd.DataFrame()
for i in range(len(a1)):
	b=a1[i][1]
	b.index=range(len(b))
	for j in range(len(b)):
		b1=int(re.findall(r'\d+',b['Mutations'][j])[0])
		b2=b['Mutations'][j][0]
		b3=b['Affected protein AC'][j]
		b4=b['Interaction AC'][j]
		b5=b['b id'][j]
		b_dic={'Affected protein AC':[b3],'location':[b1],'Original sequence':[b2],'Interaction AC':[b4],'B_ID':[b5]}
		b_dat=pd.DataFrame(b_dic)
		dat1=pd.concat([dat1,b_dat])
dat1.to_csv('Affected_protein_loc.tsv',index=False,sep='\t')



#获取序列
b=pd.read_csv('Affected_protein_loc.tsv',header=0,sep='\t')
b1=list(b['Affected protein AC'])
b2=list(set(b1))
b3=pd.DataFrame({'Affected protein AC':b2})
b3.to_csv('protein_id.txt',index=False,sep='\t')



######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################

#使用Provean计算序列位点保守性
import os
os.chdir('/data/liuyang/bishe_data/conservations/Reviewed')
import pandas as pd
a=pd.read_csv('Reviewed_completed.tsv',header=0,sep='\t')
group1=[]
for i in range(len(a)):
	temp_group=a['Affected protein AC'][i]+'_'+a['b id'][i]
	group1.append(temp_group)
group_dat=pd.DataFrame({'group':group1})
a1=pd.concat([a,group_dat],axis=1)
##################################################################################################
#首先去除重复冲突数据
a_list=list(a1.groupby(['group','Mutations','Affected protein AC']))
temp_data=pd.DataFrame()
for j in range(len(a_list)):
	if len(a_list[j][1])>1:
		temp_data=pd.concat([temp_data,a_list[j][1]])
temp_index=list(temp_data.index)

#筛选出冲突数据
a_m=a1.loc[temp_index]
a_n=a1.drop(temp_index)
a_m.index=range(len(a_m))
a_n.index=range(len(a_n))
am_list=list(a_m.groupby(['group','Mutations','Affected protein AC']))
chongtu_type1=set(['Decreasing'])
chongtu_type2=set(['No effect'])
chongtu_type3=set(['Increasing'])
chongtu_type4=set(['Disrupting'])
chongtu_type5=set(['Decreasing','Disrupting'])
chongtu_type_all=[chongtu_type1,chongtu_type2,chongtu_type3,chongtu_type4,chongtu_type5]
chongtu_data=pd.DataFrame()
for p in range(len(am_list)):
	temp_list=list(am_list[p][1]['type'])
	temp_set=set(temp_list)
	if temp_set not in chongtu_type_all:
		chongtu_data=pd.concat([chongtu_data,am_list[p][1]])
chongtu_index=list(chongtu_data.index)
am_chongtu=a_m.loc[chongtu_index]
am_zc=a_m.drop(chongtu_index)
zc_data=pd.concat([a_n,am_zc])
am_chongtu.to_csv('chongtu_data.tsv',index=False,sep='\t')
zc_data.to_csv('Reviewed.tsv',index=False,sep='\t')


#################################################################################
#使用provean

#安装各种软件和库，先去下载

#设置环境变量
#设置三个软件的路径，数据库路径，安装路径prefix（因为没有bin安装权限）
./configure PSIBLAST=/data/liuyang/bishe_data/conservations/Reviewed/ncbi-blast-2.4.0+/bin CDHIT=/data/liuyang/bishe_data/conservations/Reviewed/cd-hit-v4.6.8-2017-1208 BLASTDBCMD=/data/liuyang/bishe_data/conservations/Reviewed/ncbi-blast-2.4.0+/bin BLAST_DB=/data/liuyang/bishe_data/conservations/Reviewed/nr --prefix=/data/liuyang/bishe_data/conservations/Reviewed


#################################################################################################################################################
#Reviewed
#首先获取蛋白质fasta序列
import os
os.chdir('/data/liuyang/bishe_data/conservations/Reviewed')
import pandas as pd
a=pd.read_csv('Reviewed.tsv',header=0,sep='\t')
a1=a.drop_duplicates(subset=['Affected protein AC'],keep='first')
p_id=a1['Affected protein AC']
p_id.to_csv('p_id.tsv',index=False,sep='\t')

#将序列分别拆分成一个一个的fasta文件
os.chdir('/data/liuyang/bishe_data/conservations/Reviewed/bin/Reviewed')
from Bio import SeqIO
# 打开输入fasta文件
input_file = "Reviewed_seq.fasta"
# 读取输入fasta文件中的序列，并将每个序列写入单独的文件
for record in SeqIO.parse(input_file, "fasta"):
    output_file = "{}.fasta".format(record.id[record.id.index('|')+1:record.id.index('|',record.id.index('|')+1)])
    with open(output_file, "w") as output_handle:
        SeqIO.write(record, output_handle, "fasta")
print("拆分完毕！")

#构建突变文件
#条件：将所有类型氨基酸残基突变成gly甘氨酸，查看score得分，小于-2.5
import os
os.chdir('/data/liuyang/bishe_data/conservations/Reviewed/bin/Reviewed')
import pandas as pd
a=pd.read_csv('Reviewed.tsv',header=0,sep='\t')
path='/data/liuyang/bishe_data/conservations/Reviewed/bin/Reviewed'
dir_all=os.listdir(path)
#函数change_to_gly，将所有突变的突变后残基换成甘氨酸GLY
def change_to_gly(x):
	y=x[:-1]+'G'
	return y
for i in dir_all:
	if i[-5:] == 'fasta':
		a_i=a[a['Affected protein AC']==i[:i.index('.')]]
		a_i_mut=list(a_i['Mutations'])
		mut_i=[]
		for j in range(len(a_i_mut)):
			mut_j=change_to_gly(a_i_mut[j])
			mut_i.append(mut_j)
			mut_i_1=list(set(mut_i))
		mut_i_dat=pd.DataFrame({'mutations':mut_i_1})
		mut_i_dat.to_csv(i[:i.index('.')]+'.var',index=False,header=0)

#运行provean.sh来获取结果	
import os
os.chdir('/data/liuyang/bishe_data/conservations/Reviewed/bin/Reviewed')
import pandas as pd	
path='/data/liuyang/bishe_data/conservations/Reviewed/bin/Reviewed'
dir_all=os.listdir(path)
for i in dir_all:
	if i[-5:] == 'fasta':
		os.system('provean.sh -q '+i+' -v '+i[:i.index('.')]+'.var'+' >'+i[:i.index('.')]+'.out'+' --num_threads 20')
		print(dir_all.index(i))

#由于之前的数据统计错误，因此重新获取provean score
'''
#从out格式的输出文件中获取score信息，小于-2.5为保守残基
import os
os.chdir('/data/liuyang/bishe_data/conservations/Reviewed/bin/Reviewed')
import pandas as pd	
path='/data/liuyang/bishe_data/conservations/Reviewed/bin/Reviewed'
dir_all=os.listdir(path)
out_all=[x for x in dir_all if x[-3:]=='out']
scores=[]
a=pd.read_csv('Reviewed.tsv',header=0,sep='\t')
#获取scores分值
for i in range(len(a)):
	a_p=a['Affected protein AC'][i]
	a_m=a['Mutations'][i][:-1]+'G'
	b=pd.read_csv(a_p+'.out',header=0,sep='\t')
	b1=b.T
	score_i=b1[a_m][0]
	scores.append(score_i)
scores_dat=pd.DataFrame({'provean_score':scores})
a1=pd.concat([a,scores_dat],axis=1)

#判断该位点是否保守
conservations=[]
import re
for j in range(len(a1)):
	if float(a1['provean_score'][j]) < -2.5:
		conservations.append(1)
	else:
		conservations.append(0)
conservations_dat=pd.DataFrame({'conservations':conservations})
a2=pd.concat([a1,conservations_dat],axis=1)
a2.to_csv('Reviewed_provean_score.tsv',index=False,sep='\t')
'''

#四种类型突变数据的绘制密度分布图
import os
os.chdir('/data/liuyang/bishe_data/conservations')
a=pd.read_csv('Reviewed_provean_score.tsv',header=0,sep='\t')
a1=a[a['type']=='Disrupting']
a1.index=range(len(a1))
a2=a[a['type']=='Decreasing']
a2.index=range(len(a2))
a3=a[a['type']=='Increasing']
a3.index=range(len(a3))
a4=a[a['type']=='No effect']
a4.index=range(len(a4))

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='ticks')
g = sns.kdeplot(a1['provean_score'],shade=False,color='red')
g = sns.kdeplot(a2['provean_score'],shade=False,color='orange')
g = sns.kdeplot(a3['provean_score'],shade=False,color='blue')
g = sns.kdeplot(a4['provean_score'],shade=False,color='green')

plt.legend(['Disrupting', 'Decreasing', 'Increasing','No effect'])
plt.title('Density plot of 4 types Reviewed mutations data')

# 显示图形
plt.show()

##########################################################################################################################

import os
os.chdir('/data/liuyang/bishe_data/conservations')
a2=pd.read_csv('Reviewed_provean_score.tsv',header=0,sep='\t')
#绘制保守性分数折线图，结果不清晰
scores_1=list(a2['provean_score'])
scores_2=[float(x) for x in scores_1]
index1=list(range(len(scores_2)))
import matplotlib.pyplot as plt
plt.plot(index1,scores_2)
plt.xlabel('Mutation sites index')
plt.ylabel('provean conservation scores')
plt.title('Provean conservation scores distribution in Reviewed data')



#绘制每种类型的柱状图
import os
os.chdir('/data/liuyang/bishe_data/conservations')
import pandas as pd	
x=pd.read_csv('Reviewed_provean_score.tsv',header=0,sep='\t')
a=x[x['type']=='Disrupting']
b=x[x['type']=='Decreasing']
c=x[x['type']=='Increasing']
d=x[x['type']=='No effect']
a1=len(a[a['conservations']==1])
a2=len(a[a['conservations']==0])
b1=len(b[b['conservations']==1])
b2=len(b[b['conservations']==0])
c1=len(c[c['conservations']==1])
c2=len(c[c['conservations']==0])
d1=len(d[d['conservations']==1])
d2=len(d[d['conservations']==0])
aa=[a1,a2]
bb=[b1,b2]
cc=[c1,c2]
dd=[d1,d2]
xlabs=['conservative','not conservative']
import matplotlib.pyplot as plt
width=0.2
x1=[1,3]
x2=[1.2,3.2]
x3=[1.4,3.4]
x4=[1.6,3.6]
plt.bar(x1, aa,width=width, label='Disrupting',color='red')  # 绘制第一组数据的柱形
plt.bar(x2, bb, width, label='Decreasing',color='orange')  # 绘制第二组数据的柱形
plt.bar(x3, cc, width, label='Increasing',color='blue')  # 绘制第二组数据的柱形
plt.bar(x4, dd, width, label='No effect',color='green')  # 绘制第二组数据的柱形
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
# 添加标签和标题
plt.xlabel('type')
plt.ylabel('numbers')
plt.title('conservative numbers Reviewed data')
plt.xticks([1.3,3.3],xlabs)#设置横坐标
for m,n in zip(x1,aa):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=7)
for m,n in zip(x2,bb):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=7)
for m,n in zip(x3,cc):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=7)
for m,n in zip(x4,dd):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=7)


#绘制柱状图
a2_con=len(a2[a2['conservations']==1])
a2_uncon=len(a2[a2['conservations']==0])
a2_y=[a2_con,a2_uncon]
a2_x=['conservative','not conservative']
color_1=['red','blue']
import matplotlib.pyplot as plt
plt.bar(a2_x,a2_y,color=color_1)
plt.xlabel('Type')
plt.ylabel('Numbers')
plt.title('Conservation type in Reviewed data')
for m,n in zip(a2_x,a2_y):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=7)

#绘制饼图
data=[a2_con/len(a2),a2_uncon/len(a2)]
label=['conservative','not conservative']
color=['orange','green']
import matplotlib.pyplot as plt
plt.pie(data,labels=label,colors=color,shadow=True,autopct="(%1.1f%%)")
# 添加标签和标题
plt.title('Conservation type in Reviewed data')
# 保存图片
plt.savefig('Reviewed_Disrupting.png')









#################################################################################################################################################
#Unreviewed
#首先获取蛋白质fasta序列
import os
os.chdir('/data/liuyang/bishe_data/conservations/Unreviewed')
import pandas as pd
a=pd.read_csv('Unreviewed.tsv',header=0,sep='\t')
a1=a.drop_duplicates(subset=['Affected protein AC'],keep='first')
p_id=a1['Affected protein AC']
p_id.to_csv('Unr_p_id.tsv',index=False,sep='\t')

#将序列分别拆分成一个一个的fasta文件
#将序列分别拆分成一个一个的fasta文件
os.chdir('/data/liuyang/bishe_data/conservations/Reviewed/bin/Unreviewed')
from Bio import SeqIO
# 打开输入fasta文件
input_file = "Unreviewed_seq.fasta"
# 读取输入fasta文件中的序列，并将每个序列写入单独的文件
for record in SeqIO.parse(input_file, "fasta"):
    output_file = "{}.fasta".format(record.id[record.id.index('|')+1:record.id.index('|',record.id.index('|')+1)])
    with open(output_file, "w") as output_handle:
        SeqIO.write(record, output_handle, "fasta")
print("拆分完毕！")

#构建突变文件
#条件：将所有类型氨基酸残基突变y甘氨酸，查看score得分，小于-2.5
import os
os.chdir('/data/liuyang/Provean/Provean_all')
import pandas as pd
a=pd.read_csv('Reviewed.tsv',header=0,sep='\t')
path='/data/liuyang/bishe_data/conservations/Reviewed/bin/Unreviewed'
dir_all=os.listdir(path)
#函数change_to_gly，将所有突变的突变后残基换成甘氨酸GLY
def change_to_gly(x):
	y=x[:-1]+'G'
	return y
for i in dir_all:
	if i[-5:] == 'fasta':
		a_i=a[a['Affected protein AC']==i[:i.index('.')]]
		a_i_mut=list(a_i['Mutations'])
		mut_i=[]
		for j in range(len(a_i_mut)):
			mut_j=change_to_gly(a_i_mut[j])
			mut_i.append(mut_j)
			mut_i_1=list(set(mut_i))
		mut_i_dat=pd.DataFrame({'mutations':mut_i_1})
		mut_i_dat.to_csv(i[:i.index('.')]+'.var',index=False,header=0)

#运行provean.sh来获取结果	
import os
os.chdir('/data/liuyang/bishe_data/conservations/Reviewed/bin/Unreviewed')
import pandas as pd	
path='/data/liuyang/bishe_data/conservations/Reviewed/bin/Unreviewed'
dir_all=os.listdir(path)
for i in dir_all:
	if i[-5:] == 'fasta':
		os.system('/usr/local/bin/provean.sh -q '+i+' -v '+i[:i.index('.')]+'.var'+' >'+i[:i.index('.')]+'.out'+' --num_threads 10')
		print(dir_all.index(i))

#从out格式的输出文件中获取score信息，小于-2.5为保守残基
import os
os.chdir('/data/liuyang/bishe_data/conservations/Reviewed/bin/Unreviewed')
import pandas as pd	
path='/data/liuyang/bishe_data/conservations/Reviewed/bin/Unreviewed'
dir_all=os.listdir(path)
out_all=[x for x in dir_all if x[-3:]=='out']
scores=[]
a=pd.read_csv('Unreviewed.tsv',header=0,sep='\t')
#获取scores分值
for i in range(len(a)):
	a_p=a['Affected protein AC'][i]
	a_m=a['Mutations'][i][:-1]+'G'
	b=pd.read_csv(a_p+'.out',header=0,sep='\t')
	b1=b.T
	score_i=b1[a_m][0]
	scores.append(score_i)
scores.dat=pd.DataFrame({'provean_score':scores})
a1=pd.concat([a,scores_dat],axis=1)

#判断该位点是否保守
conservations=[]
import re
for j in range(len(a1)):
	if float(a1['provean_score'][j]) < -2.5:
		conservations.append(1)
	else:
		conservations.append(0)
conservations_dat=pd.DataFrame({'conservations':conservations})
a2=pd.concat([a1,conservations_dat],axis=1)
a2.to_csv('Unreviewed_provean_score.tsv',index=False,sep='\t')

#绘制保守性分数折线图，结果不清晰
scores_1=list(a2['provean_score'])
scores_2=[float(x) for x in scores_1]
index1=list(range(len(scores_2)))
import matplotlib.pyplot as plt
plt.plot(index1,scores_2)
plt.xlabel('Mutation sites index')
plt.ylabel('provean conservation scores')
plt.title('Provean conservation scores distribution in Unreviewed data')









#绘制柱状图
a2_con=len(a2[a2['conservations']==1])
a2_uncon=len(a2[a2['conservations']==0])
a2_y=[a2_con,a2_uncon]
a2_x=['conservative','not conservative']
color_1=['red','blue']
import matplotlib.pyplot as plt
plt.bar(a2_x,a2_y,color=color_1)
plt.xlabel('Type')
plt.ylabel('Numbers')
plt.title('Conservation type in Unreviewed data')
for m,n in zip(a2_x,a2_y):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=7)

#绘制饼图
data=[a2_con/len(a2),a2_uncon/len(a2)]
label=['conservative','not conservative']
color=['orange','green']
import matplotlib.pyplot as plt
plt.pie(data,labels=label,colors=color,shadow=True,autopct="(%1.1f%%)")
# 添加标签和标题
plt.title('Conservation type in Unreviewed data')
# 保存图片
plt.savefig('Reviewed_Disrupting.png')




'''
conda create -n consurf1 python==3.8.8
conda activate consurf1

conda deactivate
conda remove -n consurf1 --all
'''

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
#获取受影响蛋白质的结构，计算溶剂可及表面积，判断是在核心还是表面

#从alphafold-v4获取结构
a=pd.read_csv('Affected_protein_loc.tsv',header=0,sep='\t')
a1=list(a['Affected protein AC'])
a2=list(set(a1))
a2_dat=pd.DataFrame({'Affected protein AC':a2})
a2_dat.to_csv('to_pdb.txt',index=False,sep='\t')

#使用Alphafold2数据库搜索所有人类蛋白预测的pdb结构
#使用trRosetta获取其他蛋白的预测pdb结构


import requests
import json

# 读取PDB ID列表文件
a=pd.read_csv('to_pdb.txt',header=0,sep='\t')

#alphafold
#批量获取Affected protein的pdb结构

import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import PDB
import requests

import os
os.chdir('/data/liuyang/bishe_data/alphafold/Reviewed') 

a=pd.read_csv('to_pdb.txt',header=0,sep='\t')
IDs=list(a['Affected protein AC'])

headers={'User-Agent':'Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/104.0.0.0 Safari/537.36'}

no_pdb=[]
yes_pdb=[]
for i in IDs:
	url='https://alphafold.ebi.ac.uk/files/AF-'+i+'-F1-model_v4'+'.pdb'
	file_i=requests.get(url,headers=headers,verify=False)
	r=file_i.text.splitlines()
	if r[0][1] == '?':
		no_pdb.append(i)
	else:
		yes_pdb.append(i)
		with open(i+'.pdb','w') as file:
			for line in r:
				file.write(line)
				file.write('\n')


#确定无法预测结构的数据有多少
import os
os.chdir('/data/liuyang/bishe_data/colab1')
import pandas as pd
#Reviewed
a1=pd.read_csv('Reviewed.tsv',header=0,sep='\t')
no_pdb1=pd.read_csv('Reviewed_big',header=None)
no_pdb1_list=list(no_pdb1[0])
no_pdb1_data=a1[a1['Affected protein AC'].isin(no_pdb1_list)]
#最短序列长度
no_pdb1_data.index=range(len(no_pdb1_data))
length1=[]
for i in range(len(no_pdb1_data)):
	length_i=len(no_pdb1_data['A sequence'][i])
	length1.append(length_i)
min_length=sorted(length1)

#Unreviewed
a2=pd.read_csv('Unreviewed.tsv',header=0,sep='\t')
no_pdb2=pd.read_csv('Unreviewed_big',header=None)
no_pdb2_list=list(no_pdb2[0])
no_pdb2_data=a2[a2['Affected protein AC'].isin(no_pdb2_list)]
#最短序列长度
no_pdb2_data.index=range(len(no_pdb2_data))
length2=[]
for j in range(len(no_pdb2_data)):
	length_j=len(no_pdb2_data['A sequence'][j])
	length2.append(length_j)
min_length2=sorted(length2)


###############################################################################
#计算dssp
###############################################################################

##############
#Reviewed
#alphafold
import os
os.chdir('/data/liuyang/bishe_data/alphafold/Reviewed')
import pandas as pd

path1='/data/liuyang/bishe_data/alphafold/Reviewed'
path1_data=os.listdir(path1)
path1_pdb=[x for x in path1_data if x[-3:]=='pdb']

#将生成命令直接输到命令行
def get_HHblits(i):
    os.system('mkdssp -i '+i+' -o '+i[:i.index('.')]+'.dssp')
#每20个，跑一次
from joblib import Parallel,delayed
Parallel(n_jobs = 20)(delayed(get_HHblits)(i) for i in path1_pdb )

#colabfold
os.chdir('/data/liuyang/bishe_data/colab1/Reviewed')
path2='/data/liuyang/bishe_data/colab1/Reviewed'
path2_data=os.listdir(path2)
for i in path2_data:
	path2_i='/data/liuyang/bishe_data/colab1/Reviewed/'+i
	path2_i_data=os.listdir(path2_i)
	i_pdb=[x for x in path2_i_data if x[-3:]=='pdb'][0]
	os.system('mkdssp -i '+'/data/liuyang/bishe_data/colab1/Reviewed/'+i+'/'+i_pdb+' -o '+'/data/liuyang/bishe_data/colab1/Reviewed/'+i_pdb[:i_pdb.index('_')]+'.dssp')

#统计分析，绘图
'''
Location：突变在蛋白质中的位置，
COR 和 SUR 分别表示突变位于蛋白质内部和表面。
如果一个残基在蛋白质与在延展三肽中的溶液可及表面积的比值小于 0.2，
则定义该残基位于蛋白质内部，否则该残基位于蛋白质表面。
残基在蛋白质的溶液可及表面积就是DSSP算的
残基在延展三肽中的溶液可及表面积就是残基本身的表面积
'''
import os
os.chdir('/data/liuyang/bishe_data/alphafold/Reviewed')
import pandas as pd
import re
a=pd.read_csv('Reviewed.tsv',header=0,sep='\t')
path1='/data/liuyang/bishe_data/alphafold/Reviewed'
path1_data=os.listdir(path1)
alphafold_pro=[x[:x.index('.')] for x in path1_data if x[-4:]=='dssp']
path2='/data/liuyang/bishe_data/colab1/Reviewed'
path2_data=os.listdir(path2)
colabfold_pro=[y[:y.index('.')] for y in path2_data if y[-4:]=='dssp']
os.chdir('/data/liuyang/bishe_data/colab1')
#只保留有alphafold和colabfold预测结构的数据
c=pd.read_csv('Reviewed_big',header=None)
no_pdb_pro=list(c[0])
a1=a[a['Affected protein AC'].isin(alphafold_pro)]
a2=a[a['Affected protein AC'].isin(colabfold_pro)]
a_pro=pd.concat([a1,a2])
a_pro.index=range(len(a_pro))
#三肽中每种残基的溶剂可及表面积
map_surface = {'A':118.1,'R':256.0,'N':165.5,'D':158.7,'C':146.1,'Q':193.2,'E':186.2,'G':88.1,'H':202.5,'I':181.0,'L':193.1,'K':225.8,'M':203.4,'F':222.8,'P':146.8,'S':129.8,'T':152.5,'W':266.3,'Y':236.8,'V':164.5}
#函数cal_acc，从dssp计算结果文件中，获取突变位点的溶剂可及表面积与三肽溶液溶剂可及表面积的比值
#进而判断残基是核心或表面
def cal_acc(file_name_i,loc):
	start_line=0
	start_2=0
	with open(file_name_i) as f:
		for i in f:
			if i[:3].strip() !='#':
				start_line=start_line+1
			else:
				break
	with open(file_name_i) as f:
		for j in f:
			if start_2>start_line:
				if j[5:10].strip()==loc:
					res = j[13]
					res_acc = j[35:39].strip()
					res_acc_1 = float(res_acc)/map_surface[res]
			start_2=start_2+1
	return res_acc_1

dssp_acc=[]
for i in range(len(a_pro)):
	if a_pro['Affected protein AC'][i] in alphafold_pro:
		file_name_i='/data/liuyang/bishe_data/alphafold/Reviewed/'+a_pro['Affected protein AC'][i]+'.dssp'
		loc_i=int(re.findall(r'\d+',a_pro['Mutations'][i])[0])
		dssp_i=cal_acc(file_name_i,str(loc_i))
		dssp_acc.append(dssp_i)
	elif a_pro['Affected protein AC'][i] in colabfold_pro:
		file_name_i='/data/liuyang/bishe_data/colab1/Reviewed/'+a_pro['Affected protein AC'][i]+'.dssp'
		loc_i=int(re.findall(r'\d+',a_pro['Mutations'][i])[0])
		dssp_i=cal_acc(file_name_i,str(loc_i))
		dssp_acc.append(dssp_i)

acc_dat=pd.DataFrame({'dssp_acc':dssp_acc})
a_dssp=pd.concat([a_pro,acc_dat],axis=1)
#添加类型
residue_location=[]
for p in range(len(a_dssp)):
	if a_dssp['dssp_acc'][p]<0.2:
		residue_location.append('Core')
	else:
		residue_location.append('Surface')
residue_dat=pd.DataFrame({'residue_location':residue_location})
a_result=pd.concat([a_dssp,residue_dat],axis=1)
os.chdir('/data/liuyang/bishe_data/alphafold')
a_result.to_csv('Reviewed_dssp.tsv',index=False,sep='\t')

#绘制统计图
#绘制acc分数折线图，结果不清晰
import os
os.chdir('/data/liuyang/bishe_data/alphafold')
import pandas as pd
a2=pd.read_csv('Reviewed_dssp.tsv',header=0,sep='\t')
scores_1=list(a2['dssp_acc'])
scores_2=[float(x) for x in scores_1]
index1=list(range(len(scores_2)))
import matplotlib.pyplot as plt
plt.plot(index1,scores_2)
plt.xlabel('Mutation sites index')
plt.ylabel('dssp acc scores')
plt.title('dssp acc distribution in Reviewed data')

#绘制柱状图
a2_core=len(a2[a2['residue_location']=='Core'])
a2_surface=len(a2[a2['residue_location']=='Surface'])
a2_y=[a2_core,a2_surface]
a2_x=['Core','Surface']
color_1=['red','blue']
import matplotlib.pyplot as plt
plt.bar(a2_x,a2_y,color=color_1)
plt.xlabel('Type')
plt.ylabel('Numbers')
plt.title('Residue location in Unreviewed data')
for m,n in zip(a2_x,a2_y):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=7)

#绘制饼图
data=[a2_core/len(a2),a2_surface/len(a2)]
label=['Core','surface']
color=['orange','green']
import matplotlib.pyplot as plt
plt.pie(data,labels=label,colors=color,shadow=True,autopct="(%1.1f%%)")
# 添加标签和标题
plt.title('Residue location in Unreviewed data')




##############
#Unreviewed
#alphafold
import os
os.chdir('/data/liuyang/bishe_data/alphafold/Unreviewed')
import pandas as pd

path1='/data/liuyang/bishe_data/alphafold/Unreviewed'
path1_data=os.listdir(path1)
path1_pdb=[x for x in path1_data if x[-3:]=='pdb']
#将生成命令直接输到命令行
def get_HHblits(i):
    os.system('mkdssp -i '+i+' -o '+i[:i.index('.')]+'.dssp')
#每20个，跑一次
from joblib import Parallel,delayed
Parallel(n_jobs = 20)(delayed(get_HHblits)(i) for i in path1_pdb )

#colabfold
os.chdir('/data/liuyang/bishe_data/colab1/Unreviewed')
path2='/data/liuyang/bishe_data/colab1/Unreviewed'
path2_data=os.listdir(path2)
for i in path2_data:
	path2_i='/data/liuyang/bishe_data/colab1/Unreviewed/'+i
	path2_i_data=os.listdir(path2_i)
	i_pdb=[x for x in path2_i_data if x[-3:]=='pdb'][0]
	os.system('mkdssp -i '+'/data/liuyang/bishe_data/colab1/Unreviewed/'+i+'/'+i_pdb+' -o '+'/data/liuyang/bishe_data/colab1/Unreviewed/'+i_pdb[:i_pdb.index('_')]+'.dssp')


#统计分析
import os
os.chdir('/data/liuyang/bishe_data/alphafold/Unreviewed')
import pandas as pd
import re
a=pd.read_csv('Unreviewed.tsv',header=0,sep='\t')
path1='/data/liuyang/bishe_data/alphafold/Unreviewed'
path1_data=os.listdir(path1)
alphafold_pro=[x[:x.index('.')] for x in path1_data if x[-4:]=='dssp']
path2='/data/liuyang/bishe_data/colab1/Unreviewed'
path2_data=os.listdir(path2)
colabfold_pro=[y[:y.index('.')] for y in path2_data if y[-4:]=='dssp']
os.chdir('/data/liuyang/bishe_data/colab1')
#只保留有alphafold和colabfold预测结构的数据
c=pd.read_csv('Unreviewed_big',header=None)
no_pdb_pro=list(c[0])
a1=a[a['Affected protein AC'].isin(alphafold_pro)]
a2=a[a['Affected protein AC'].isin(colabfold_pro)]
a_pro=pd.concat([a1,a2])
a_pro.index=range(len(a_pro))
#三肽中每种残基的溶剂可及表面积
map_surface = {'A':118.1,'R':256.0,'N':165.5,'D':158.7,'C':146.1,'Q':193.2,'E':186.2,'G':88.1,'H':202.5,'I':181.0,'L':193.1,'K':225.8,'M':203.4,'F':222.8,'P':146.8,'S':129.8,'T':152.5,'W':266.3,'Y':236.8,'V':164.5}
#函数cal_acc，从dssp计算结果文件中，获取突变位点的溶剂可及表面积与三肽溶液溶剂可及表面积的比值
#进而判断残基是核心或表面
def cal_acc(file_name_i,loc):
	start_line=0
	start_2=0
	with open(file_name_i) as f:
		for i in f:
			if i[:3].strip() !='#':
				start_line=start_line+1
			else:
				break
	with open(file_name_i) as f:
		for j in f:
			if start_2>start_line:
				if j[5:10].strip()==loc:
					res = j[13]
					res_acc = j[35:39].strip()
					res_acc_1 = float(res_acc)/map_surface[res]
			start_2=start_2+1
	return res_acc_1

dssp_acc=[]
for i in range(len(a_pro)):
	if a_pro['Affected protein AC'][i] in alphafold_pro:
		file_name_i='/data/liuyang/bishe_data/alphafold/Unreviewed/'+a_pro['Affected protein AC'][i]+'.dssp'
		loc_i=int(re.findall(r'\d+',a_pro['Mutations'][i])[0])
		dssp_i=cal_acc(file_name_i,str(loc_i))
		dssp_acc.append(dssp_i)
	elif a_pro['Affected protein AC'][i] in colabfold_pro:
		file_name_i='/data/liuyang/bishe_data/colab1/Unreviewed/'+a_pro['Affected protein AC'][i]+'.dssp'
		loc_i=int(re.findall(r'\d+',a_pro['Mutations'][i])[0])
		dssp_i=cal_acc(file_name_i,str(loc_i))
		dssp_acc.append(dssp_i)

acc_dat=pd.DataFrame({'dssp_acc':dssp_acc})
a_dssp=pd.concat([a_pro,acc_dat],axis=1)
#添加类型
residue_location=[]
for p in range(len(a_dssp)):
	if a_dssp['dssp_acc'][p]<0.2:
		residue_location.append('Core')
	else:
		residue_location.append('Surface')
residue_dat=pd.DataFrame({'residue_location':residue_location})
a_result=pd.concat([a_dssp,residue_dat],axis=1)
os.chdir('/data/liuyang/bishe_data/alphafold')
a_result.to_csv('Unreviewed_dssp.tsv',index=False,sep='\t')

#绘制统计图
#绘制acc分数折线图，结果不清晰
import os
os.chdir('/data/liuyang/bishe_data/alphafold')
import pandas as pd
a2=pd.read_csv('Unreviewed_dssp.tsv',header=0,sep='\t')
scores_1=list(a2['dssp_acc'])
scores_2=[float(x) for x in scores_1]
index1=list(range(len(scores_2)))
import matplotlib.pyplot as plt
plt.plot(index1,scores_2)
plt.xlabel('Mutation sites index')
plt.ylabel('dssp acc scores')
plt.title('dssp acc distribution in Reviewed data')

#绘制柱状图
a2_core=len(a2[a2['residue_location']=='Core'])
a2_surface=len(a2[a2['residue_location']=='Surface'])
a2_y=[a2_core,a2_surface]
a2_x=['Core','Surface']
color_1=['red','blue']
import matplotlib.pyplot as plt
plt.bar(a2_x,a2_y,color=color_1)
plt.xlabel('Type')
plt.ylabel('Numbers')
plt.title('Residue location in Unreviewed data')
for m,n in zip(a2_x,a2_y):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=7)

#绘制饼图
data=[a2_core/len(a2),a2_surface/len(a2)]
label=['Core','surface']
color=['orange','green']
import matplotlib.pyplot as plt
plt.pie(data,labels=label,colors=color,shadow=True,autopct="(%1.1f%%)")
# 添加标签和标题
plt.title('Residue location in Unreviewed data')

'''
#去除冲突数据
import os
os.chdir('/data/liuyang/bishe_data/alphafold/Unreviewed') 
import pandas as pd
a=pd.read_csv('Unreviewed.tsv',header=0,sep='\t')
a_pro=list(set(list(a['Affected protein AC'])))
#alphafold结构
path1='/data/liuyang/bishe_data/alphafold/Unreviewed'
list_data1=os.listdir(path1)
path1_pro=[x[:-4] for x in list_data1 if x[-3:]=='pdb']
rm1_pro=[y for y in path1_pro if y not in a_pro]
#colabfold结构
path2='/data/liuyang/bishe_data/colab1/Unreviewed'
list_data2=os.listdir(path2)
path2_pro=[m[:m.index('_')] for m in list_data2]
rm2_pro=[n for n in path2_pro if n not in a_pro]
#无法预测结构的数据
os.chdir('/data/liuyang/bishe_data/colab1')
no_pdb=pd.read_csv('Unreviewed_big',header=None,sep='\t')
no_pdb_pro=list(no_pdb[0])
rm3_pro=[p for p in no_pdb_pro if p not in a_pro]
'''

##################################################################################################
##################################################################################################
##################################################################################################
#Unreviewed

import os
os.chdir('/data/liuyang/bishe_data/alphafold/Unreviewed')
import pandas as pd
import re
a=pd.read_csv('Unreviewed_completed.tsv',header=0,sep='\t')
#按Affected protein进行分组
a1=list(a.groupby('Affected protein AC'))
dat1=pd.DataFrame()
for i in range(len(a1)):
	b=a1[i][1]
	b.index=range(len(b))
	for j in range(len(b)):
		b1=int(re.findall(r'\d+',b['Mutations'][j])[0])
		b2=b['Mutations'][j][0]
		b3=b['Affected protein AC'][j]
		b4=b['Interaction AC'][j]
		b5=b['b id'][j]
		b_dic={'Affected protein AC':[b3],'location':[b1],'Original sequence':[b2],'Interaction AC':[b4],'B_ID':[b5]}
		b_dat=pd.DataFrame(b_dic)
		dat1=pd.concat([dat1,b_dat])
dat1.to_csv('Affected_protein_loc.tsv',index=False,sep='\t')



#获取序列
b=pd.read_csv('Affected_protein_loc.tsv',header=0,sep='\t')
b1=list(b['Affected protein AC'])
b2=list(set(b1))
b3=pd.DataFrame({'Affected protein AC':b2})
b3.to_csv('protein_id.txt',index=False,sep='\t')




'''
conda create -n consurf1 python==3.8.8
conda activate consurf1

conda deactivate
conda remove -n consurf1 --all
'''

######################################################################################################################################################
#从alphafold-v4获取结构
a=pd.read_csv('Affected_protein_loc.tsv',header=0,sep='\t')
a1=list(a['Affected protein AC'])
a2=list(set(a1))
a2_dat=pd.DataFrame({'Affected protein AC':a2})
a2_dat.to_csv('to_pdb.txt',index=False,sep='\t')

#使用Alphafold2数据库搜索所有人类蛋白预测的pdb结构
#使用trRosetta获取其他蛋白的预测pdb结构


import requests
import json
import os
os.chdir('/data/liuyang/bishe_data/alphafold/Unreviewed') 
#alphafold
#批量获取Affected protein的pdb结构
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import PDB
import requests

# 读取PDB ID列表文件
a=pd.read_csv('to_pdb.txt',header=0,sep='\t')
IDs=list(a['Affected protein AC'])

headers={'User-Agent':'Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/104.0.0.0 Safari/537.36'}

no_pdb=[]
yes_pdb=[]
for i in IDs:
	url='https://alphafold.ebi.ac.uk/files/AF-'+i+'-F1-model_v4'+'.pdb'
	file_i=requests.get(url,headers=headers,verify=False)
	r=file_i.text.splitlines()
	if r[0][1] == '?':
		no_pdb.append(i)
	else:
		yes_pdb.append(i)
		with open(i+'.pdb','w') as file:
			for line in r:
				file.write(line)
				file.write('\n')





'''

######################################################################################################################################################
#获取trRosetta预测结构
#使用trRosetta环境

#Reviewed
import os
os.chdir('/data/liuyang/bishe_data/trRosetta/step/Reviewed')
import pandas as pd


from Bio import SeqIO
# 打开输入fasta文件
input_file = "Reviewed.fasta"

# 读取输入fasta文件中的序列，并将每个序列写入单独的文件
for record in SeqIO.parse(input_file, "fasta"):
    output_file = "{}.fasta".format(record.id[record.id.index('|')+1:record.id.index('|',record.id.index('|')+1)])
    with open(output_file, "w") as output_handle:
        SeqIO.write(record, output_handle, "fasta")
print("拆分完毕！")


#使用HH-suite3
hhblits -i Q9RA63.fasta -d pfam -oa3m Q9RA63.a3m

#批量处理
import os
os.chdir('/data/liuyang/bishe_data/trRosetta/step/Reviewed')

path='/data/liuyang/bishe_data/trRosetta/step/Reviewed'
a=os.listdir(path)
x=[]
for i in a:
	if i[-5:] !='fasta':
		x.append(i)
y=[j for j in a if j not in x]

#将生成命令直接输到命令行
def get_HHblits(idd):
    os.system('hhblits -i '+idd+' -d pfam -oa3m '+idd[:idd.index('.')]+'.a3m')

#每10个，跑一次
from joblib import Parallel,delayed
Parallel(n_jobs = 10)(delayed(get_HHblits)(idd) for idd in y )

#结果获得a3m格式的多序列比对文件319个



#安装trRosetta
# download package
git clone https://github.com/gjoni/trRosetta
cd trRosetta
# download pre-trained network
wget https://files.ipd.uw.edu/pub/trRosetta/model2019_07.tar.bz2
tar xf model2019_07.tar.bz2


#计算距离和方向
python /data/liuyang/bishe_data/pfamA/trRosetta/network/predict.py  -m /data/liuyang/bishe_data/pfamA/trRosetta/model2019_07 P27709.a3m P27709.npz

python /data/liuyang/bishe_data/pfamA/trRosetta/network/predict.py -m /data/liuyang/bishe_data/pfamA/trRosetta/model2019_07 P27709.a3m P27709.npz
#批量处理
import os
os.chdir('/data/liuyang/bishe_data/trRosetta/step/Reviewed/trRosetta')
import pandas as pd
path='/data/liuyang/bishe_data/trRosetta/step/Reviewed/trRosetta'
a=os.listdir(path)
x=[]
for i in a:
	if i[-3:] !='a3m':
		x.append(i)
y=[j for j in a if j not in x]
z=[]
for p in a:
	if p[-3:]=='npz':
		z.append(p[:p.index('.')])
y1=[q for q in y if q[:q.index('.')] not in z]
#使用sh批处理
import pandas as pd
y_codes=[]
for p in y:
	y_code='python /data/liuyang/bishe_data/trRosetta/step/Reviewed/trRosetta/network/predict.py  -m /data/liuyang/bishe_data/trRosetta/step/Reviewed/trRosetta/model2019_07 '+p+' '+p[:p.index('.')]+'.npz'
	y_codes.append(y_code)
y_dat=pd.DataFrame({'y_codes':y_codes})
y_dat.to_csv('trRoestta_code1.txt',index=False,sep='\t')


#修改成批处理文件，最前面加#!/bin/bash
#source trRoestta_code.sh






#####################################################################################################
#####################################################################################################
#####################################################################################################

#Unreviewed

import os
os.chdir('/data/liuyang/bishe_data/trRosetta/step/Unreviewed')
import pandas as pd


from Bio import SeqIO
# 打开输入fasta文件
input_file = "Unreviewed.fasta"

# 读取输入fasta文件中的序列，并将每个序列写入单独的文件
for record in SeqIO.parse(input_file, "fasta"):
    output_file = "{}.fasta".format(record.id[record.id.index('|')+1:record.id.index('|',record.id.index('|')+1)])
    with open(output_file, "w") as output_handle:
        SeqIO.write(record, output_handle, "fasta")
print("拆分完毕！")


#使用HH-suite3
hhblits -i Q9RA63.fasta -d pfam -oa3m Q9RA63.a3m

#批量处理
import os
os.chdir('/data/liuyang/bishe_data/trRosetta/step/Unreviewed')

path='/data/liuyang/bishe_data/trRosetta/step/Unreviewed'
a=os.listdir(path)
x=[]
for i in a:
	if i[-5:] !='fasta':
		x.append(i)
y=[j for j in a if j not in x]

#将生成命令直接输到命令行
def get_HHblits(idd):
    os.system('hhblits -i '+idd+' -d pfam -oa3m '+idd[:idd.index('.')]+'.a3m')
#每10个，跑一次
from joblib import Parallel,delayed
Parallel(n_jobs = 10)(delayed(get_HHblits)(idd) for idd in y )


#计算方向和距离
import os
os.chdir('/data/liuyang/bishe_data/trRosetta/step/Unreviewed')
import pandas as pd
path='/data/liuyang/bishe_data/trRosetta/step/Unreviewed'
a=os.listdir(path)
x=[]
for i in a:
	if i[-3:] !='a3m':
		x.append(i)
y=[j for j in a if j not in x]

#使用sh批处理
import pandas as pd
y_codes=[]
for p in y:
	y_code='python /data/liuyang/bishe_data/trRosetta/step/Unreviewed/trRosetta/network/predict.py  -m /data/liuyang/bishe_data/trRosetta/step/Unreviewed/trRosetta/model2019_07 '+p+' '+p[:p.index('.')]+'.npz'
	y_codes.append(y_code)
y_dat=pd.DataFrame({'y_codes':y_codes})
y_dat.to_csv('trRoestta_code.txt',index=False,sep='\t')

'''



############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################

#预测binding site
#将所有蛋白质的pdb都转移到py文件的路径下
#Reviewed
import os
os.chdir('/data/liuyang/bishe_data/scannet/ScanNet-main')
import pandas as pd
path1='/data/liuyang/bishe_data/scannet/ScanNet-main'
path1_list=os.listdir(path1)
path1_pdb=[x for x in path1_list if x[-3:]=='pdb']
codes=[]
for i in path1_pdb:
	code_i='python predict_bindingsites.py '+i+' --noMSA'
	codes.append(code_i)
codes_dat=pd.DataFrame({'codes_i':codes})
codes_dat.to_csv('Reviewed_scannet.txt',index=False,sep='\t')
#将txt文件，修改成sh批处理文件运行
'''
def get_HHblits(idd):
    os.system('python predict_bindingsites.py '+idd+' --noMSA')
#每1个，跑一次
from joblib import Parallel,delayed
Parallel(n_jobs = 1)(delayed(get_HHblits)(idd) for idd in path1_pdb)
'''

import os
os.chdir('/data/liuyang/bishe_data/scannet/ScanNet-main')
import pandas as pd
path1='/data/liuyang/bishe_data/colab1/Reviewed'
path1_data=os.listdir(path1)
pro_path1=[x for x in path1_data if x[-4:] != 'dssp']
codes1=[]
for i in pro_path1:
	path2='/data/liuyang/bishe_data/colab1/Reviewed'+i
	path2_data=os.listdir(path2)
	pro_pdb=[y for y in path2_data if y[-3:]=='pdb'][0]
	codes1.append('python predict_bindingsites.py /data/liuyang/bishe_data/colab1/Reviewed/'+i+'/'+pro_pdb+' --noMSA')
codes1_dat=pd.DataFrame({'codes_i':codes1})
codes1_dat.to_csv('Reviewed_colabfold.txt',index=False,sep='\t')

python predict_bindingsites.py /data/liuyang/bishe_data/colab1/Reviewed/B3FK34_9b172/B3FK34_9b172_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb --noMSA

#统计学分析
import os
os.chdir('/data/liuyang/bishe_data/scannet/ScanNet-main')
import pandas as pd
import re
a=pd.read_csv('Reviewed1.tsv',header=0,sep='\t')
path1='/data/liuyang/bishe_data/scannet/ScanNet-main/predictions'
path1_data=os.listdir(path1)
probability_data=[]
for i in range(len(a)):
	loc1=int(re.findall(r'\d+',a['Mutations'][i])[0])
	for j in path1_data:
		if j[:j.index('_')]==a['Affected protein AC'][i]:
			path2='/data/liuyang/bishe_data/scannet/ScanNet-main/predictions/'+j
			path2_data=os.listdir(path2)
			file_i=[x for x in path2_data if x[-3:]=='csv'][0]
			file=pd.read_csv(path2+'/'+file_i,header=0,sep=',')
			probability_i=list(file[file['Residue Index']==loc1]['Binding site probability'])[0]
			probability_data.append(probability_i)



###############################################################################################################
#Unreviewed

import os
os.chdir('/data/liuyang/bishe_data/scannet/ScanNet-main')
import pandas as pd
path1='/data/liuyang/bishe_data/colab1/Unreviewed'
path1_data=os.listdir(path1)
pro_path1=[x for x in path1_data if x[-4:] != 'dssp']
codes1=[]
for i in pro_path1:
	path2='/data/liuyang/bishe_data/colab1/Unreviewed'+i
	path2_data=os.listdir(path2)
	pro_pdb=[y for y in path2_data if y[-3:]=='pdb'][0]
	codes1.append('python predict_bindingsites.py /data/liuyang/bishe_data/colab1/Unreviewed/'+i+'/'+pro_pdb+' --noMSA')
codes1_dat=pd.DataFrame({'#!/bin/bash':codes1})
codes1_dat.to_csv('Unreviewed_colabfold.sh',index=False)

path2='/data/liuyang/bishe_data/alphafold/Unreviewed'
path2_data=os.listdir(path2)
pro_path2=[x for x in path2_data if x[-3:]=='pdb']
codes2=[]
for j in pro_path2:
	codes2.append('python predict_bindingsites.py /data/liuyang/bishe_data/alphafold/Unreviewed/'+j+' --noMSA')
codes2_dat=pd.DataFrame({'#!/bin/bash':codes2})
codes2.to_csv('Unreviewed_alphafold.sh',index=False)






###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
#获取结构信息
#Reviewed
#使用pdbe的restful接口，获取结构
import os
os.chdir('/data/liuyang/bishe_data/get_pdb/Reviewed')
import pandas as pd

import requests
import json
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import PDB

a=pd.read_csv('Reviewed_completed.tsv',header=0,sep='\t')
a1=list(a['Affected protein AC'])
IDs=list(set(a1))

def get_message(x):
	for k in list(x.keys()):
		x[k]=[x[k]]
	x_dat=pd.DataFrame(x)
	return x_dat

headers={'User-Agent':'Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/104.0.0.0 Safari/537.36'}
no_pdb=[]
yes_pdb=[]
result1=pd.DataFrame()
for i in IDs:
	url='https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures/'+i
	file_i=requests.get(url,headers=headers,verify=False)
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

result1.to_csv('best_structures.tsv',index=False,sep='\t')





import os
os.chdir('/data/liuyang/bishe_data/get_pdb/Reviewed')
import requests
import json
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import PDB
headers={'User-Agent':'Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/104.0.0.0 Safari/537.36'}

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


a1=pd.read_csv('best_structure.tsv',header=0,sep='\t')
a=a1.drop_duplicates(subset=['Uniprot_id','end','chain_id','pdb_id','start','unp_end','coverage','unp_start','resolution','experimental_method','tax_id'],keep='first')
a.index=range(len(a))
no_pro=[]
duo=[]
data1=pd.DataFrame()
for j in range(len(a)):
	a_j=a['pdb_id'][j]
	url='https://www.ebi.ac.uk/pdbe/api/mappings/all_isoforms/'+a_j
	file_j=requests.get(url,headers=headers,verify=False)
	r=json.loads(file_j.content.decode())
	if r=={}:
		no_pro.append(a_j)
	else:
		b=get_message(r)
		data1=pd.concat([data1,b])
data1.to_csv('Reviewed_SIFTS_mapping.tsv',index=False,sep='\t')

#获取所需结构
import os
os.chdir('/data/liuyang/bishe_data/get_pdb/Reviewed/map')
import pandas as pd
import math
a=pd.read_csv('Reviewed.tsv',header=0,sep='\t')
b=pd.read_csv('result1.tsv',header=0,sep='\t')
b1=pd.concat([b['uniprot_id'],b['pdb_id'],b['coverage'],b['resolution'],b['experimental_method']],axis=1)
b2=b1.drop_duplicates(subset=['uniprot_id','pdb_id','coverage','resolution','experimental_method'],keep='first')

c=pd.read_csv('Reviewed_SIFTS_mapping.tsv',header=0,sep='\t')
data_c=[]
entity_groups=[]
for i in range(len(c)):
	entity_group=[]
	for j in range(435):
		entity_i=c['entity_id_'+str(j)][i]
		entity_group.append(entity_i)
	entity_group1=[x for x in entity_group if math.isnan(x) == False]
	entity_group2=set(entity_group1)
	entity_groups.append(entity_group2)
entity_i_dat=pd.DataFrame({'entity_group':entity_groups})
c_n=pd.concat([c,entity_i_dat],axis=1)
#判断是否带有其他元件
type1=[]
for j in range(len(c_n)):
	if len(c_n['entity_group'][j])>2: #将大于二聚体的标记为1
		type1.append(1)
	elif c_n['entity_group'][j]=={1} or c_n['entity_group'][j]=={1,2}: #将正确二聚体标记为0
		type1.append(0)
	else: #将含有DNA等配体的标记为2
		type1.append(2)
type_dat=pd.DataFrame({'ligand':type1})
c_m=pd.concat([c_n,type_dat],axis=1)
c_m.to_csv('pdb_ligand.tsv',index=False,sep='\t')



#将pdb数据与突变数据整合起来
a=pd.read_csv('Reviewed.tsv',header=0,sep='\t')
b=pd.read_csv('pdb_ligand.tsv',header=0,sep='\t')
#将a，b中的相互作用蛋白都放在集合中，保存在一个新列中
#a
proteins=[]
for i in range(len(a)):
	pro_group=[]
	a_i_1=a['Affected protein AC'][i]
	a_i_2=a['b id'][i]
	pro_group.append(a_i_1)
	pro_group.append(a_i_2)
	pro_set=set(pro_group)
	proteins.append(pro_set)
protein_dat=pd.DataFrame({'protein_group':proteins})
a1=pd.concat([a,protein_dat],axis=1)
#b,b中最多有435种蛋白质的列
proteins1=[]
for j in range(len(b)):
	pro_group=[]
	for p in range(435):
		b_p=b['Uniprot_id_'+str(p)][j]
		pro_group.append(b_p)
	pro_group1=[x for x in pro_group if pd.isna(x) == False]
	pro_set=set(pro_group1)
	proteins1.append(pro_set)
protein_dat1=pd.DataFrame({'protein_group':proteins1})
b1=pd.concat([b,protein_dat1],axis=1)

a1.to_csv('Reviewed_group.tsv',index=False,sep='\t')
b1.to_csv('Reviewed_SIFTS_group.tsv',index=False,sep='\t')



#pdb复合物中，受突变影响蛋白质Affected protein AC不应该为多条链上都有？
#因此做标记
a1=pd.read_csv('Reviewed_group.tsv',header=0,sep='\t')
b1=pd.read_csv('Reviewed_SIFTS_group.tsv',header=0,sep='\t')
protein_str=[]
for m in range(len(a1)):
	temp_m=sorted(list(a1['protein_group'][m]))
	temp_n=''
	for n in len(temp_m):
		temp_n=
for i in range(len(b1)):
	for j in b1['']




#首先确定每个复合物中相同的蛋白，uniprot起始终止是不是一样
data_all1=pd.DataFrame()
for m in range(len(a1)):
	for n in range(len(b1)):
		if eval(a1['protein_group'][m]).issubset(eval(b1['protein_group'][n])) == True:
			x=a1[m:m+1]
			x.index=range(len(x))
			y=b1[n:n+1]
			y.index=range(len(y))
			data_m=pd.concat([x,y],axis=1)
			data_all1=pd.concat([data_all1,data_m])
	print(m)
#太慢
list1=[]
for p in range(len(a1)):
	list1.append(tuple(a1['protein_group'][p]))
list1_dat=pd.DataFrame({'list_group':list1})
a2=pd.concat([a1,list1_dat],axis=1)
list2=[]
for q in range(len(b1)):
	list2.append(tuple(b1['protein_group'][q]))
list2_dat=pd.DataFrame({'list_group':list2})
b2=pd.concat([b1,list2_dat],axis=1)

a3=list(a2.groupby('list_group'))
b3=list(b2.groupby('list_group'))

data_i=pd.DataFrame()
for i in range(len(a3)):
	for j in range(len(b3)):
		if set(a3[i][0]).issubset(b3[j][0]):
			temp_b=b3[j][1]
			temp_b.index=range(len(temp_b))
			temp_a=a3[i][1]
			temp_a.index=range(len(temp_a))
			for m in range(len(temp_a)):
				temp_a_m=temp_a[m:m+1]
				temp_a_ms=pd.concat([temp_a_m]*len(temp_b))
				temp_a_ms.index=range(len(temp_a_ms))
				temp_data_i=pd.concat([temp_a_ms,temp_b],axis=1)
				data_i=pd.concat([data_i,temp_data_i])
	print(i)


#合并
data_all=pd.DataFrame()
is_duo=[]
import re
for m in range(len(a1)):
	for n in range(len(b1)):
		if a1['protein_group'][m].issubset(b1['protein_group'][n]) == True:
			num1=0
			nums=[]
			for p in range(435):
				if b1['Uniprot_id_'+str(p)][n]==a1['Affected protein AC'][m]:
					num1=num1+1
					nums.append(p)
			#确定Affected protein在复合物中只存在于一条链上
			if num1 == 1:
				start1=b1['unp_start_'+str(nums[0])][n]
				end1=b1['unp_end_'+str(nums[0])][n]
				loc1=int(re.findall(r'\d+',a1['Mutations'][m])[0])
				if loc1 >= start1 and loc1 <= end1:
					a_m=a1[m:m+1]
					a_m.index=range(len(a_m))
					b_n=b1[n:n+1]
					b_n.index=range(len(b_n))
					data_m=pd.concat([a_m,b_n],axis=1)
					data_all=pd.concat([data_all,data_m])
	print(m)
data_all.to_csv('Reviewed_map.tsv',index=False,sep='\t')










#将突变位点从uniprot序列，map到pdb序列上
import os
os.chdir('/data/liuyang/bishe_data/get_pdb/Reviewed/map')
import pandas as pd
a=pd.read_csv('Reviewed_map.tsv',header=0,sep='\t')
for i in range(len(a)):
	for j in range(435):
		if a['uniprot_id_'+str(j)][i]==a['Affected protein AC'][i]:
			mut_chain=a['chain_id_'+str(j)][i]





'''
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
				m_i_dic={'Uniprot_id_'+index1:i,'identifier_'+index1:[m_i['identifier']],'chain_id_'+index1:[m_i['mappings'][p]['chain_id']],'unp_end_'+index1:[m_i['mappings'][p]['unp_end']],
				'unp_start_'+index1:[m_i['mappings'][p]['unp_start']],'pdb_start_'+index1:[m_i['mappings'][p]['pdb_start']],'pdb_end_'+index1:[m_i['mappings'][p]['pdb_end']],
				'is_canonical_'+index1:[m_i['mappings'][p]['is_canonical']]}
				m_i_dat=pd.DataFrame(m_i_dic)
				dat1=pd.concat([dat1,m_i_dat],axis=1)
				index2=index2+1
		else:
			index1=str(index2)
			m_i_dic={'Uniprot_id_'+index1:i,'identifier_'+index1:[m_i['identifier']],'chain_id_'+index1:[m_i['mappings'][0]['chain_id']],'unp_end_'+index1:[m_i['mappings'][0]['unp_end']],
			'unp_start_'+index1:[m_i['mappings'][0]['unp_start']],'pdb_start_'+index1:[m_i['mappings'][0]['pdb_start']],'pdb_end_'+index1:[m_i['mappings'][0]['pdb_end']],
			'is_canonical_'+index1:[m_i['mappings'][0]['is_canonical']]}
			m_i_dat=pd.DataFrame(m_i_dic)
			dat1=pd.concat([dat1,m_i_dat],axis=1)
			index2=index2+1
	dat2=pd.concat([pdb_dat,dat1],axis=1)
	if index1==1:
		return dat2
	else:
		return None


a1=pd.read_csv('best_structure.tsv',header=0,sep='\t')
a=a1.drop_duplicates(subset=['Uniprot_id','end','chain_id','pdb_id','start','unp_end','coverage','unp_start','resolution','experimental_method','tax_id'],keep='first')
a.index=range(len(a))
no_pro=[]
duo=[]
data1=pd.DataFrame()
for j in range(len(a)):
	a_j=a['pdb_id'][j]
	url='https://www.ebi.ac.uk/pdbe/api/mappings/all_isoforms/'+a_j
	file_j=requests.get(url,headers=headers,verify=False)
	r=json.loads(file_j.content.decode())
	if r=={}:
		no_pro.append(a_j)
	else:
		b=get_message(r)
		data1=pd.concat([data1,b])
'''




#################################################################################################################################################
#Unreviewed
#使用pdbe的restful接口，获取结构
import os
os.chdir('/data/liuyang/bishe_data/get_pdb/Unreviewed')
import pandas as pd

import requests
import json
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import PDB

a=pd.read_csv('Unreviewed_completed.tsv',header=0,sep='\t')
a1=list(a['Affected protein AC'])
IDs=list(set(a1))

def get_message(x):
	for k in list(x.keys()):
		x[k]=[x[k]]
	x_dat=pd.DataFrame(x)
	return x_dat

headers={'User-Agent':'Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/104.0.0.0 Safari/537.36'}
no_pdb=[]
yes_pdb=[]
result1=pd.DataFrame()
for i in IDs:
	url='https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures/'+i
	file_i=requests.get(url,headers=headers,verify=False)
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

result1.to_csv('best_structures.tsv',index=False,sep='\t')



#筛选：去除实验方法不是cryo-EM、NMR、X-RAY的数据-------------------去除分辨率大于3的pdb结构
























sha1:6a6933be9da3:90eb47ad952c305ecac037906d6f788c2096c596

