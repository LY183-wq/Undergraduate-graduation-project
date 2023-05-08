# -*- coding: utf-8 -*-
# @Author: Liu Yang
# @Date:   2023-03-16 14:40:25
# @Last Modified by:   Liu Yang
# @Last Modified time: 2023-03-16 15:31:51

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################


#绘制饼状图
import os
os.chdir('/data/liuyang/bishe_data/data_tongji')
import pandas as pd

#####################################################################################################################################################
#Reviewed
a=pd.read_csv('Reviewed_right.tsv',header=0,sep='\t')
a1=a[a['type']=='Disrupting']
a2=a[a['type']=='Decreasing']
a3=a[a['type']=='Increasing']
a4=a[a['type']=='No effect']
data=[len(a1)/len(a),len(a2)/len(a),len(a3)/len(a),len(a4)/len(a)]
label=['Disrupting','Decreasing','Increasing','No effect']
color=['red','orange','blue','green']
import matplotlib.pyplot as plt
plt.pie(data,labels=label,colors=color,shadow=True,autopct="(%1.1f%%)")
# 添加标签和标题
# 保存图片
plt.savefig('Reviewed_pie.png')








import os
os.chdir('/data/liuyang/PNG')
import pandas as pd
import math

#####################################################################################################################################################
#Reviewed
a=pd.read_csv('Reviewed_right.tsv',header=0,sep='\t')
a1=a[a['type']=='Disrupting']
a2=a[a['type']=='Decreasing']
a3=a[a['type']=='Increasing']
a4=a[a['type']=='No effect']

#####################################################################################
#Disrupting
index1=[]
index_1=[]
index_a1=[]
a1_len=[]
a1_list=list(a1.groupby('Interaction AC'))
for i in range(len(a1_list)):
	len_num=len(a1_list[i][1])
	index1.append(len_num)
	index_1=index_1+list(a1_list[i][1].index)
	a1_len.append(a1_list[i][0])
data_int=pd.DataFrame({'Interaction AC':a1_len})
data_num=pd.DataFrame({'numbers':index1})
data1=pd.concat([data_int,data_num],axis=1)

nums=list(set(index1))
mut_nums=[]
nums_j=[]
for j in nums:
	num_j=len(data1[data1['numbers']==j])
	nums_j.append(num_j)
	mut_num=j
	mut_nums.append(mut_num)
mut_dat=pd.DataFrame({'mut_numbers':mut_nums})
int_dat=pd.DataFrame({'Interaction numbers':nums_j})
x=pd.concat([mut_dat,int_dat],axis=1)

#由于数据差距太原，因此将横坐标取log
mut_nums1=[]
for p in range(len(mut_nums)):
	temp_p=math.log(mut_nums[p],2)
	mut_nums1.append(temp_p)
#绘制统计柱状图
import matplotlib.pyplot as plt
plt.figure(figsize=(15,12))
plt.bar(mut_nums1,nums_j,width=0.5)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
# 添加标签和标题
plt.xlabel('Log2(Mutation numbers)',fontsize=40)
plt.ylabel('Interaction numbers',fontsize=40)
# 保存图片
plt.savefig('Reviewed_Disrupting.png',dpi=500)


#绘制折线图
import matplotlib.pyplot as plt
plt.plot(mut_nums,nums_j)
plt.xlabel('Mutation numbers')
plt.ylabel('Interaction numbers')
plt.title('Disrupting mutation numbers in Reviewed data')




#####################################################################################
#Decreasing
index1=[]
index_1=[]
index_a1=[]
a1_len=[]
a1_list=list(a2.groupby('Interaction AC'))
for i in range(len(a1_list)):
	len_num=len(a1_list[i][1])
	index1.append(len_num)
	index_1=index_1+list(a1_list[i][1].index)
	a1_len.append(a1_list[i][0])
data_int=pd.DataFrame({'Interaction AC':a1_len})
data_num=pd.DataFrame({'numbers':index1})
data1=pd.concat([data_int,data_num],axis=1)

nums=list(set(index1))
mut_nums=[]
nums_j=[]
for j in nums:
	num_j=len(data1[data1['numbers']==j])
	nums_j.append(num_j)
	mut_num=j
	mut_nums.append(mut_num)
mut_dat=pd.DataFrame({'mut_numbers':mut_nums})
int_dat=pd.DataFrame({'Interaction numbers':nums_j})
x=pd.concat([mut_dat,int_dat],axis=1)

#绘制统计柱状图
import matplotlib.pyplot as plt

# 绘制柱状图
plt.bar(mut_nums,nums_j)

# 添加标签和标题
plt.xlabel('Mutation numbers')
plt.ylabel('Interaction numbers')
plt.title('Decreasing mutation numbers in Reviewed data')

# 保存图片
plt.savefig('Reviewed_Decreasing.png')



#####################################################################################

#Increasing
index1=[]
index_1=[]
index_a1=[]
a1_len=[]
a1_list=list(a3.groupby('Interaction AC'))
for i in range(len(a1_list)):
	len_num=len(a1_list[i][1])
	index1.append(len_num)
	index_1=index_1+list(a1_list[i][1].index)
	a1_len.append(a1_list[i][0])
data_int=pd.DataFrame({'Interaction AC':a1_len})
data_num=pd.DataFrame({'numbers':index1})
data1=pd.concat([data_int,data_num],axis=1)

nums=list(set(index1))
mut_nums=[]
nums_j=[]
for j in nums:
	num_j=len(data1[data1['numbers']==j])
	nums_j.append(num_j)
	mut_num=j
	mut_nums.append(mut_num)
mut_dat=pd.DataFrame({'mut_numbers':mut_nums})
int_dat=pd.DataFrame({'Interaction numbers':nums_j})
x=pd.concat([mut_dat,int_dat],axis=1)

#绘制统计柱状图
import matplotlib.pyplot as plt

# 绘制柱状图
plt.bar(mut_nums,nums_j)

# 添加标签和标题
plt.xlabel('Mutation numbers')
plt.ylabel('Interaction numbers')
plt.title('Increasing mutation numbers in Reviewed data')

# 保存图片
plt.savefig('Reviewed_Increasing.png')


#####################################################################################

#No effect
index1=[]
index_1=[]
index_a1=[]
a1_len=[]
a1_list=list(a4.groupby('Interaction AC'))
for i in range(len(a1_list)):
	len_num=len(a1_list[i][1])
	index1.append(len_num)
	index_1=index_1+list(a1_list[i][1].index)
	a1_len.append(a1_list[i][0])
data_int=pd.DataFrame({'Interaction AC':a1_len})
data_num=pd.DataFrame({'numbers':index1})
data1=pd.concat([data_int,data_num],axis=1)

nums=list(set(index1))
mut_nums=[]
nums_j=[]
for j in nums:
	num_j=len(data1[data1['numbers']==j])
	nums_j.append(num_j)
	mut_num=j
	mut_nums.append(mut_num)
mut_dat=pd.DataFrame({'mut_numbers':mut_nums})
int_dat=pd.DataFrame({'Interaction numbers':nums_j})
x=pd.concat([mut_dat,int_dat],axis=1)

#绘制统计柱状图
import matplotlib.pyplot as plt

# 绘制柱状图
plt.bar(mut_nums,nums_j)

# 添加标签和标题
plt.xlabel('Mutation numbers')
plt.ylabel('Interaction numbers')
plt.title('No effect mutation numbers in Reviewed data')

# 保存图片
plt.savefig('Reviewed_No_effect.png')

#####################################################################################################################################################
#Unreviewed
a=pd.read_csv('Unreviewed_completed.tsv',header=0,sep='\t')
a1=a[a['type']=='Disrupting']
a2=a[a['type']=='Decreasing']
a3=a[a['type']=='Increasing']
a4=a[a['type']=='No effect']

#####################################################################################
#Disrupting
index1=[]
index_1=[]
index_a1=[]
a1_len=[]
a1_list=list(a1.groupby('Interaction AC'))
for i in range(len(a1_list)):
	len_num=len(a1_list[i][1])
	index1.append(len_num)
	index_1=index_1+list(a1_list[i][1].index)
	a1_len.append(a1_list[i][0])
data_int=pd.DataFrame({'Interaction AC':a1_len})
data_num=pd.DataFrame({'numbers':index1})
data1=pd.concat([data_int,data_num],axis=1)

nums=list(set(index1))
mut_nums=[]
nums_j=[]
for j in nums:
	num_j=len(data1[data1['numbers']==j])
	nums_j.append(num_j)
	mut_num=j
	mut_nums.append(mut_num)
mut_dat=pd.DataFrame({'mut_numbers':mut_nums})
int_dat=pd.DataFrame({'Interaction numbers':nums_j})
x=pd.concat([mut_dat,int_dat],axis=1)

#绘制统计柱状图
import matplotlib.pyplot as plt

# 绘制柱状图
plt.bar(mut_nums,nums_j)

# 添加标签和标题
plt.xlabel('Mutation numbers')
plt.ylabel('Interaction numbers')
plt.title('Disrupting mutation numbers in Unreviewed data')

# 保存图片
plt.savefig('Unreviewed_Disrupting.png')

#####################################################################################
#Decreasing
index1=[]
index_1=[]
index_a1=[]
a1_len=[]
a1_list=list(a2.groupby('Interaction AC'))
for i in range(len(a1_list)):
	len_num=len(a1_list[i][1])
	index1.append(len_num)
	index_1=index_1+list(a1_list[i][1].index)
	a1_len.append(a1_list[i][0])
data_int=pd.DataFrame({'Interaction AC':a1_len})
data_num=pd.DataFrame({'numbers':index1})
data1=pd.concat([data_int,data_num],axis=1)

nums=list(set(index1))
mut_nums=[]
nums_j=[]
for j in nums:
	num_j=len(data1[data1['numbers']==j])
	nums_j.append(num_j)
	mut_num=j
	mut_nums.append(mut_num)
mut_dat=pd.DataFrame({'mut_numbers':mut_nums})
int_dat=pd.DataFrame({'Interaction numbers':nums_j})
x=pd.concat([mut_dat,int_dat],axis=1)

#绘制统计柱状图
import matplotlib.pyplot as plt

# 绘制柱状图
plt.bar(mut_nums,nums_j)

# 添加标签和标题
plt.xlabel('Mutation numbers')
plt.ylabel('Interaction numbers')
plt.title('Decreasing mutation numbers in Unreviewed data')

# 保存图片
plt.savefig('Unreviewed_Decreasing.png')



#####################################################################################

#Increasing
index1=[]
index_1=[]
index_a1=[]
a1_len=[]
a1_list=list(a3.groupby('Interaction AC'))
for i in range(len(a1_list)):
	len_num=len(a1_list[i][1])
	index1.append(len_num)
	index_1=index_1+list(a1_list[i][1].index)
	a1_len.append(a1_list[i][0])
data_int=pd.DataFrame({'Interaction AC':a1_len})
data_num=pd.DataFrame({'numbers':index1})
data1=pd.concat([data_int,data_num],axis=1)

nums=list(set(index1))
mut_nums=[]
nums_j=[]
for j in nums:
	num_j=len(data1[data1['numbers']==j])
	nums_j.append(num_j)
	mut_num=j
	mut_nums.append(mut_num)
mut_dat=pd.DataFrame({'mut_numbers':mut_nums})
int_dat=pd.DataFrame({'Interaction numbers':nums_j})
x=pd.concat([mut_dat,int_dat],axis=1)

#绘制统计柱状图
import matplotlib.pyplot as plt

# 绘制柱状图
plt.bar(mut_nums,nums_j)

# 添加标签和标题
plt.xlabel('Mutation numbers')
plt.ylabel('Interaction numbers')
plt.title('Increasing mutation numbers in Unreviewed data')

# 保存图片
plt.savefig('Unreviewed_Increasing.png')


#####################################################################################

#No effect
index1=[]
index_1=[]
index_a1=[]
a1_len=[]
a1_list=list(a4.groupby('Interaction AC'))
for i in range(len(a1_list)):
	len_num=len(a1_list[i][1])
	index1.append(len_num)
	index_1=index_1+list(a1_list[i][1].index)
	a1_len.append(a1_list[i][0])
data_int=pd.DataFrame({'Interaction AC':a1_len})
data_num=pd.DataFrame({'numbers':index1})
data1=pd.concat([data_int,data_num],axis=1)

nums=list(set(index1))
mut_nums=[]
nums_j=[]
for j in nums:
	num_j=len(data1[data1['numbers']==j])
	nums_j.append(num_j)
	mut_num=j
	mut_nums.append(mut_num)
mut_dat=pd.DataFrame({'mut_numbers':mut_nums})
int_dat=pd.DataFrame({'Interaction numbers':nums_j})
x=pd.concat([mut_dat,int_dat],axis=1)

#绘制统计柱状图
import matplotlib.pyplot as plt

# 绘制柱状图
plt.bar(mut_nums,nums_j)

# 添加标签和标题
plt.xlabel('Mutation numbers')
plt.ylabel('Interaction numbers')
plt.title('No effect mutation numbers in Unreviewed data')

# 保存图片
plt.savefig('Unreviewed_No_effect.png')


#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################

#将所有类型的突变全部绘制在一张图中
import os
os.chdir('/data/liuyang/bishe_data/data_tongji')
import pandas as pd

#####################################################################################################################################################
#Reviewed
a=pd.read_csv('Reviewed_completed.tsv',header=0,sep='\t')
a1=a[a['type']=='Disrupting']
a2=a[a['type']=='Decreasing']
a3=a[a['type']=='Increasing']
a4=a[a['type']=='No effect']

#####################################################################################
#Disrupting
index1=[]
index_1=[]
index_a1=[]
a1_len=[]
a1_list=list(a1.groupby('Interaction AC'))
for i in range(len(a1_list)):
	len_num=len(a1_list[i][1])
	index1.append(len_num)
	index_1=index_1+list(a1_list[i][1].index)
	a1_len.append(a1_list[i][0])
data_int=pd.DataFrame({'Interaction AC':a1_len})
data_num=pd.DataFrame({'numbers':index1})
data1=pd.concat([data_int,data_num],axis=1)
nums=list(set(index1))
mut_nums=[]
nums_j=[]
for j in nums:
	num_j=len(data1[data1['numbers']==j])
	nums_j.append(num_j)
	mut_num=j
	mut_nums.append(mut_num)
mut_dat=pd.DataFrame({'mut_numbers':mut_nums})
int_dat=pd.DataFrame({'Interaction numbers':nums_j})
x1=pd.concat([mut_dat,int_dat],axis=1)
#
#Decreasing
index1=[]
index_1=[]
index_a1=[]
a1_len=[]
a1_list=list(a2.groupby('Interaction AC'))
for i in range(len(a1_list)):
	len_num=len(a1_list[i][1])
	index1.append(len_num)
	index_1=index_1+list(a1_list[i][1].index)
	a1_len.append(a1_list[i][0])
data_int=pd.DataFrame({'Interaction AC':a1_len})
data_num=pd.DataFrame({'numbers':index1})
data1=pd.concat([data_int,data_num],axis=1)
nums=list(set(index1))
mut_nums=[]
nums_j=[]
for j in nums:
	num_j=len(data1[data1['numbers']==j])
	nums_j.append(num_j)
	mut_num=j
	mut_nums.append(mut_num)
mut_dat=pd.DataFrame({'mut_numbers':mut_nums})
int_dat=pd.DataFrame({'Interaction numbers':nums_j})
x2=pd.concat([mut_dat,int_dat],axis=1)
#
#Increasing
index1=[]
index_1=[]
index_a1=[]
a1_len=[]
a1_list=list(a3.groupby('Interaction AC'))
for i in range(len(a1_list)):
	len_num=len(a1_list[i][1])
	index1.append(len_num)
	index_1=index_1+list(a1_list[i][1].index)
	a1_len.append(a1_list[i][0])
data_int=pd.DataFrame({'Interaction AC':a1_len})
data_num=pd.DataFrame({'numbers':index1})
data1=pd.concat([data_int,data_num],axis=1)

nums=list(set(index1))
mut_nums=[]
nums_j=[]
for j in nums:
	num_j=len(data1[data1['numbers']==j])
	nums_j.append(num_j)
	mut_num=j
	mut_nums.append(mut_num)
mut_dat=pd.DataFrame({'mut_numbers':mut_nums})
int_dat=pd.DataFrame({'Interaction numbers':nums_j})
x3=pd.concat([mut_dat,int_dat],axis=1)
#
#No effect
index1=[]
index_1=[]
index_a1=[]
a1_len=[]
a1_list=list(a4.groupby('Interaction AC'))
for i in range(len(a1_list)):
	len_num=len(a1_list[i][1])
	index1.append(len_num)
	index_1=index_1+list(a1_list[i][1].index)
	a1_len.append(a1_list[i][0])
data_int=pd.DataFrame({'Interaction AC':a1_len})
data_num=pd.DataFrame({'numbers':index1})
data1=pd.concat([data_int,data_num],axis=1)

nums=list(set(index1))
mut_nums=[]
nums_j=[]
for j in nums:
	num_j=len(data1[data1['numbers']==j])
	nums_j.append(num_j)
	mut_num=j
	mut_nums.append(mut_num)
mut_dat=pd.DataFrame({'mut_numbers':mut_nums})
int_dat=pd.DataFrame({'Interaction numbers':nums_j})
x4=pd.concat([mut_dat,int_dat],axis=1)

#构建全部数据
y1=pd.merge(x1,x2,how='outer',on='mut_numbers')
y1.rename(columns={'Interaction numbers_x':'Disrupting_mut_nums','Interaction numbers_y':'Decreasing_mut_nums'},inplace=True)
y2=pd.merge(y1,x3,how='outer',on='mut_numbers')
y3=pd.merge(y2,x4,how='outer',on='mut_numbers')
y3.rename(columns={'Interaction numbers_x':'Increasing_mut_nums','Interaction numbers_y':'No_effect_mut_nums'},inplace=True)
y3.fillna(0,inplace=True)
mut_nums=list(y3['mut_numbers'])
Disrupting_data=list(y3['Disrupting_mut_nums'])
Decreasing_data=list(y3['Decreasing_mut_nums'])
Increasing_data=list(y3['Increasing_mut_nums'])
No_effect_data=list(y3['No_effect_mut_nums'])

#绘制柱状图

import matplotlib.pyplot as plt
plt.bar(mut_nums,Disrupting_data,color='red',label='Disrupting')
plt.bar(mut_nums,Decreasing_data,color='orange',label='Decreasing')
plt.bar(mut_nums,Increasing_data,color='blue',label='Increasing')
plt.bar(mut_nums,No_effect_data,color='green',label='No effect')
plt.legend()

# 添加标签和标题
plt.xlabel('Mutation numbers')
plt.ylabel('Interaction numbers')
plt.title('Mutation numbers in Reviewed data')
# 保存图片
plt.savefig('Reviewed_Disrupting.png')


#绘制折线图
import matplotlib.pyplot as plt
plt.plot(mut_nums,Disrupting_data,color='red',label='Disrupting')
plt.plot(mut_nums,Decreasing_data,color='orange',linestyle='--',label='Decreasing')
plt.plot(mut_nums,Increasing_data,color='blue',linestyle=':',label='Increasing')
plt.plot(mut_nums,No_effect_data,color='green',linestyle='-',label='No effect')
plt.legend()
plt.xlabel('Mutation numbers')
plt.ylabel('Interaction numbers')
plt.title('Disrupting mutation numbers in Reviewed data')



total_width,n=0.2,4
width=total_width/n
import numpy as np
x=np.array(mut_nums1)
x=x-width*2

import matplotlib.pyplot as plt
plt.bar(x,Disrupting_data,color='red',width=width,label='Disrupting')
plt.bar(x+width,Decreasing_data,width=width,color='orange',label='Decreasing')
plt.bar(x+2*width,Increasing_data,width=width,color='blue',label='Increasing')
plt.bar(x+3*width,No_effect_data,width=width,color='green',label='No effect')
plt.legend()

# 添加标签和标题
plt.xlabel('log2(Mutation numbers)')
plt.ylabel('Interaction numbers')
plt.title('Mutation numbers in Reviewed data')


########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################

#统计map结构过程中删除多少数据
import os
os.chdir('/data/liuyang/bishe_data/get_pdb/Reviewed')
import pandas as pd
no_id=pd.read_csv('no_pdb_ids.txt',header=0,sep='\t')
no_id_list=list(no_id['no_pdb_id'])
a=pd.read_csv('Reviewed_completed.tsv',header=0,sep='\t')
a1=a[a['Affected protein AC'].isin(no_id_list)]
no_index=list(a1.index)
a2=a.drop(no_index)
#删7772个，剩余31457个



###########################################################################################################
#ceshi
#将数据分成几组跑
import os
os.chdir('/data/liuyang/bishe_data/get_pdb/Reviewed/ceshi')
import pandas as pd
x1=pd.read_csv('n_result1.tsv',header=0,sep='\t')
x2=pd.read_csv('n1.tsv',header=0,sep='\t')
x3=pd.read_csv('n2.tsv',header=0,sep='\t')
x4=pd.read_csv('n3.tsv',header=0,sep='\t')
a=pd.concat([x1,x2,x3,x4])
a.index=range(len(a))

#entity_id代表的是复合物分子中，元件的个数
#判断有多少个蛋白质参与，除了蛋白质单体双方之外，如果还有其他配体分子参与复合物，那么它们的entity_id就不是1和2
data_a=[]
for i in range(len(a)):
	a1=a['entity_id_0'][i]
	a2=a['entity_id_1'][i]
	temp_a=[a1,a2]
	data_a.append(temp_a)
data_dat=pd.DataFrame({'entity_group':data_a})
a_n=pd.concat([a,data_dat],axis=1)

type1=[]
for j in range(len(a_n)):
	if a_n['entity_group'][j]==[1,1] or a_n['entity_group'][j]==[2,1] or a_n['entity_group'][j]==[1,2]:
		type1.append(0)
	else:
		type1.append(1)
type_dat=pd.DataFrame({'ligand':type1})
a_m=pd.concat([a_n,type_dat],axis=1)
#判断是否有其他元件参与，通过ligand的0和1来判断，0代表的是无其他元件参与，1代表有其他元件参与。

#处理结构信息b
b=pd.read_csv('result1.tsv',header=0,sep='\t')
b1=pd.concat([b['uniprot_id'],b['pdb_id'],b['coverage'],b['resolution'],b['experimental_method']],axis=1)
b2=b1.drop_duplicates(subset=['uniprot_id','pdb_id','coverage','resolution','experimental_method'],keep='first')

#查看重复
b3=b1[b1.duplicated(subset=['uniprot_id','pdb_id','coverage','resolution','experimental_method'],keep='first')]

data_result=pd.merge(a_m,b1,how='inner',on='pdb_id')
data_result1=data_result.drop_duplicates(subset=['pdb_id','resolution','experimental_method'],keep='first')

#将横坐标取log



