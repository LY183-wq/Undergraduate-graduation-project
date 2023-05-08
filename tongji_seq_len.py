
#统计蛋白质中链长分布
#Reviewed

import os
os.chdir('/data/liuyang/bishe_data/tongji_seq_length')
import pandas as pd
#对整体序列长度分布进行统计
a=pd.read_csv('Reviewed_right.tsv',header=0,sep='\t')#46534
pro_A_length=[]
for i in range(len(a)):
	temp_length=len(a['A sequence'][i])
	pro_A_length.append(temp_length)
pro_A_dat=pd.DataFrame({'pro_A_length':pro_A_length})
a1=pd.concat([a,pro_A_dat],axis=1)

#绘制全部蛋白链长的密度分布图
import seaborn as sns
a2=a1.drop_duplicates(subset=['Affected protein AC'],keep='first')
a3=pd.concat([a2['Affected protein AC'],a2['pro_A_length']])
sns.set(style='ticks')
g = sns.FacetGrid(a3, height=4)
g.map(sns.kdeplot, 'pro_A_length', shade=True)
g.add_legend()




#统计未预测结构的部分数据
b=pd.read_csv('/data/liuyang/bishe_data/colab1/Reviewed_big',header=None)
b_list=list(b[0])
b1=a[a['Affected protein AC'].isin(b_list)]
b1.index=range(len(b1))
pro_B_length=[]
for j in range(len(b1)):
	temp_length1=len(b1['A sequence'][j])
	pro_B_length.append(temp_length1)
pro_B_dat=pd.DataFrame({'pro_A_length':pro_B_length})
b2=pd.concat([b1,pro_B_dat],axis=1)





#获取结构的pLDDT打分
#要求中位数大于70，突变位点大于70
import os
os.chdir('/data/liuyang/bishe_data/predict_structures')
import pandas as pd
import statistics

path1='/data/liuyang/bishe_data/colab1/all'
path1_data=os.listdir(path1)
#确定哪些是colabfold预测的结果（文件夹）
colab_data=[]
for i in path1_data:
    # 判断是否为文件夹，如果是，则输出文件夹的名称
    if os.path.isdir(os.path.join(path1, i)):
        colab_data.append(i)

#打开文件夹，获取其中的结果json文件，从这个文件中提取出pLDDT结果
def get_pLDDT_colab(x):
	import json
	with open(x) as f:
		content=f.read()
		data=json.loads(content)
		plDDT_data=data['plddt']
	return plDDT_data


#判断alphafold结构的pLDDT中位数是否大于70
def get_pLDDT_AlphaFold(pdb_file_path):
    plddt = []
    import statistics
    with open(pdb_file_path) as f:
        model_plddt=[]
        locs=[]
        for line in f:
            if line.startswith("ATOM"):
                plddt_value = float(line[61:66].strip())
                model_plddt.append(plddt_value)
                loc=int(line[23:28].strip())
                locs.append(loc)
            if line.startswith("ENDMDL"):
                break
    loc_dat=pd.DataFrame({'loc':locs})
    plddt_dat=pd.DataFrame({'plddt':model_plddt})
    x=pd.concat([loc_dat,plddt_dat],axis=1)
    x1=x.drop_duplicates(subset=['loc'],keep='first')
    y=list(x1['plddt'])
    y1=statistics.median(y)
    return y1



a=pd.read_csv('Reviewed_right.tsv',header=0,sep='\t')
b=pd.read_csv('Unreviewed_right.tsv',header=0,sep='\t')
a_list=list(a['Affected protein AC'])
b_list=list(b['Affected protein AC'])

bad_pdb=[]    #存放pLDDT中位数小于70的模型结构
bad_pdb1=[]
bad_pdb2=[]
bad_pdb3=[]

alpha_a=[x for x in a_list if x in a_data]
alpha_b=[y for y in b_list if y in a_data]
for i in range(len(b)):
    if b['Affected protein AC'][i] in alpha_b:
        file_i='/data/liuyang/bishe_data/alphafold/Unreviewed/'+b['Affected protein AC'][i]+'.pdb'
        plddt_median=get_pLDDT_AlphaFold(file_i)
        if plddt_median < 70:
            bad_pdb.append(b['Affected protein AC'][i])

for i in range(len(a)):
    if a['Affected protein AC'][i] in alpha_a:
        file_i='/data/liuyang/bishe_data/alphafold/Reviewed/'+a['Affected protein AC'][i]+'.pdb'
        plddt_median=get_pLDDT_AlphaFold(file_i)
        if plddt_median < 70:
            bad_pdb1.append(b['Affected protein AC'][i])

#删除colabfold效果差的
colab_a=[x[:x.index('_')] for x in colab_data if x[:x.index('_')] in a_list]
colab_b=[y[:y.index('_')] for y in colab_data if y[:y.index('_')] in b_list]
for j in colab_data:
	path2='/data/liuyang/bishe_data/colab1/all/'+j
	path2_data=os.listdir(path2)
	file_json=[y for y in path2_data if 'scores_rank' in y][0]
	file_name=path2+'/'+file_json
	plddt_j=get_pLDDT(file_name)
	plddt_median=statistics.median(plddt_j)
	if plddt_median<70:
		bad_pdb.append(j)
good_pdb=[x for x in colab_data if x not in bad_pdb]
good_pdb1=[y[:y.index('_')] for y in good_pdb]
bad_pdb1=[z[:z.index('_')] for z in bad_pdb]

a_good=a[a['Affected protein AC'].isin(good_pdb1)]
a_bad=a[a['Affected protein AC'].isin(bad_pdb1)]
b_good=b[b['Affected protein AC'].isin(good_pdb1)]
b_bad=b[b['Affected protein AC'].isin(bad_pdb1)]

import re
plddt_scores=[]    #存放所有数据的pLDDT
for m in range(len(a)):
	if a['Affected protein AC'][m] in good_pdb1:
		loc_m=int(re.findall(r'\d+',a['Mutations'][m])[0])
		pdb_file=[x for x in good_pdb if a['Affected protein AC'][m] in x][0]
		path_m='/data/liuyang/bishe_data/colab1/all/'+pdb_file
		path_m_data=os.listdir(path_m)
		file_m_name=[y for y in path_m_data if 'scores_rank' in y][0]
		file_m_json=path_m+'/'+file_m_name
		m_plddt=get_pLDDT(file_m_json)
		plddt_loc=m_plddt[loc_m-1]
		plddt_scores.append(plddt_loc)
	else:
		plddt_scores.append('no')

plddt_dat=pd.DataFrame({'plddt':plddt_scores})
a_plddt=pd.concat([a,plddt_dat],axis=1)

#Unreviewed,b
import re
plddt_scores=[]    #存放所有数据的pLDDT
for m in range(len(b)):
	if b['Affected protein AC'][m] in good_pdb1:
		loc_m=int(re.findall(r'\d+',b['Mutations'][m])[0])
		pdb_file=[x for x in good_pdb if b['Affected protein AC'][m] in x][0]
		path_m='/data/liuyang/bishe_data/colab1/all/'+pdb_file
		path_m_data=os.listdir(path_m)
		file_m_name=[y for y in path_m_data if 'scores_rank' in y][0]
		file_m_json=path_m+'/'+file_m_name
		m_plddt=get_pLDDT(file_m_json)
		plddt_loc=m_plddt[loc_m-1]
		plddt_scores.append(plddt_loc)
	else:
		plddt_scores.append('no')

plddt_dat=pd.DataFrame({'plddt':plddt_scores})
b_plddt=pd.concat([b,plddt_dat],axis=1)

os.chdir('/data/liuyang/bishe_data/colab1/all')
a_plddt.to_csv('Reviewed_plddt.tsv',index=False,sep='\t')
b_plddt.to_csv('Unreviewed_plddt.tsv',index=False,sep='\t')

a_bad_index=list(a[a['Affected protein AC'].isin(bad_pdb1)].index)
a_r=a.drop(a_bad_index)

b_bad_index=list(b[b['Affected protein AC'].isin(bad_pdb1)].index)
b_r=b.drop(b_bad_index)

a_r.to_csv('Reviewed_r.tsv',index=False,sep='\t')
b_r.to_csv('Unreviewed_r.tsv',index=False,sep='\t')

#绘制密度分布图

#保守性绘制图
import os
os.chdir('/data/liuyang/bishe_data/conservations')
import pandas as pd	
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
sns.set(style='white',font_scale=1.5)
g = sns.kdeplot(a1['provean_score'],shade=False,color='red')
g = sns.kdeplot(a2['provean_score'],shade=False,color='orange')
g = sns.kdeplot(a3['provean_score'],shade=False,color='blue')
g = sns.kdeplot(a4['provean_score'],shade=False,color='green')
plt.legend(['Disrupting', 'Decreasing', 'Increasing','No effect'],fontsize=15,bbox_to_anchor=(1.4, 1),loc='upper right')
plt.title('Density plot of 4 types Reviewed mutations Provean scores data',fontsize=15)
threshold = -2.5
plt.axvline(x=threshold, color='r', linestyle='--')

# 显示图形
plt.show()

#四类柱状图
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
fig = plt.figure()
plt.bar(x1, aa,width=width, label='Disrupting',color='red')  # 绘制第一组数据的柱形
plt.bar(x2, bb, width, label='Decreasing',color='orange')  # 绘制第二组数据的柱形
plt.bar(x3, cc, width, label='Increasing',color='blue')  # 绘制第二组数据的柱形
plt.bar(x4, dd, width, label='No effect',color='green')  # 绘制第二组数据的柱形
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
# 添加标签和标题
plt.xlabel('type',fontsize=15)
plt.ylabel('numbers',fontsize=15)
plt.title('conservative numbers Reviewed data',fontsize=15)
plt.xticks([1.3,3.3],xlabs)#设置横坐标
for m,n in zip(x1,aa):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=10)
for m,n in zip(x2,bb):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=10)
for m,n in zip(x3,cc):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=10)
for m,n in zip(x4,dd):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=10)







#对此部分可用结构数据，计算dssp溶剂可及表面积，以及scannet的probability
os.chdir('/data/liuyang/bishe_data/colab1/all')
path2='/data/liuyang/bishe_data/colab1/all'
path2_data=os.listdir(path2)
for i in good_pdb:
	path2_i='/data/liuyang/bishe_data/colab1/all/'+i
	path2_i_data=os.listdir(path2_i)
	i_pdb=[x for x in path2_i_data if x[-3:]=='pdb'][0]
	os.system('mkdssp -i '+'/data/liuyang/bishe_data/colab1/all/'+i+'/'+i_pdb+' -o '+'/data/liuyang/bishe_data/colab1/all/'+i_pdb[:i_pdb.index('_')]+'.dssp')
#生成61个dssp格式的结果文件

#重新统计溶剂可及表面积，这里使用蛋白质有结构的数据
import os
os.chdir('/data/liuyang/bishe_data')
import pandas as pd
import re
a=pd.read_csv('Unreviewed_pdb.tsv',header=0,sep='\t')
path1='/data/liuyang/bishe_data/alphafold/Reviewed'
path1_data=os.listdir(path1)
alphafold_pro=[x[:x.index('.')] for x in path1_data if x[-4:]=='dssp']
path2='/data/liuyang/bishe_data/colab1/all'
path2_data=os.listdir(path2)
colabfold_pro=[y[:y.index('.')] for y in path2_data if y[-4:]=='dssp']
path3='/data/liuyang/bishe_data/alphafold/Unreviewed'
path3_data=os.listdir(path3)
alphafold_pro1=[x[:x.index('.')] for x in path3_data if x[-4:]=='dssp']
os.chdir('/data/liuyang/bishe_data/colab1')
#只保留有alphafold和colabfold预测结构的数据
pre_pro=list(set(alphafold_pro+colabfold_pro+alphafold_pro1))
a_pro=a[a['Affected protein AC'].isin(pre_pro)]
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
	elif a_pro['Affected protein AC'][i] in alphafold_pro1:
		file_name_i='/data/liuyang/bishe_data/alphafold/Unreviewed/'+a_pro['Affected protein AC'][i]+'.dssp'
		loc_i=int(re.findall(r'\d+',a_pro['Mutations'][i])[0])
		dssp_i=cal_acc(file_name_i,str(loc_i))
		dssp_acc.append(dssp_i)
	elif a_pro['Affected protein AC'][i] in colabfold_pro:
		file_name_i='/data/liuyang/bishe_data/colab1/all/'+a_pro['Affected protein AC'][i]+'.dssp'
		loc_i=int(re.findall(r'\d+',a_pro['Mutations'][i])[0])
		dssp_i=cal_acc(file_name_i,str(loc_i))
		dssp_acc.append(dssp_i)

acc_dat=pd.DataFrame({'dssp_acc':dssp_acc})
a_dssp=pd.concat([a_pro,acc_dat],axis=1)
#添加类型，比值小于0.2的为核心
residue_location=[]
for p in range(len(a_dssp)):
	if a_dssp['dssp_acc'][p]<0.2:
		residue_location.append('Core')
	else:
		residue_location.append('Surface')
residue_dat=pd.DataFrame({'residue_location':residue_location})
a_result=pd.concat([a_dssp,residue_dat],axis=1)
os.chdir('/data/liuyang/bishe_data/alphafold')
a_result.to_csv('Reviewed_dssp1.tsv',index=False,sep='\t')

#统计分析
#密度分布图
import os
os.chdir('/data/liuyang/bishe_data/alphafold')
a=pd.read_csv('Reviewed_dssp1.tsv',header=0,sep='\t')
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
sns.set(style='white',font_scale=1.5)
g = sns.kdeplot(a1['dssp_acc'],shade=False,color='red')
g = sns.kdeplot(a2['dssp_acc'],shade=False,color='orange')
g = sns.kdeplot(a3['dssp_acc'],shade=False,color='blue')
g = sns.kdeplot(a4['dssp_acc'],shade=False,color='green')
plt.legend(['Disrupting', 'Decreasing', 'Increasing','No effect'],fontsize=15,bbox_to_anchor=(1.4, 1),loc='upper right')
plt.title('Density plot of 4 types Reviewed mutations data',fontsize=15)
threshold = 0.2
plt.axvline(x=threshold, color='r', linestyle='--')
# 显示图形
plt.show()

#四类柱状图
import os
os.chdir('/data/liuyang/bishe_data/alphafold')
x=pd.read_csv('Reviewed_dssp1.tsv',header=0,sep='\t')
a=x[x['type']=='Disrupting']
b=x[x['type']=='Decreasing']
c=x[x['type']=='Increasing']
d=x[x['type']=='No effect']
a1=len(a[a['residue_location']=='Core'])
a2=len(a[a['residue_location']=='Surface'])
b1=len(b[b['residue_location']=='Core'])
b2=len(b[b['residue_location']=='Surface'])
c1=len(c[c['residue_location']=='Core'])
c2=len(c[c['residue_location']=='Surface'])
d1=len(d[d['residue_location']=='Core'])
d2=len(d[d['residue_location']=='Surface'])
aa=[a1,a2]
bb=[b1,b2]
cc=[c1,c2]
dd=[d1,d2]
xlabs=['Core','Surface']
import matplotlib.pyplot as plt
width=0.2
x1=[1,3]
x2=[1.2,3.2]
x3=[1.4,3.4]
x4=[1.6,3.6]
fig = plt.figure()
plt.bar(x1, aa,width=width, label='Disrupting',color='red')  # 绘制第一组数据的柱形
plt.bar(x2, bb, width, label='Decreasing',color='orange')  # 绘制第二组数据的柱形
plt.bar(x3, cc, width, label='Increasing',color='blue')  # 绘制第二组数据的柱形
plt.bar(x4, dd, width, label='No effect',color='green')  # 绘制第二组数据的柱形
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
# 添加标签和标题
plt.xlabel('type',fontsize=15)
plt.ylabel('numbers',fontsize=15)
plt.title('Reviewed location data',fontsize=15)
plt.xticks([1.3,3.3],xlabs)#设置横坐标
for m,n in zip(x1,aa):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=10)
for m,n in zip(x2,bb):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=10)
for m,n in zip(x3,cc):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=10)
for m,n in zip(x4,dd):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=10)












#计算scannet的probability
import os
os.chdir('/data/liuyang/bishe_data/scannet/ScanNet-main')
import pandas as pd
codes1=[]
for i in good_pdb:
	path2='/data/liuyang/bishe_data/colab1/all/'+i
	path2_data=os.listdir(path2)
	pro_pdb=[y for y in path2_data if y[-3:]=='pdb'][0]
	codes1.append('python predict_bindingsites.py /data/liuyang/bishe_data/colab1/all/'+i+'/'+pro_pdb+' --noMSA')
codes1_dat=pd.DataFrame({'codes_i':codes1})
codes1_dat.to_csv('chongxin_colab.txt',index=False,sep='\t')
#将第一行改成#!/bin/bash，将后缀改为sh，批处理

import os
os.chdir('/data/liuyang/bishe_data')
import pandas as pd
import re
a=pd.read_csv('Unreviewed_pdb.tsv',header=0,sep='\t')
path1='/data/liuyang/bishe_data/alphafold/Reviewed'
path1_data=os.listdir(path1)
alphafold_pro=[x[:x.index('.')] for x in path1_data if x[-4:]=='dssp']
path2='/data/liuyang/bishe_data/colab1/all'
path2_data=os.listdir(path2)
colabfold_pro=[y[:y.index('.')] for y in path2_data if y[-4:]=='dssp']
path3='/data/liuyang/bishe_data/alphafold/Unreviewed'
path3_data=os.listdir(path3)
alphafold_pro1=[x[:x.index('.')] for x in path3_data if x[-4:]=='dssp']
os.chdir('/data/liuyang/bishe_data/colab1')
#只保留有alphafold和colabfold预测结构的数据
pre_pro=list(set(alphafold_pro+colabfold_pro+alphafold_pro1))
a_pro=a[a['Affected protein AC'].isin(pre_pro)]
a_pro.index=range(len(a_pro))



path1='/data/liuyang/bishe_data/scannet/ScanNet-main/predictions'
path1_data=os.listdir(path1)
probability_data=[]
for i in range(len(a_pro)):
	loc1=int(re.findall(r'\d+',a_pro['Mutations'][i])[0])
	for j in path1_data:
		if j[:j.index('_')]==a_pro['Affected protein AC'][i]:
			path2='/data/liuyang/bishe_data/scannet/ScanNet-main/predictions/'+j
			path2_data=os.listdir(path2)
			file_i=[x for x in path2_data if x[-3:]=='csv'][0]
			file=pd.read_csv(path2+'/'+file_i,header=0,sep=',')
			probability_i=list(file[file['Residue Index']==loc1]['Binding site probability'])[0]
			probability_data.append(probability_i)

proba_dat=pd.DataFrame({'probability':probability_data})
ax=pd.concat([a_pro,proba_dat],axis=1)

#判断是否为binding site，1表示是binding-site
binding_site=[]
for i in range(len(ax)):
	if ax['probability'][i]>0.5:
		binding_site.append(1)
	else:
		binding_site.append(0)
binding_dat=pd.DataFrame({'binding-site':binding_site})
ay=pd.concat([ax,binding_dat],axis=1)
os.chdir('/data/liuyang/bishe_data')
ay.to_csv('Unreviewed_scannet.tsv',index=False,sep='\t')
############################################################



#统计分析
#密度分布图
import os
os.chdir('/data/liuyang/bishe_data')
a=pd.read_csv('Reviewed_scannet.tsv',header=0,sep='\t')
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
sns.set(style='white',font_scale=1.5)
g = sns.kdeplot(a1['probability'],shade=False,color='red')
g = sns.kdeplot(a2['probability'],shade=False,color='orange')
g = sns.kdeplot(a3['probability'],shade=False,color='blue')
g = sns.kdeplot(a4['probability'],shade=False,color='green')
plt.legend(['Disrupting', 'Decreasing', 'Increasing','No effect'],fontsize=15,bbox_to_anchor=(1.4, 1),loc='upper right')
plt.title('Density plot of 4 types Reviewed mutations probability data',fontsize=15)
threshold = 0.5
plt.axvline(x=threshold, color='r', linestyle='--')
# 显示图形
plt.show()

#四类柱状图
import os
os.chdir('/data/liuyang/bishe_data')
x=pd.read_csv('Reviewed_scannet.tsv',header=0,sep='\t')
a=x[x['type']=='Disrupting']
b=x[x['type']=='Decreasing']
c=x[x['type']=='Increasing']
d=x[x['type']=='No effect']
a1=len(a[a['binding-site']==1])
a2=len(a[a['binding-site']==0])
b1=len(b[b['binding-site']==1])
b2=len(b[b['binding-site']==0])
c1=len(c[c['binding-site']==1])
c2=len(c[c['binding-site']==0])
d1=len(d[d['binding-site']==1])
d2=len(d[d['binding-site']==0])
aa=[a1,a2]
bb=[b1,b2]
cc=[c1,c2]
dd=[d1,d2]
xlabs=['binding-site','Not binding-site']
import matplotlib.pyplot as plt
width=0.2
x1=[1,3]
x2=[1.2,3.2]
x3=[1.4,3.4]
x4=[1.6,3.6]
fig = plt.figure()
plt.bar(x1, aa,width=width, label='Disrupting',color='red')  # 绘制第一组数据的柱形
plt.bar(x2, bb, width, label='Decreasing',color='orange')  # 绘制第二组数据的柱形
plt.bar(x3, cc, width, label='Increasing',color='blue')  # 绘制第二组数据的柱形
plt.bar(x4, dd, width, label='No effect',color='green')  # 绘制第二组数据的柱形
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
# 添加标签和标题
plt.xlabel('type',fontsize=15)
plt.ylabel('numbers',fontsize=15)
plt.title('binding-site numbers Reviewed data',fontsize=15)
plt.xticks([1.3,3.3],xlabs)#设置横坐标
for m,n in zip(x1,aa):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=10)
for m,n in zip(x2,bb):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=10)
for m,n in zip(x3,cc):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=10)
for m,n in zip(x4,dd):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=10)


#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################

#绘图
import os
os.chdir('/data/liuyang/PNG')
import pandas as pd	
a=pd.read_csv('Reviewed_provean_score.tsv',header=0,sep='\t')
ax=a.rename(columns={'provean_score':'Provean'})
a=ax
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
plt.figure(figsize=(15,12))
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
sns.set(style='white',font_scale=5)
g = sns.kdeplot(a1['Provean'],shade=False,color='red')
g = sns.kdeplot(a2['Provean'],shade=False,color='orange')
g = sns.kdeplot(a3['Provean'],shade=False,color='blue')
g = sns.kdeplot(a4['Provean'],shade=False,color='green')
threshold = -2.5
plt.axvline(x=threshold, color='r', linestyle='--')
ax = plt.gca()
line = ax.lines[0]
line.set_linewidth(5)
line = ax.lines[1]
line.set_linewidth(5)
line = ax.lines[2]
line.set_linewidth(5)
line = ax.lines[3]
line.set_linewidth(5)
plt.savefig("Reviewed_provean1.png",dpi=500)


import os
os.chdir('/data/liuyang/PNG')
import pandas as pd	
a=pd.read_csv('Unreviewed_provean_score.tsv',header=0,sep='\t')
ax=a.rename(columns={'provean_score':'Provean'})
a=ax
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
plt.figure(figsize=(15,12))
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
sns.set(style='white',font_scale=5)
g = sns.kdeplot(a1['Provean'],shade=False,color='red')
g = sns.kdeplot(a2['Provean'],shade=False,color='orange')
g = sns.kdeplot(a3['Provean'],shade=False,color='blue')
g = sns.kdeplot(a4['Provean'],shade=False,color='green')
plt.legend(['Disrupting', 'Decreasing', 'Increasing','No effect'],fontsize=28,bbox_to_anchor=(1, 1),loc='upper right')
threshold = -2.5
plt.axvline(x=threshold, color='r', linestyle='--')
ax = plt.gca()
line = ax.lines[0]
line.set_linewidth(5)
line = ax.lines[1]
line.set_linewidth(5)
line = ax.lines[2]
line.set_linewidth(5)
line = ax.lines[3]
line.set_linewidth(5)

plt.savefig("Unreviewed_provean1.png",dpi=500)


import os
os.chdir('/data/liuyang/PNG')
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
xlabs=['Conserved residue','Non-conserved residue']
import matplotlib.pyplot as plt
width=0.2
x1=[1,3]
x2=[1.2,3.2]
x3=[1.4,3.4]
x4=[1.6,3.6]
fig = plt.figure()
plt.figure(figsize=(30,13))
plt.bar(x1, aa,width=width, label='Disrupting',color='red')  # 绘制第一组数据的柱形
plt.bar(x2, bb, width, label='Decreasing',color='orange')  # 绘制第二组数据的柱形
plt.bar(x3, cc, width, label='Increasing',color='blue')  # 绘制第二组数据的柱形
plt.bar(x4, dd, width, label='No effect',color='green')  # 绘制第二组数据的柱形
# 添加标签和标题
plt.ylabel('Numbers of mutation',fontsize=40)
plt.xticks([1.3,3.3],xlabs,fontsize=40)#设置横坐标
for m,n in zip(x1,aa):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)
for m,n in zip(x2,bb):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)
for m,n in zip(x3,cc):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)
for m,n in zip(x4,dd):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)

plt.savefig('Reviewed_provean2.png',dpi=500)




import os
os.chdir('/data/liuyang/PNG')
import pandas as pd	
x=pd.read_csv('Unreviewed_provean_score.tsv',header=0,sep='\t')
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
xlabs=['Conserved residue','Non-conserved residue']
import matplotlib.pyplot as plt
width=0.2
x1=[1,3]
x2=[1.2,3.2]
x3=[1.4,3.4]
x4=[1.6,3.6]
fig = plt.figure()
plt.figure(figsize=(32,12))
plt.bar(x1, aa,width=width, label='Disrupting',color='red')  # 绘制第一组数据的柱形
plt.bar(x2, bb, width, label='Decreasing',color='orange')  # 绘制第二组数据的柱形
plt.bar(x3, cc, width, label='Increasing',color='blue')  # 绘制第二组数据的柱形
plt.bar(x4, dd, width, label='No effect',color='green')  # 绘制第二组数据的柱形
plt.legend(bbox_to_anchor=(0.7, 1), loc='upper left')
# 添加标签和标题
plt.ylabel('Numbers of mutation',fontsize=40)
plt.xticks([1.3,3.3],xlabs,fontsize=40)#设置横坐标
for m,n in zip(x1,aa):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)
for m,n in zip(x2,bb):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)
for m,n in zip(x3,cc):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)
for m,n in zip(x4,dd):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)

plt.savefig('Unreviewed_provean2.png',dpi=500)



##################################################################################################################################

import os
os.chdir('/data/liuyang/PNG')
a=pd.read_csv('Reviewed_dssp1.tsv',header=0,sep='\t')
ax=a.rename(columns={'dssp_acc':'Dssp'})
a=ax
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
plt.figure(figsize=(15,12))
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
sns.set(style='white',font_scale=5)
g = sns.kdeplot(a1['Dssp'],shade=False,color='red')
g = sns.kdeplot(a2['Dssp'],shade=False,color='orange')
g = sns.kdeplot(a3['Dssp'],shade=False,color='blue')
g = sns.kdeplot(a4['Dssp'],shade=False,color='green')
threshold = 0.2
plt.axvline(x=threshold, color='r', linestyle='--')
ax = plt.gca()
line = ax.lines[0]
line.set_linewidth(5)
line = ax.lines[1]
line.set_linewidth(5)
line = ax.lines[2]
line.set_linewidth(5)
line = ax.lines[3]
line.set_linewidth(5)
plt.savefig("Reviewed_dssp1.png",dpi=500)

import os
os.chdir('/data/liuyang/PNG')
a=pd.read_csv('Unreviewed_dssp1.tsv',header=0,sep='\t')
ax=a.rename(columns={'dssp_acc':'Dssp'})
a=ax
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
plt.figure(figsize=(15,12))
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
sns.set(style='white',font_scale=5)
g = sns.kdeplot(a1['Dssp'],shade=False,color='red')
g = sns.kdeplot(a2['Dssp'],shade=False,color='orange')
g = sns.kdeplot(a3['Dssp'],shade=False,color='blue')
g = sns.kdeplot(a4['Dssp'],shade=False,color='green')
plt.legend(['Disrupting', 'Decreasing', 'Increasing','No effect'],fontsize=40,bbox_to_anchor=(1, 1),loc='upper right')
threshold = 0.2
plt.axvline(x=threshold, color='r', linestyle='--')
ax = plt.gca()
line = ax.lines[0]
line.set_linewidth(5)
line = ax.lines[1]
line.set_linewidth(5)
line = ax.lines[2]
line.set_linewidth(5)
line = ax.lines[3]
line.set_linewidth(5)
plt.savefig("Unreviewed_dssp1.png",dpi=500)


import os
os.chdir('/data/liuyang/PNG')
x=pd.read_csv('Reviewed_dssp1.tsv',header=0,sep='\t')
a=x[x['type']=='Disrupting']
b=x[x['type']=='Decreasing']
c=x[x['type']=='Increasing']
d=x[x['type']=='No effect']
a1=len(a[a['residue_location']=='Core'])
a2=len(a[a['residue_location']=='Surface'])
b1=len(b[b['residue_location']=='Core'])
b2=len(b[b['residue_location']=='Surface'])
c1=len(c[c['residue_location']=='Core'])
c2=len(c[c['residue_location']=='Surface'])
d1=len(d[d['residue_location']=='Core'])
d2=len(d[d['residue_location']=='Surface'])
aa=[a1,a2]
bb=[b1,b2]
cc=[c1,c2]
dd=[d1,d2]
xlabs=['Core','Surface']
import matplotlib.pyplot as plt
width=0.2
x1=[1,3]
x2=[1.2,3.2]
x3=[1.4,3.4]
x4=[1.6,3.6]
fig = plt.figure()
plt.figure(figsize=(30,12))
plt.bar(x1, aa,width=width, label='Disrupting',color='red')  # 绘制第一组数据的柱形
plt.bar(x2, bb, width, label='Decreasing',color='orange')  # 绘制第二组数据的柱形
plt.bar(x3, cc, width, label='Increasing',color='blue')  # 绘制第二组数据的柱形
plt.bar(x4, dd, width, label='No effect',color='green')  # 绘制第二组数据的柱形
# 添加标签和标题
plt.ylabel('Numbers of mutations',fontsize=40)
plt.xticks([1.3,3.3],xlabs)#设置横坐标
for m,n in zip(x1,aa):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)
for m,n in zip(x2,bb):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)
for m,n in zip(x3,cc):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)
for m,n in zip(x4,dd):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)
plt.savefig("Reviewed_dssp2.png",dpi=500)


import os
os.chdir('/data/liuyang/PNG')
x=pd.read_csv('Unreviewed_dssp1.tsv',header=0,sep='\t')
a=x[x['type']=='Disrupting']
b=x[x['type']=='Decreasing']
c=x[x['type']=='Increasing']
d=x[x['type']=='No effect']
a1=len(a[a['residue_location']=='Core'])
a2=len(a[a['residue_location']=='Surface'])
b1=len(b[b['residue_location']=='Core'])
b2=len(b[b['residue_location']=='Surface'])
c1=len(c[c['residue_location']=='Core'])
c2=len(c[c['residue_location']=='Surface'])
d1=len(d[d['residue_location']=='Core'])
d2=len(d[d['residue_location']=='Surface'])
aa=[a1,a2]
bb=[b1,b2]
cc=[c1,c2]
dd=[d1,d2]
xlabs=['Core','Surface']
import matplotlib.pyplot as plt
width=0.2
x1=[1,3]
x2=[1.2,3.2]
x3=[1.4,3.4]
x4=[1.6,3.6]
fig = plt.figure()
plt.figure(figsize=(30,15))
plt.bar(x1, aa,width=width, label='Disrupting',color='red')  # 绘制第一组数据的柱形
plt.bar(x2, bb, width, label='Decreasing',color='orange')  # 绘制第二组数据的柱形
plt.bar(x3, cc, width, label='Increasing',color='blue')  # 绘制第二组数据的柱形
plt.bar(x4, dd, width, label='No effect',color='green')  # 绘制第二组数据的柱形
plt.legend(bbox_to_anchor=(0.1, 1), loc='upper left',fontsize=50)
# 添加标签和标题
plt.ylabel('Numbers of mutations',fontsize=40)
plt.xticks([1.3,3.3],xlabs)#设置横坐标
for m,n in zip(x1,aa):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)
for m,n in zip(x2,bb):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)
for m,n in zip(x3,cc):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)
for m,n in zip(x4,dd):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)
plt.savefig("Unreviewed_dssp2.png",dpi=500)


#################################################################################################################################################################
import os
os.chdir('/data/liuyang/PNG')
a=pd.read_csv('Reviewed_scannet.tsv',header=0,sep='\t')
ax=a.rename(columns={'probability':'Probability'})
a=ax
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
plt.figure(figsize=(15,12))
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
sns.set(style='white',font_scale=5)
g = sns.kdeplot(a1['Probability'],shade=False,color='red')
g = sns.kdeplot(a2['Probability'],shade=False,color='orange')
g = sns.kdeplot(a3['Probability'],shade=False,color='blue')
g = sns.kdeplot(a4['Probability'],shade=False,color='green')
threshold = 0.5
plt.axvline(x=threshold, color='r', linestyle='--')
ax = plt.gca()
line = ax.lines[0]
line.set_linewidth(5)
line = ax.lines[1]
line.set_linewidth(5)
line = ax.lines[2]
line.set_linewidth(5)
line = ax.lines[3]
line.set_linewidth(5)
plt.savefig('Reviewed_scannet1.png',dpi=500)

import os
os.chdir('/data/liuyang/PNG')
a=pd.read_csv('Unreviewed_scannet.tsv',header=0,sep='\t')
ax=a.rename(columns={'probability':'Probability'})
a=ax
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
plt.figure(figsize=(15,12))
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
sns.set(style='white',font_scale=5)
g = sns.kdeplot(a1['Probability'],shade=False,color='red')
g = sns.kdeplot(a2['Probability'],shade=False,color='orange')
g = sns.kdeplot(a3['Probability'],shade=False,color='blue')
g = sns.kdeplot(a4['Probability'],shade=False,color='green')
plt.legend(['Disrupting', 'Decreasing', 'Increasing','No effect'],fontsize=40,bbox_to_anchor=(1, 1),loc='upper right')
threshold = 0.5
plt.axvline(x=threshold, color='r', linestyle='--')
ax = plt.gca()
line = ax.lines[0]
line.set_linewidth(5)
line = ax.lines[1]
line.set_linewidth(5)
line = ax.lines[2]
line.set_linewidth(5)
line = ax.lines[3]
line.set_linewidth(5)
plt.savefig('Unreviewed_scannet1.png',dpi=500)

import os
os.chdir('/data/liuyang/PNG')
x=pd.read_csv('Reviewed_scannet.tsv',header=0,sep='\t')
a=x[x['type']=='Disrupting']
b=x[x['type']=='Decreasing']
c=x[x['type']=='Increasing']
d=x[x['type']=='No effect']
a1=len(a[a['binding-site']==1])
a2=len(a[a['binding-site']==0])
b1=len(b[b['binding-site']==1])
b2=len(b[b['binding-site']==0])
c1=len(c[c['binding-site']==1])
c2=len(c[c['binding-site']==0])
d1=len(d[d['binding-site']==1])
d2=len(d[d['binding-site']==0])
aa=[a1,a2]
bb=[b1,b2]
cc=[c1,c2]
dd=[d1,d2]
xlabs=['Binding-site','Not binding-site']
import matplotlib.pyplot as plt
width=0.2
x1=[1,3]
x2=[1.2,3.2]
x3=[1.4,3.4]
x4=[1.6,3.6]
fig = plt.figure()
plt.figure(figsize=(30,12))
plt.bar(x1, aa,width=width, label='Disrupting',color='red')  # 绘制第一组数据的柱形
plt.bar(x2, bb, width, label='Decreasing',color='orange')  # 绘制第二组数据的柱形
plt.bar(x3, cc, width, label='Increasing',color='blue')  # 绘制第二组数据的柱形
plt.bar(x4, dd, width, label='No effect',color='green')  # 绘制第二组数据的柱形

# 添加标签和标题
plt.ylabel('Numbers of mutations',fontsize=40)
plt.xticks([1.3,3.3],xlabs,fontsize=40)#设置横坐标
for m,n in zip(x1,aa):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)
for m,n in zip(x2,bb):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)
for m,n in zip(x3,cc):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)
for m,n in zip(x4,dd):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)
plt.savefig('Reviewed_scannet2.png',dpi=500)



import os
os.chdir('/data/liuyang/PNG')
x=pd.read_csv('Unreviewed_scannet.tsv',header=0,sep='\t')
a=x[x['type']=='Disrupting']
b=x[x['type']=='Decreasing']
c=x[x['type']=='Increasing']
d=x[x['type']=='No effect']
a1=len(a[a['binding-site']==1])
a2=len(a[a['binding-site']==0])
b1=len(b[b['binding-site']==1])
b2=len(b[b['binding-site']==0])
c1=len(c[c['binding-site']==1])
c2=len(c[c['binding-site']==0])
d1=len(d[d['binding-site']==1])
d2=len(d[d['binding-site']==0])
aa=[a1,a2]
bb=[b1,b2]
cc=[c1,c2]
dd=[d1,d2]
xlabs=['Binding-site','Not binding-site']
import matplotlib.pyplot as plt
width=0.2
x1=[1,3]
x2=[1.2,3.2]
x3=[1.4,3.4]
x4=[1.6,3.6]
fig = plt.figure()
plt.figure(figsize=(30,12))
plt.bar(x1, aa,width=width, label='Disrupting',color='red')  # 绘制第一组数据的柱形
plt.bar(x2, bb, width, label='Decreasing',color='orange')  # 绘制第二组数据的柱形
plt.bar(x3, cc, width, label='Increasing',color='blue')  # 绘制第二组数据的柱形
plt.bar(x4, dd, width, label='No effect',color='green')  # 绘制第二组数据的柱形
plt.legend(bbox_to_anchor=(0.3, 1), loc='upper left',fontsize=50)
# 添加标签和标题
plt.ylabel('Numbers of mutations',fontsize=40)
plt.xticks([1.3,3.3],xlabs,fontsize=40)#设置横坐标
for m,n in zip(x1,aa):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)
for m,n in zip(x2,bb):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)
for m,n in zip(x3,cc):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)
for m,n in zip(x4,dd):   #柱子上的数字显示
	plt.text(m,n,'%.f'%n,ha='center',va='bottom',fontsize=40)
plt.savefig('Unreviewed_scannet2.png',dpi=500)

