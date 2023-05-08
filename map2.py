import os
os.chdir('/data/liuyang/bishe_data/PremPS3')
import pandas as pd

a=pd.read_csv('S4191.txt',header=0,sep='\t')

Chains=[]
for i in range(len(a)):
	chain1=[x+'_1' for x in a['Partner1'][i]]
	chain2=[y+'_1' for y in a['Partner2'][i]]
	chain_i=chain1+chain2
	chains_i='.'.join(chain_i)
	Chains.append(chains_i)
Chains_dat=pd.DataFrame({'Chains':Chains})
PDBfile=a['PDB id'].str.cat(['.pdb'] * len(a),sep='')
MutChain=a['Mutated Chain'].str.cat(['_1']*len(a))
pdb_chains=a['PDB id'].str.lower().str.cat(['_'] * len(a),sep='') + a['Mutated Chain'].str.lower()

a1=pd.concat([PDBfile,Chains_dat,MutChain,a['Mutation(s)_PDB'],pdb_chains],axis=1)
a1.columns.values[4]='pdb_chains'
a1=a1.rename(columns={'PDB id':'PDBfile','Mutated Chain':'MutChain','Mutation(s)_PDB':'Mutation_PDB'})

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

result_data3.to_csv('S4191_before_PremPS.txt',index=False,sep='\t')


import os
os.chdir('/data/liuyang')
import pandas as pd
a=pd.read_csv('S4191.txt',header=0,sep='\t')
b=pd.read_csv('S4191_before_PremPS.txt',header=0,sep='\t')
b1=b.drop_duplicates(subset=['isPI'],keep='first')
b2=pd.concat([b1['PDBfile'].str.split('.').str[0],b1['MutChain'].str.split('_').str[0],b1['isPI']],axis=1)
b3=b2.rename(columns={b2.columns[0]:'PDB id',b2.columns[1]:'Mutated Chain'})
a1=pd.merge(a,b3,how='left',on=['PDB id','Mutated Chain'])

#统计蛋白A的链长,x为数据（即上面a4），y为路径（即PremPS结果所在路径）
def get_seq(x,y):
	A_length=[]
	for i in range(len(x)):
		mutchain_i=x['Mutated Chain'][i]
		PI_i=x['isPI'][i]
		path_i=y+PI_i+'out'
		all_i=os.listdir(path_i)
		file1=x['PDB id'][i].upper()+'_'+mutchain_i+'_1.seq'
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
		mutation_i=x['Mutation(s)_PDB'][i]
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

path4='/data/liuyang/Skempi2_PremPS/'
a5=get_seq(a1,path4)
a5_qu=a5[a5['A length']=='NO seq']
qu_index=list(a5_qu.index)
a6=a5.drop(qu_index)
a6.index=range(len(a6))


path5='/data/liuyang/Skempi2_PremPS/'
a7=get_PremPS(a6,path5)


a7.to_csv('S4191_PremPS_before_result.tsv',index=False,sep='\t')


import os
os.chdir('/data/liuyang')
import pandas as pd
a=pd.read_csv('S4191_PremPS_before_result.tsv',header=0,sep='\t')
a1=a[a['PremPS']=='NO']
a2=a[a['PremPS']=='NONE PDB']
qu_index=list(a1.index)+list(a2.index)
a_result=a.drop(qu_index)
a_result.index=range(len(a_result))
a_result.to_csv('S4191_PremPS_result.tsv',index=False,sep='\t')