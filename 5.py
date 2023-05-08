# -*- coding: utf-8 -*-
# @Author: Liu Yang
# @Date:   2023-03-11 09:48:53
# @Last Modified by:   Liu Yang
# @Last Modified time: 2023-03-16 10:22:02


#################################################################################################
#处理Unreviewed数据

import os
os.chdir('D://abishe//Fifth//Unreviewed')
import pandas as pd

data=pd.read_csv('Unreviewed_resullt.tsv',header=0,sep='\t')


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
	result.to_csv('Unreviewed_shaixuan3.tsv',index=False,sep='\t')#12538


#####################################################################################################################

#查找Affected protein AC有没有isoform
shuju=pd.read_csv('Unreviewed_shaixuan3.tsv',header=0,sep='\t')#12813
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
iso_index1=[]
P_ID=list(set(P_ID))
for p in range(len(shuju)):
	iso1=shuju['Affected protein AC'][p][10:]
	if iso1 in P_ID:
		iso_index1.append(p)
data_iso=shuju.loc[iso_index1]
data_zc=shuju.drop(iso_index1)
data_iso.to_csv('Unreviewed_iso.tsv',index=False,sep='\t')#3105
data_zc.to_csv('Unreviewed_zc.tsv',index=False,sep='\t')#9433


##################################################################################################################################################33
#先将全是canoincal的Unreviewed数据map序列sequence
a=pd.read_csv('Unreviewed_zc.tsv',header=0,sep='\t')
Un_pro=[]
for i in range(len(a)):
	a1=eval(a['Partners_1'][i])[0]
	a2=eval(a['Partners_1'][i])[1]
	Un_pro.append(a1)
	Un_pro.append(a2)
Un_pro1=list(set(Un_pro))
Un_dat=pd.DataFrame({'Protein ID':Un_pro1})
Un_dat.to_csv('Unreviewed_zc_ID.txt',index=False,sep='\t')



#匹配序列
a=pd.read_csv('Unreviewed_zc.tsv',header=0,sep='\t')
#读取fasta格式的uniprot现有且审查蛋白数据,并进行处理
f=open('Unreviewed_zc_sequence.fasta')
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
ID1=P_ID

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
ID2=P_ID



seqs1=[]
seqs2=[]
A_index=[]
B_index=[]
all_index=[]
for i in range(len(a)):
	temp1=eval(a['Partners_1'][i])[0][10:]
	temp2=eval(a['Partners_1'][i])[1][10:]
	if temp1 in ID2 and temp2 not in ID2:
		A_index.append(i)
		seq1=list(all_isoform[all_isoform['isoform_ID']==temp1]['isoform_seq'])[0]
		if '-' in temp2:
			seq2=list(Uniprot[Uniprot['Protein_ID']==temp2[:temp2.index('-')]]['Protein_seq'])[0]
		else:
			seq2=list(Uniprot[Uniprot['Protein_ID']==temp2]['Protein_seq'])[0]
		seqs1.append(seq1)
		seqs2.append(seq2)
		continue
	if temp1 not in ID2 and temp2 in ID2:
		B_index.append(i)
		if '-' in temp1:
			seq1=list(Uniprot[Uniprot['Protein_ID']==temp1[:temp1.index('-')]]['Protein_seq'])[0]
		else:
			seq1=list(Uniprot[Uniprot['Protein_ID']==temp1]['Protein_seq'])[0]
		seq2=list(all_isoform[all_isoform['isoform_ID']==temp2]['isoform_seq'])[0]
		seqs1.append(seq1)
		seqs2.append(seq2)
		continue
	if temp1 in ID2 and temp2 in ID2:
		all_index.append(i)
		seq1=list(all_isoform[all_isoform['isoform_ID']==temp1]['isoform_seq'])[0]
		seq2=list(all_isoform[all_isoform['isoform_ID']==temp2]['isoform_seq'])[0]
		seqs1.append(seq1)
		seqs2.append(seq2)
	else:
		if '-' in temp1:
			seq1=list(Uniprot[Uniprot['Protein_ID']==temp1[:temp1.index('-')]]['Protein_seq'])[0]
		else:
			seq1=list(Uniprot[Uniprot['Protein_ID']==temp1]['Protein_seq'])[0]
		if '-' in temp2:
			seq2=list(Uniprot[Uniprot['Protein_ID']==temp2[:temp2.index('-')]]['Protein_seq'])[0]
		else:
			seq2=list(Uniprot[Uniprot['Protein_ID']==temp2]['Protein_seq'])[0]
		seqs1.append(seq1)
		seqs2.append(seq2)
seq1_dat=pd.DataFrame({'A sequence':seqs1})
seq2_dat=pd.DataFrame({'B sequence':seqs2})
data=pd.concat([a,seq1_dat,seq2_dat],axis=1)
data.to_csv('Unreviewed_zc_seq.tsv',index=False,sep='\t')

#保存有效数据

data1=pd.read_csv('Unreviewed_zc_seq.tsv',header=0,sep='\t')
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
g=data1['A sequence']
h=data1['B sequence']
data2=pd.concat([d,c,mut_data,e,f,b,g,h],axis=1)
data2.to_csv('Unreviewed_zc_result.tsv',index=False,sep='\t')




################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################


#将isoform匹配到canonical上
a=pd.read_csv('Unreviewed_iso.tsv',header=0,sep='\t')
for i in range(len(a)):
	a['Affected protein AC'][i]=a['Affected protein AC'][i][10:]
a1=a['Affected protein AC']
b=list(set(a1))
iso=pd.DataFrame({'iso ID':b})
iso.to_csv('iso_id.txt',index=False,sep='\t')


#读取fasta格式的uniprot现有且审查蛋白数据,并进行处理
f=open('iso_sequence.fasta')
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
iso_data=pd.DataFrame(data_uniprot)
f.close()


iso_id=pd.read_csv('iso_id.txt',header=0,sep='\t')
iso_list=list(iso_id['iso ID'])

data1=iso_data[iso_data['Protein_ID'].isin(iso_list)]
data1.to_csv('isoform_seq.txt',index=False,sep='\t')





#提取fasta
seq=''
a=0
f=open('iso_sequence.fasta')
filename = 'isoform_seq.fasta'
with open(filename, 'a') as file:
	for line in f:
		if line[0]=='>' and seq=='':
			header=line
		elif line[0] !='>':
			seq=seq+line
		elif line[0] =='>' and seq !='':
			if header[header.index('|')+1:header.index('|',header.index('|')+1)] in iso_list:
				file.write(header+seq)
				seq=''
				header=line
				a=a+1
			else:
				seq=''
				header=line
f.close()
file.close()


#将isoform和canonical分别拆分成一个一个的fasta文件
from Bio import SeqIO
# 打开输入fasta文件
input_file = "isoform_seq.fasta"

# 读取输入fasta文件中的序列，并将每个序列写入单独的文件
for record in SeqIO.parse(input_file, "fasta"):
    output_file = "{}.fasta".format(record.id)
    with open(output_file, "w") as output_handle:
        SeqIO.write(record, output_handle, "fasta")
print("拆分完毕！")

########################################################

#使用clustalo做批量的序列比对
#clustalo --profile1 D3Z7P3-2.fasta --profile2 D3Z7P3.fasta --seqtype Protein --infmt fa --outfile result1.fasta
'''
O00159
O00187
clustalo --profile1 O00159-2.fasta --profile2 O00159.fasta --seqtype Protein --infmt fa --outfile result2.fasta
clustalo --profile1 O00187-2.fasta --profile2 O00187.fasta --seqtype Protein --infmt fa --outfile result3.fasta

clustalo --profile1 P11532-5.fasta --profile2 P11532.fasta --seqtype Protein --infmt fa --outfmt=clu --outfile result10
'''
import os
import pandas as pd
os.chdir('D://abishe//Fifth//Unreviewed//map//clustal-omega-1.2.2-win64')
path='D://abishe//Fifth//Unreviewed//map//clustal-omega-1.2.2-win64'
a=os.listdir(path)
code_list=[]
for i in range(len(a)):
	for j in range(len(a)):
		if 'fasta' in a[i] and 'fasta' in a[j]:
			if '-' in a[i] and a[i][:a[i].index('-')]==a[j][:a[j].index('.')]:
				code_i='clustalo --profile1 '+a[j]+' --profile2 '+a[i]+' --seqtype Protein --infmt fa --outfmt clu --outfile '+'result_'+a[i][:a[i].index('.')]+'.clu'
				code_list.append(code_i)
code_dat=pd.DataFrame({'code':code_list})
code_dat.to_csv('clustalo_code.txt',index=False,sep='\t')


'''
#验证
cuo_index=[]
for p in range(len(code_list)):
	temp1=code_list[p].split(' ')[2][:code_list[p].split(' ')[2].index('-')]
	temp2=code_list[p].split(' ')[4][:code_list[p].split(' ')[4].index('.')]
	if temp1 !=temp2:
		cuo_index.append(p)
'''


#获取isoform序列和canonical序列
import os
import pandas as pd
os.chdir('D://abishe//Fifth//Unreviewed//map')
a=pd.read_csv('Unreviewed_iso.tsv',header=0,sep='\t')
f=open('isoform_seq.fasta')
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

f=open('canonical_seq.fasta')
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

isos=[]
cans=[]
for i in range(len(a)):
	iso_id=a['Affected protein AC'][i][10:]
	iso_seq=list(all_isoform[all_isoform['isoform_ID']==iso_id]['isoform_seq'])[0]
	isos.append(iso_seq)
	can_id=a['Affected protein AC'][i][10:a['Affected protein AC'][i].index('-')]
	can_seq=list(Uniprot[Uniprot['Protein_ID']==can_id]['Protein_seq'])[0]
	cans.append(can_seq)
isos_dat=pd.DataFrame({'isoform seq':isos})
cans_dat=pd.DataFrame({'canonical seq':cans})
data=pd.concat([a,isos_dat,cans_dat],axis=1)
data.to_csv('Unreviewed_iso1.tsv',index=False,sep='\t')




#############################################################################################################################################

#学姐方法
#x是ID_pair
#将文件处理成c序列-i序列-map情况的dataframe，先保存下来

#去除dataframe某列等于某值并重新编行号
def function1(x,y,z):
    #x为dataframe,y为列名，z为包含需要去除的行对应的值的列表
    x = x[-x[y].isin(z)]
    return x.reset_index(drop=True)

def map_df(x):
    try:
        df2 = pd.read_csv(x,skip_blank_lines = False)
        #去掉nan的行
        a = df2
        b = 'CLUSTAL O(1.2.2) multiple sequence alignment'
        c = [NaN]
        df2 = function1(a,b,c)
        canonnical = []
        isoform = []
        map_result = []
        for i in range(int(len(df2['CLUSTAL O(1.2.2) multiple sequence alignment'])/3)):
            index1 = 3*i
            index2 = 3*i+1
            index3 = 3*i+2
            len1 = len(df2['CLUSTAL O(1.2.2) multiple sequence alignment'][0])-60
            len2 = len(df2['CLUSTAL O(1.2.2) multiple sequence alignment'][index1])-len1
            str1 = df2['CLUSTAL O(1.2.2) multiple sequence alignment'][index1][-len2:]
            canonnical.append(str1)
            str2 = df2['CLUSTAL O(1.2.2) multiple sequence alignment'][index2][-len2:]
            isoform.append(str2)
            str3 = df2['CLUSTAL O(1.2.2) multiple sequence alignment'][index3][-len2:]
            map_result.append(str3)
        canonnical_seq1 = ''.join(canonnical)
        isoform_seq1 = ''.join(isoform)
        map_result1 = ''.join(map_result)
        canonnical_seq = list(canonnical_seq1)
        isoform_seq = list(isoform_seq1)
        map_result = list(map_result1)
        df3 = pd.DataFrame()
        df3['canonnical_seq'] = canonnical_seq
        df3['isoform_seq'] = isoform_seq
        df3['map_result'] = map_result
        name = 'map1_'+x[:x.index('.')] +'.txt'
        df3.to_csv(name,index = 0, sep ='\t')
    except:
        return x

#为map上的序列残基编号
#x为ID_pair
def sort_seq(x):
    df_res = pd.read_csv(x,sep = '\t')
    amino_acid = 'ABCDEFGHIKLMNPQRSTUVWXYZ'
    amino_acid_list = list(amino_acid)
    number1 = 0
    number2 = 0
    num1 = []
    num2 = []
    res_pos1 = []
    res_pos2 = []
    for i in range(len(df_res['canonnical_seq'])):
        res1 = df_res['canonnical_seq'][i]
        if res1 in amino_acid_list:
            number1 = number1 + 1
            num1.append(str(number1))
            res_pos1.append(res1+str(number1))
        else:
            not_res = 'None'
            num1.append(not_res)
            res_pos1.append(not_res)
        res2 = df_res['isoform_seq'][i]
        if res2 in amino_acid_list:
            number2 = number2 + 1
            num2.append(str(number2))
            res_pos2.append(res2+str(number2))
        else:
            not_res = 'None'
            num2.append(not_res)
            res_pos2.append(not_res)
    df_res['canonnical_number'] = num1
    df_res['isoform_number'] = num2
    df_res['canonnical_res_pos'] = res_pos1
    df_res['isoform_res_pos'] = res_pos2
    df_res = df_res[['canonnical_seq','canonnical_number','canonnical_res_pos','isoform_seq','isoform_number','isoform_res_pos','map_result']]
    name = 'map2_'+x[5:x.index('.')]+'.txt'
    df_res.to_csv(name,index = 0, sep ='\t')
    return df_res



import os
import pandas as pd
from numpy import NaN
os.chdir('D://abishe//Fifth//Unreviewed//map//result1')
path='D://abishe//Fifth//Unreviewed//map//result1'
a=os.listdir(path)
for i in a:
	map_df(i)

#将结果文件转移到另一个文件夹result2中
os.chdir('D://abishe//Fifth//Unreviewed//map//result2')
path='D://abishe//Fifth//Unreviewed//map//result2'
b=os.listdir(path)
for j in b:
	sort_seq(j)

########################################################
#将序列比对结果保存在map_result文件夹中
#将突变位点进行替换
import os
import pandas as pd
os.chdir('D://abishe//Fifth//Unreviewed//map//map_result')
data=pd.read_csv('Unreviewed_iso1.tsv',header=0,sep='\t')
index1=[]
can_locs=[]
can_oris=[]
for i in range(len(data)):
	pro=data['Affected protein AC'][i][10:]
	loc=data['Feature range(s)'][i].split('-')[0]
	ori=data['Original sequence'][i]
	name='map2_result_'+pro+'.txt'
	a=pd.read_csv(name,header=0,sep='\t')
	a['isoform_number']=a['isoform_number'].astype(str)
	a['canonnical_number']=a['canonnical_number'].astype(str)
	can_loc=list(a[a['isoform_number']==loc]['canonnical_number'])[0]
	can_ori=list(a[a['isoform_number']==loc]['canonnical_seq'])[0]
	if ori != can_ori:
		index1.append(i)
	can_locs.append(can_loc)
	can_oris.append(can_ori)
data1=data.drop('Feature range(s)',axis=1)
locs_dat=pd.DataFrame({'Feature range(s)':can_locs})
data2=pd.concat([data1,locs_dat],axis=1)
data3=data2.drop('isoform seq',axis=1)
data3.to_csv('map_result.tsv',index=False,sep='\t')

'''
#验证
data.index=range(len(data))
index1=[]
for j in range(len(data)):
	loc=int(data['Feature range(s)'][j].split('-')[0])
	if data['isoform seq'][j][loc-1]!=data['Original sequence'][j]:
		index1.append(j)
'''

'''
#验证
a=pd.read_csv('Unreviewed_iso1.tsv',header=0,sep='\t')
index1=[]
for i in range(len(a)):
	a1=a['isoform seq'][i]
	a2=a['canonical seq'][i]
	if len(a1) > len(a2):
		index1.append(i)
#292个，isoform长度大于canonical
'''



'''
################################################################
#寻找字符在字符串中第几次出现的函数
def xunzhao(str1,strs,nums):
	for i in range(len(strs)):
		if strs[i]==str1 and strs[:i+1].count(str1)==nums:
			return i
#将突变位点map到canonical上
import os
import pandas as pd
os.chdir('D://abishe//Fifth//Unreviewed//map//result')
a=pd.read_csv('Unreviewed_iso1.tsv',header=0,sep='\t')
for j in range(len(a)):
	a['Affected protein AC'][j]=a['Affected protein AC'][j][10:]
a_seqs=[]
for i in range(len(a)):
	a_i=a['Affected protein AC'][i]
	file_name='result_'+a_i+'.fasta'
	f=open(file_name)
	seq={}
	for line in f:
		if line.startswith('>'):
			name=line.replace('\n','').split()[0]
			seq[name]=''
		else:
			seq[name]+=line.replace('\n','').strip()
	P_ID=[]
	for x in seq.keys():
		start=x.index('|')
		end=x.index('|',x.index('|')+1)
		P_ID.append(x[start+1:end])
	P_seq=[]
	for y in seq.values():
		P_seq.append(y)
	isoform_dic={'isoform_ID':P_ID,'isoform_seq':P_seq}
	all_isoform=pd.DataFrame(isoform_dic)
	f.close()
	can_name=all_isoform['isoform_ID'][1]
	a['Affected protein AC'][i]=can_name
	iso_loc=int(a['Feature range(s)'][i].split('-')[0])
	iso_ori=a['Original sequence'][i]
	iso_seq=a['isoform seq'][i]
	iso_num=iso_seq[:iso_loc].count(iso_ori)#第iso_num次出现iso_ori氨基酸
	can_seq=all_isoform['isoform_seq'][0]
	can_loc=xunzhao(iso_ori,can_seq,iso_num)+1
	a['Feature range(s)'][i]=can_loc
	a_seq=all_isoform['isoform_seq'][1]
	a_seqs.append(a_seq)
aseq_dat=pd.DataFrame({'A sequence':a_seqs})
data=pd.concat([a,aseq_dat],axis=1)
data1=data.drop('isoform seq',axis=1)
data1.to_csv('Unreviewed_map_can.tsv',index=False,sep='\t')



#验证，添加canonical的sequence序列
a=pd.read_csv('Unreviewed_map_can.tsv',header=0,sep='\t')

index3=[]
for j in range(len(a)):
	a1=a['Feature range(s)'][j]
	a4=a['A sequence'][j][:a1]
	if '-' in a4:
		a_num=a4.count('-')
		a_temp=a1-a_num
		a['Feature range(s)'][j]=a_temp

import os
import pandas as pd
os.chdir('D://abishe//Fifth//Unreviewed')
f=open('canonical_seq.fasta')
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
a1=a.drop('A sequence',axis=1)
seqs=[]
for p in range(len(a1)):
	temp_seq=list(Uniprot[Uniprot['Protein_ID']==a1['Affected protein AC'][p]]['Protein_seq'])[0]
	seqs.append(temp_seq)
seqs_dat=pd.DataFrame({'A sequence':seqs})
data=pd.concat([a1,seqs_dat],axis=1)
data.to_csv('map_result.tsv',index=False,sep='\t')
'''

#############################################################################################################
#添加另一条序列
data=pd.read_csv('map_result.tsv',header=0,sep='\t')
for i in range(len(data)):
	a1=data['Affected protein AC'][i][10:data['Affected protein AC'][i].index('-')]
	data['Affected protein AC'][i]=a1
temp_pro1=[]
for q in range(len(data)):
	temp_group=eval(data['Partners_1'][q])
	temp_pro=data['Affected protein AC'][q]
	if '-' in temp_group[0] and '-' not in temp_group[1]:
		if temp_group[0][10:temp_group[0].index('-')] ==temp_pro:
			temp_pro1.append(temp_group[1][10:])
		elif temp_group[1][10:] ==temp_pro:
			temp_pro1.append(temp_group[0][10:])

	elif '-' in temp_group[1] and '-' not in temp_group[0]:
		if temp_group[0][10:] ==temp_pro:
			temp_pro1.append(temp_group[1][10:])
		elif temp_group[1][10:temp_group[1].index('-')] ==temp_pro:
			temp_pro1.append(temp_group[0][10:])

	elif '-' not in temp_group[1] and '-' not in temp_group[0]:
		if temp_group[0][10:] ==temp_pro:
			temp_pro1.append(temp_group[1][10:])
		elif temp_group[1][10:] ==temp_pro:
			temp_pro1.append(temp_group[0][10:])

	elif '-' in temp_group[1] and '-' in temp_group[0]:
		if temp_group[0][10:temp_group[0].index('-')] ==temp_pro:
			temp_pro1.append(temp_group[1][10:])
		elif temp_group[1][10:temp_group[1].index('-')] ==temp_pro:
			temp_pro1.append(temp_group[0][10:])
pro1_dat=pd.DataFrame({'B ID':temp_pro1})
data2=pd.concat([data,pro1_dat],axis=1)
data2.to_csv('append_seq.tsv',index=False,sep='\t')
temp_pro2=list(set(temp_pro1))
pro2_dat=pd.DataFrame({'pro_B':temp_pro2})
pro2_dat.to_csv('B_pro_ID.txt',index=False,sep='\t')

#去uniprot上的ID-mapping，找到它们的序列，匹配

a=pd.read_csv('append_seq.tsv',header=0,sep='\t')
f=open('b_seq.fasta')
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
b_seq=[]
for i in range(len(a)):
	if '-' in a['B ID'][i]:
		temp_seq=list(Uniprot[Uniprot['Protein_ID']==a['B ID'][i][:a['B ID'][i].index('-')]]['Protein_seq'])[0]
		b_seq.append(temp_seq)
	else:
		temp_seq=list(Uniprot[Uniprot['Protein_ID']==a['B ID'][i]]['Protein_seq'])[0]
		b_seq.append(temp_seq)
bseq_dat=pd.DataFrame({'B sequence':b_seq})
data1=a.drop('B ID',axis=1)
data=pd.concat([data1,bseq_dat],axis=1)
data.to_csv('map_result1.tsv',index=False,sep='\t')

#保存有效数据

data1=pd.read_csv('map_result1.tsv',header=0,sep='\t')
data1.rename(columns={'canonical seq':'A sequence'},inplace=True)
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
g=data1['A sequence']
h=data1['B sequence']
data2=pd.concat([d,c,mut_data,e,f,b,g,h],axis=1)
data2.to_csv('Unreviewed_iso_result.tsv',index=False,sep='\t')




'''
a=pd.read_csv('map_result.tsv',header=0,sep='\t')
index2=[]
for i in range(len(a)):
	a1=a['Feature range(s)'][i]
	a2=a['Original sequence'][i]
	a3=a['A sequence'][i]
	if a3[a1-1] !=a2:
		index2.append(i)
'''


########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################



#检查突变位点与原序列一致性，正常数据的A和B的序列是不确定的
import os
os.chdir('D://abishe//Fifth//Unreviewed')
import pandas as pd
import re
data=pd.read_csv('Unreviewed_zc_result.tsv',header=0,sep='\t')
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
data_1.to_csv('Unreviewed_zc_cuowu.tsv',index=False,sep='\t')
data1.to_csv('Unreviewed_zc_zhengque.tsv',index=False,sep='\t')


#isoform数据，A序列为isoform匹配后的序列，B序列为剩余参与者B的canonical序列
import os
os.chdir('D://abishe//Fifth//Unreviewed')
import pandas as pd
import re
b=pd.read_csv('Unreviewed_iso_result.tsv',header=0,sep='\t')
index_cuo2=[]
index_cuo3=[]
for j in range(len(b)):
	temp2=b['A sequence'][j]
	if 'None' in b['Mutations'][j]:
		index_cuo3.append(j)
	else:
		b1=int(re.findall(r'\d+',b['Mutations'][j])[0])
		b2=b['Mutations'][j][0]
		temp2=b['A sequence'][j]
		if temp2[b1-1] !=b2:
			index_cuo2.append(j)
data_2=b.loc[index_cuo2]
data_3=b.loc[index_cuo3]
data_4=pd.concat([data_2,data_3])
data2=b.drop(index_cuo2)
data3=data2.drop(index_cuo3)
data_4.to_csv('Unreviewed_iso_cuowu.tsv',index=False,sep='\t')
data3.to_csv('Unreviewed_iso_zhengque.tsv',index=False,sep='\t')


import os
import pandas as pd
os.chdir('D://abishe//Fifth//Unreviewed')

a=pd.read_csv('Unreviewed_zc_zhengque.tsv',header=0,sep='\t')
b=pd.read_csv('Unreviewed_iso_zhengque.tsv',header=0,sep='\t')
data=pd.concat([a,b])
data.to_csv('Unreviewed_zhengque.tsv',index=False,sep='\t')

###########################################################################################

a=pd.read_csv('Unreviewed_zc_cuowu.tsv',header=0,sep='\t')
b=pd.read_csv('Unreviewed_iso_cuowu.tsv',header=0,sep='\t')
shuju_3=pd.concat([a,b])
shuju_3.to_csv('Unreviewed_cuowu.tsv',index=False,sep='\t')


'''
a=pd.read_csv('Unreviewed_zc_zhengque.tsv',header=0,sep='\t')
b=pd.read_csv('Unreviewed_iso_zhengque.tsv',header=0,sep='\t')
shuju_3=pd.concat([a,b])
print(len(shuju_3))
print(len(shuju_3.drop_duplicates(subset=['Affected protein AC','Mutations'],keep='first')))
print(len(shuju_3.drop_duplicates(subset=['Affected protein AC'],keep='first')))
print(len(shuju_3.drop_duplicates(subset=['Interaction AC'],keep='first')))
print(len(shuju_3.drop_duplicates(subset=['Partners_1'],keep='first')))
'''

#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#去除冲突数据

#Unreviewed
import os
os.chdir('D://abishe//Fifth//Unreviewed')
import pandas as pd

data=pd.read_csv('Unreviewed_zhengque.tsv',header=0,sep='\t')
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

data_result1.to_csv('Unreviewed_last1.tsv',index=False,sep='\t')

data_chuli.to_csv('Unreviewed_data_chuli.tsv',index=False,sep='\t')



a=pd.read_csv('Unreviewed_data_chuli.tsv',header=0,sep='\t')
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
data1=pd.read_csv('Unreviewed_last1.tsv',header=0,sep='\t')
data3=pd.concat([data1,data2])
data3.to_csv('Unreviewed_last_result.tsv',index=False,sep='\t')


'''
shuju_3=pd.read_csv('Unreviewed_last_result.tsv',header=0,sep='\t')
print(len(shuju_3))
print(len(shuju_3.drop_duplicates(subset=['Affected protein AC','Mutations'],keep='first')))
print(len(shuju_3.drop_duplicates(subset=['Affected protein AC'],keep='first')))
print(len(shuju_3.drop_duplicates(subset=['Interaction AC'],keep='first')))
print(len(shuju_3.drop_duplicates(subset=['Partners_1'],keep='first')))
'''



#按相互作用进行统计
#一对相互作用中，有多少个disrupting，有多少个decreasing

import os
import pandas as pd
os.chdir('D://abishe')
#Reviewed
a=pd.read_csv('Reviewed_last_result.tsv',header=0,sep='\t')

a1=a[a['type']=='Disrupting']
a2=a[a['type']=='Decreainsg']

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

a_data=a.loc[index_1]
a_data.to_csv('Reviewed_Disrupting_data.tsv',index=False,sep='\t')
a_partners=a_data['Partners_1']
a_int=a_data['Interaction AC']
a_1=pd.concat([a_int,a_partners],axis=1)
a_2=a_1.drop_duplicates(subset=['Interaction AC'],keep='first')
data2=pd.merge(data1,a_2,how='inner',on='Interaction AC')
data2.to_csv('Reviewed_Disrupting_nums.tsv',index=False,sep='\t')


#Decreasing
index2=[]
index_2=[]
index_a2=[]
a2_len=[]
a2_list=list(a2.groupby('Interaction AC'))
for i in range(len(a2_list)):
	len_num=len(a2_list[i][1])
	index2.append(len_num)
	index_2=index_2+list(a2_list[i][1].index)
	a2_len.append(a2_list[i][0])
data_int=pd.DataFrame({'Interaction AC':a2_len})
data_num=pd.DataFrame({'numbers':index2})
data1=pd.concat([data_int,data_num],axis=1)

a_data=a.loc[index_2]
a_data.to_csv('Reviewed_Decreasing_data.tsv',index=False,sep='\t')
a_partners=a_data['Partners_1']
a_int=a_data['Interaction AC']
a_2=pd.concat([a_int,a_partners],axis=1)
a_3=a_2.drop_duplicates(subset=['Interaction AC'],keep='first')
data2=pd.merge(data1,a_3,how='inner',on='Interaction AC')
data2.to_csv('Reviewed_Decreasing_nums.tsv',index=False,sep='\t')
##################################################################################################################


#Unreviewed
import os
import pandas as pd
os.chdir('D://abishe')
#Reviewed
a=pd.read_csv('Unreviewed_last_result.tsv',header=0,sep='\t')

a1=a[a['type']=='Disrupting']
a2=a[a['type']=='Decreasing']

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

a_data=a.loc[index_1]
a_data.to_csv('Unreviewed_Disrupting_data.tsv',index=False,sep='\t')
a_partners=a_data['Partners_1']
a_int=a_data['Interaction AC']
a_1=pd.concat([a_int,a_partners],axis=1)
a_2=a_1.drop_duplicates(subset=['Interaction AC'],keep='first')
data2=pd.merge(data1,a_2,how='inner',on='Interaction AC')
data2.to_csv('Unreviewed_Disrupting_nums.tsv',index=False,sep='\t')


#Decreasing
index2=[]
index_2=[]
index_a2=[]
a2_len=[]
a2_list=list(a2.groupby('Interaction AC'))
for i in range(len(a2_list)):
	len_num=len(a2_list[i][1])
	index2.append(len_num)
	index_2=index_2+list(a2_list[i][1].index)
	a2_len.append(a2_list[i][0])
data_int=pd.DataFrame({'Interaction AC':a2_len})
data_num=pd.DataFrame({'numbers':index2})
data1=pd.concat([data_int,data_num],axis=1)

a_data=a.loc[index_2]
a_data.to_csv('Unreviewed_Decreasing_data.tsv',index=False,sep='\t')
a_partners=a_data['Partners_1']
a_int=a_data['Interaction AC']
a_2=pd.concat([a_int,a_partners],axis=1)
a_3=a_2.drop_duplicates(subset=['Interaction AC'],keep='first')
data2=pd.merge(data1,a_3,how='inner',on='Interaction AC')
data2.to_csv('Unreviewed_Decreasing_nums.tsv',index=False,sep='\t')


#统计数据
shuju_3=pd.read_csv('Reviewed_Disrupting_data.tsv',header=0,sep='\t')

shuju_3=pd.read_csv('Reviewed_Decreasing_data.tsv',header=0,sep='\t')

shuju_3=pd.read_csv('Unreviewed_Disrupting_data.tsv',header=0,sep='\t')

shuju_3=pd.read_csv('Unreviewed_Decreasing_data.tsv',header=0,sep='\t')

print(len(shuju_3))
print(len(shuju_3.drop_duplicates(subset=['Affected protein AC','Mutations'],keep='first')))
print(len(shuju_3.drop_duplicates(subset=['Affected protein AC'],keep='first')))
print(len(shuju_3.drop_duplicates(subset=['Interaction AC'],keep='first')))
print(len(shuju_3.drop_duplicates(subset=['Partners_1'],keep='first')))








#处理序列信息数据
#Reviewed
import os
os.chdir('/data/liuyang/bishe_data/chuansong/tongji')
import pandas as pd

a=pd.read_csv('Reviewed_completed.tsv',header=0,sep='\t')
a1=a['Affected protein AC']
a2=a['b id']
a5=[]
a_seq=[]
for i in range(len(a)):
	temp_a1=a['Affected protein AC'][i]
	temp_a2=a['b id'][i]
	temp_a3=a['Mutations'][i]
	temp1=temp_a1+'-'+temp_a2+'_'+temp_a3+'_wt'
	a3=a['A sequence'][i]
	a4=a['B sequence'][i]
	temp2=a3+','+a4
	a5.append(temp1)
	a_seq.append(temp2)
a5_dat=pd.DataFrame({'complex id':a5})
a_dat=pd.DataFrame({'sequence':a_seq})
a6=pd.concat([a5_dat,a_dat],axis=1)
a6.to_csv('ref_Reviewed.txt',index=False,sep='\t')


def get_mutation(seq,loc,res):
	seq1=seq[:loc]+res+seq[loc+1:]
	return seq1

import re
index1=[]
for j in range(len(a)):
	temp_loc=int(re.findall(r'\d+',a['Mutations'][j])[0])
	temp_ref=a['Mutations'][j][0]
	temp_res=a['Mutations'][j][-1]
	temp_seq=a['A sequence'][j]
	temp_seq1=get_mutation(temp_seq,temp_loc-1,temp_res)
	a['A sequence'][j]=temp_seq1


a9=[]
b_seq=[]
for p in range(len(a)):
	temp_a1=a['Affected protein AC'][p]
	temp_a2=a['b id'][p]
	temp_a3=a['Mutations'][p]
	temp1=temp_a1+'-'+temp_a2+'_'+temp_a3+'_mut'
	a7=a['A sequence'][p]
	a8=a['B sequence'][p]
	temp2=a7+','+a8
	a9.append(temp1)
	b_seq.append(temp2)
a9_dat=pd.DataFrame({'complex id':a9})
b_dat=pd.DataFrame({'sequence':b_seq})
a10=pd.concat([a9_dat,b_dat],axis=1)
a10.to_csv('res_Reviewed.txt',index=False,sep='\t')
a_data=pd.concat([a6,a10])
a_data.to_csv('Reviewed_seqs.txt',index=False,sep='\t')


#Unreviewed
import os
os.chdir('/data/liuyang/bishe_data/chuansong/tongji')
import pandas as pd

a=pd.read_csv('Unreviewed_completed.tsv',header=0,sep='\t')
a1=a['Affected protein AC']
a2=a['b id']
a5=[]
a_seq=[]
for i in range(len(a)):
	temp_a1=a['Affected protein AC'][i]
	temp_a2=a['b id'][i]
	temp_a3=a['Mutations'][i]
	temp1=temp_a1+'-'+temp_a2+'_'+temp_a3+'_wt'
	a3=a['A sequence'][i]
	a4=a['B sequence'][i]
	temp2=a3+','+a4
	a5.append(temp1)
	a_seq.append(temp2)
a5_dat=pd.DataFrame({'complex id':a5})
a_dat=pd.DataFrame({'sequence':a_seq})
a6=pd.concat([a5_dat,a_dat],axis=1)
a6.to_csv('ref_Unreviewed.txt',index=False,sep='\t')


def get_mutation(seq,loc,res):
	seq1=seq[:loc]+res+seq[loc+1:]
	return seq1

import re
index1=[]
for j in range(len(a)):
	temp_loc=int(re.findall(r'\d+',a['Mutations'][j])[0])
	temp_ref=a['Mutations'][j][0]
	temp_res=a['Mutations'][j][-1]
	temp_seq=a['A sequence'][j]
	temp_seq1=get_mutation(temp_seq,temp_loc-1,temp_res)
	a['A sequence'][j]=temp_seq1


a9=[]
b_seq=[]
for p in range(len(a)):
	temp_a1=a['Affected protein AC'][p]
	temp_a2=a['b id'][p]
	temp_a3=a['Mutations'][p]
	temp1=temp_a1+'-'+temp_a2+'_'+temp_a3+'_mut'
	a7=a['A sequence'][p]
	a8=a['B sequence'][p]
	temp2=a7+','+a8
	a9.append(temp1)
	b_seq.append(temp2)
a9_dat=pd.DataFrame({'complex id':a9})
b_dat=pd.DataFrame({'sequence':b_seq})
a10=pd.concat([a9_dat,b_dat],axis=1)
a10.to_csv('res_Unreviewed.txt',index=False,sep='\t')
a_data=pd.concat([a6,a10])
a_data.to_csv('Unreviewed_seqs.txt',index=False,sep='\t')


###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
#分成一个的

#Reviewed
import os
os.chdir('/data/liuyang/bishe_data/chuansong/tongji')
import pandas as pd

a=pd.read_csv('Reviewed_completed.tsv',header=0,sep='\t')
a1=a['Affected protein AC']
a2=a['b id']
a3=a['A sequence']
a4=a['B sequence']

a5=pd.concat([a1,a3],axis=1)
a5.rename(columns={'A sequence': 'sequence'}, inplace=True)

a6=pd.concat([a2,a4],axis=1)
a6.rename(columns={'b id': 'Affected protein AC'}, inplace=True)
a6.rename(columns={'B sequence': 'sequence'}, inplace=True)

data_r=pd.concat([a5,a6])
data_r.to_csv('ref_Reviewed.txt',index=False,sep='\t')


def get_mutation(seq,loc,res):
	seq1=seq[:loc]+res+seq[loc+1:]
	return seq1

import re
index1=[]
for j in range(len(a)):
	temp_loc=int(re.findall(r'\d+',a['Mutations'][j])[0])
	temp_ref=a['Mutations'][j][0]
	temp_res=a['Mutations'][j][-1]
	temp_seq=a['A sequence'][j]
	temp_seq1=get_mutation(temp_seq,temp_loc-1,temp_res)
	a['A sequence'][j]=temp_seq1


for p in range(len(a)):
	temp_a1=a['Affected protein AC'][p]
	temp_a3=a['Mutations'][p]
	temp1=temp_a1+'_'+temp_a3+'_mut'
	a['Affected protein AC'][p]=temp1
a7=a['Affected protein AC']
a8=a['A sequence']
a9=pd.concat([a7,a8],axis=1)
a9.rename(columns={'A sequence': 'sequence'}, inplace=True)
a9.to_csv('res_Reviewed.txt',index=False,sep='\t')


b1=pd.read_csv('ref_Reviewed.txt',header=0,sep='\t')
b2=pd.read_csv('res_Reviewed.txt',header=0,sep='\t')
b3=pd.concat([b1,b2])
b4=b3.drop_duplicates(subset=['Affected protein AC'],keep='first')
b4.to_csv('Reviewed_seqs.txt',index=False,sep='\t')


#########################################################################################################################
#Unreviewed
import os
os.chdir('/data/liuyang/bishe_data/chuansong/tongji')
import pandas as pd

a=pd.read_csv('Unreviewed_completed.tsv',header=0,sep='\t')
a1=a['Affected protein AC']
a2=a['b id']
a3=a['A sequence']
a4=a['B sequence']

a5=pd.concat([a1,a3],axis=1)
a5.rename(columns={'A sequence': 'sequence'}, inplace=True)

a6=pd.concat([a2,a4],axis=1)
a6.rename(columns={'b id': 'Affected protein AC'}, inplace=True)
a6.rename(columns={'B sequence': 'sequence'}, inplace=True)

data_r=pd.concat([a5,a6])
data_r.to_csv('ref_Unreviewed.txt',index=False,sep='\t')


def get_mutation(seq,loc,res):
	seq1=seq[:loc]+res+seq[loc+1:]
	return seq1

import re
index1=[]
for j in range(len(a)):
	temp_loc=int(re.findall(r'\d+',a['Mutations'][j])[0])
	temp_ref=a['Mutations'][j][0]
	temp_res=a['Mutations'][j][-1]
	temp_seq=a['A sequence'][j]
	temp_seq1=get_mutation(temp_seq,temp_loc-1,temp_res)
	a['A sequence'][j]=temp_seq1


for p in range(len(a)):
	temp_a1=a['Affected protein AC'][p]
	temp_a3=a['Mutations'][p]
	temp1=temp_a1+'_'+temp_a3+'_mut'
	a['Affected protein AC'][p]=temp1
a7=a['Affected protein AC']
a8=a['A sequence']
a9=pd.concat([a7,a8],axis=1)
a9.rename(columns={'A sequence': 'sequence'}, inplace=True)
a9.to_csv('res_Unreviewed.txt',index=False,sep='\t')


b1=pd.read_csv('ref_Unreviewed.txt',header=0,sep='\t')
b2=pd.read_csv('res_Unreviewed.txt',header=0,sep='\t')
b3=pd.concat([b1,b2])
b4=b3.drop_duplicates(subset=['Affected protein AC'],keep='first')
b4.to_csv('Unreviewed_seqs.txt',index=False,sep='\t')