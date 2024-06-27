# conda activate enrichment
# cd /share/home/shh/HCC_snRNA-seq/23.pathway/result
# python /share/home/shh/HCC_snRNA-seq/script/23.3.pathway.integration.py Hep02
# python /share/home/shh/HCC_snRNA-seq/script/23.3.pathway.integration.py Hep09
# python /share/home/shh/HCC_snRNA-seq/script/23.3.pathway.integration.py Hep2
# python /share/home/shh/HCC_snRNA-seq/script/23.3.pathway.integration.py Hep9


import math
import sys
################################################################################
##### target
args = sys.argv
choo_subtype = args[1]
# choo_subtype = 'Hep02'
# choo_subtype = 'Hep9'


################################################################################
##### setting
dir_dataset='/share/home/shh/HCC_snRNA-seq/23.pathway/dataset/'
dir_result='/share/home/shh/HCC_snRNA-seq/23.pathway/result/'
file_diffexp=open(dir_result+'result_marker_'+choo_subtype+'.txt','r')
file_genes=open(dir_result+'result_marker_'+choo_subtype+'.txt','r')
file_GO=open(dir_result+'result_GO_'+choo_subtype+'.txt','r')
file_KEGG=open(dir_result+'result_KEGG_'+choo_subtype+'.txt','r')
file_GSEA=open(dir_result+'result_GSEA_'+choo_subtype+'.txt','r')
save_path=dir_dataset

############################################################################
dict_ID={}
dict_exp={}

# file_genes=open(dir_dataset+'result_marker_'+choo_subtype+'.txt','r')  # debug
for line in file_genes:
  if line[0]!='#':
    info=line.replace("\n", "").split('\t')  
    dict_ID[info[0]]=info[6] 

file_genes.close()

# file_diffexp=open(dir_dataset+'result_marker_'+choo_subtype+'.txt','r')  # debug
for line in file_diffexp:
  if line[0]!='#':
    info=line.replace("\n", "").split('\t')  
    SYMBOL=info[0]
    Log2FC=float(info[2])
    if SYMBOL in dict_ID.keys():
      ENTREZID=dict_ID[SYMBOL]
      if Log2FC > 0:
        dict_exp[ENTREZID]='up'
      else:
        dict_exp[ENTREZID]='down'

file_diffexp.close()



################################################################################
##### KEGG
save_KEGG=open(save_path+'KEGG_'+choo_subtype+'.txt','w')
save_KEGG.write('id'+'\t'+'category'+'\t'+'gene_num.min'+'\t'+'gene_num.max'+'\t'+'gene_num.rich'+'\t'+'-log10Pvalue'+'\t'+'up.regulated'+'\t'+'down.regulated'+'\t'+'rich.factor'+'\n')

# file_KEGG=open(dir_result+'result_KEGG_'+choo_subtype+'.txt','r')  # debug
for line in file_KEGG:
  if line[0] != '#':
    info=line.replace("\n", "").split('\t')  
    id=info[2]
    category=info[0]
    list_counts_max=info[4].split('/')
    counts=int(list_counts_max[0])
    max=int(list_counts_max[1])
    pvalue=-math.log10(float(info[6]))
    num_up=0
    num_down=0
    list_genes=info[-2].split('/')
    for gene in list_genes:
      if gene in dict_exp.keys():
        regulate=dict_exp[gene]
        if regulate == 'up':
          num_up+=1
        elif regulate == 'down':
          num_down+=1
    ratio=counts/max
    save_KEGG.write(id+'\t'+category+'\t'+'0'+'\t'+str(max)+'\t'+str(counts)+'\t'+str(pvalue)+'\t'+str(num_up)+'\t'+str(num_down)+'\t'+str(ratio)+'\n')


file_KEGG.close()
save_KEGG.close()


################################################################################
##### GSEA
save_GSEA=open(save_path+'GSEA_'+choo_subtype+'.txt','w')
save_GSEA.write('id'+'\t'+'category'+'\t'+'gene_num.min'+'\t'+'gene_num.max'+'\t'+'gene_num.rich'+'\t'+'-log10Pvalue'+'\t'+'up.regulated'+'\t'+'down.regulated'+'\t'+'rich.factor'+'\n')

# file_GSEA=open(dir_result+'result_GSEA_'+choo_subtype+'.txt','r')  # debug
for line in file_GSEA:
  if line[0] != '#':
    info=line.replace("\n", "").split('\t')  
    id=info[0]
    category=info[0]
    max=info[2]
    pvalue=-math.log10(float(info[5]))
    ratio=info[3]
    counts=0
    num_up=0
    num_down=0
    list_genes=info[-1].split('/')
    # print(len(list_genes))
    for gene in list_genes:
      counts+=1
      if gene in dict_exp.keys():
        regulate=dict_exp[gene]
        # print(regulate)
        if regulate == 'up':
          num_up+=1
        if regulate == 'down':
          num_down+=1
      print(num_down)
    save_GSEA.write(id+'\t'+category+'\t'+'0'+'\t'+str(max)+'\t'+str(counts)+'\t'+str(pvalue)+'\t'+str(num_up)+'\t'+str(num_down)+'\t'+str(ratio)+'\n')

file_GSEA.close()
save_GSEA.close()
print('finished!')



