### 0.0rawdata_QC
for i in 0.0raw_data/*.R1.fastq.gz;do
nohup fastqc $i -o 0.0rawdata_QC/ >> 0.0rawdata_QC/qc.log & 
done
multiqc .

### 0.1clean_data


for i in 0.0raw_data/*R1.fastq.gz;do 
o=${i##*/};echo $o
nohup cutadapt -j 40 -a AGATCGGAAGAGC -A AGATCGGAAGAGC --trim-n -m 50 -q 20,20 -o 0.1clean_data/$o -p 0.1clean_data/${o/R1/R2} $i ${i/R1/R2}> 0.1clean_data/logs/${o%%.R1*}.cutadapt.log & 
done

### 0.1cleandata_QC


for i in 0.1cleandata_QC/*R1.fastq.gz;do
nohup fastqc $i -o 0.1cleandata_QC/ >> 0.1cleandata_QC/qc.log &
done
multiqc .

### 1.0hicpro


nohup HiC-Pro -i data_for_hicpro/ -o 1.0hicpro -c config-hicpro.txt &
(data_for_hicpro 格式：各样本存放于单独文件夹中)


#### 02.validPairs_merge 合并结果文件

nohup HiC-Pro -i 01.bam_link/ -o 02.validPairs_merge -c config-hicpro.txt -s merge_persample &
（bam_link: hicpro运行结果hic_results中.validPairs文件）


#### 03.contact_map 计算互作matrix

nohup HiC-Pro -i 02.validPairs_merge/hic_results/data/ -o 03.contact_map -c config-hicpro.txt -s build_contact_maps &


#### 04.hic_resolution 计算理论上分辨率

####当我们用最小单位长度构建一个接触矩阵时，80%的单位长度所含的contact数量至少在1000以上，
####这种最小单位长度即为分辨率
for j in $raw_mat/*;do k=$(basename $j); for i in $raw_mat/${k}/raw/*/*.matrix;do 
o=${i##*/};o=${o%%.matrix};echo $o;
echo `wc -l $raw_mat/${k}/raw/*/${o}_abs.bed ` >> ${k}_bin_count.temp; 
awk 'BEGIN{FS="\t";OFS="\t"}($1!=$2){
A[$1]=A[$1]+$3;
A[$2]=A[$2]+$3}
END{for(i in A) print i,A[i]}' $i | sort -k2nr,2 | awk 'BEGIN{FS="\t";OFS="\t"} {n+=1} $2<1000{print "'$o'",n ;exit} 
END{print "'$o'",NR}' | uniq >> ${k}_stat.temp; done; done

for i in $raw_mat/*;do k=$(basename $i);echo $k;
paste ${k}_stat.temp ${k}_bin_count.temp | awk '{print $1,$2,$3-2,$2/($3-2)}' > ${k}_resolution.tab; done



#### 05.hicplotter


HiCPlotter=/usr/local/software/1.HiC-software/7.HiCPlotter-0.6.6/HiCPlotter.py
nohup python $HiCPlotter -f $matrix_dir/LT/raw/40000/LT_40000.matrix 
-bed $matrix_dir/LT/raw/40000/LT_40000_abs.bed -n LT_chr8 -o 05.hicplotter/40k/LT/LT_chr8 -chr chr8 -tri 1 -r 40000 -hmc 1 -ptr 1 &
#-f指定输入的Hi-C交互矩阵，-n指定交互矩阵热图的标题名称，-chr指定要画的染色质的名称，-r指定对应的分辨率，-o指定输出文件的前缀，-tri 1
声明输入的交互矩阵的格式是HiC-Pro的输出格式
#-hmc 热图的颜色 Greys(0), Reds(1), YellowToBlue(2), YellowToRed(3-default), Hot(4), BlueToRed(5).
#plotTriangular (-ptr): 绘制旋转半矩阵（下面的附加图，可以更方便地找到可能的TAD）
#-mm: an integer value for the interaction matrix heatmap scale upper-limit.
#可同时可视化多个交互 对应的-f和-n参数为空格分隔的多个参数

#加TAD
python $HiCPlotter -n YM_LT_Selp -tri 1 -bed $ConfigHP/40000_mm10.bed -r 40000 -f $iced_mat_Rmerge/LT/iced/40000/LT_40000_iced.matrix \
-pcd 1 -pcdf TAD/LT.tad -hist $di_bed/LT.bedGraph -hl Directional_Index -chr chr1 -s 4050 -e 4150 -hmc 1 -ptr 1 -o plot/YM_LT_Selp

#### 06.ice_norm


nohup HiC-Pro -i 03.contact_map/hic_results/matrix/ -o 06.ice_norm/ -c config-hicpro.txt -s ice_norm &

#### qc结果统计
比对结果统计
workdir=1.0hicpro/bowtie_results/bwt2/
for i in */*.pairstat;do o=${i##*/};o=${o%_*};echo $o >> a.temp;awk -F "\t" '$1~"Total_pairs_processed" {print $2/1000000 >> "b.temp"} \
$1~"Reported_pairs" {print $3/100 >> "c.temp"}' $i;done
paste a.temp b.temp c.temp > uniq_align_rate.tab
rm *temp

validpairs结果统计
workdir=1.0hicpro/hic_results/stats
for i in */*mergestat;do o=${i##*/};o=${o%_*};echo $o >> res.tab; \
awk -F "\t" 'NR <=4 {print $1,$2/1000000 >> "res.tab"}' $i;done

for i in */*mergestat;do o=${i##*/};o=${o%_*};echo $o; \
awk -F "\t" 'NR>1 {print $1,$2/1000000 >> "proportion.tab"}' $i;done
服务器hicpro升级 multiqc可用


### 2.0hicexplorer
#### 01.h5


#将hicPro互作matrix转为h5格式
for i in $raw_mat/*/raw/40000/*.matrix;do 
o=${i##*/};o=${o%%.*};echo $o; 
nohup hicConvertFormat --matrices $i --inputFormat hicpro --outFileName $h5_dir/${o}.h5 --outputFormat h5 --bedFileHicpro ${i/.matrix/_abs.bed} & done

#### 02.hicCorrelate

nohup hicCorrelate -m 01.h5/40k/*h5 --method=pearson --log1p --range 40000:2000000 --outFileNameHeatmap 02.hicCorrelate/40k-2M/Old_heatmap --outFileNameScatter 02.hicCorrelate/40k-2M/Old_scatterplot --plotFileFormat pdf &

#The –range parameter should be selected at a meaningful genomic scale according to, for example, the mean size of the TADs in the organism you work with.




