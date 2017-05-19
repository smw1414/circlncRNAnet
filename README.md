# Use of circlncnet in command mode

1. Download all scritps

2. Download the db files from(about 2GB)  
http://app.cgu.edu.tw/circlnc/db/db.zip  

3. Make all scipts executable   
``` chmod +x *.[Rr] ```  
and create a ouput folder  
``` mkdir output ```  

4. Perform differential expresion analysis  
``` ./deseq2.r ```  

5. Select the gene you intrerest then perform co-expression analysis  
```./correlation.r ```

6. Perfrom gene enrichment, RBP, miRNA sponge analysis  
  * lncRNA  
      + RBP ```./triple_network_lncRNA_RBP_2step.R```  
      + miRNA sponge ```./triple_network_lncRNA_sponge_2step.R`` `  
  * circRNA  
      + RBP ```./triple_network_circRNA_RBP_2step.R```  
      + miRNA sponge ```./triple_network_circRNA_sponge_2step.R```  

7. Scatter plot  
```./scatterplot.r```  

8. Co-expressed gene heatmap  
```./heatmap.R```


