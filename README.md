# Use of circlncRNAnet in command line mode

1. Download all scritps

2. Download db files from (about 2GB)  
http://app.cgu.edu.tw/circlnc/db/db.zip  

3. Make all scripts executable   
``` chmod +x *.[Rr] ```  
and create an output folder  
``` mkdir output ```  

4. Prepare your own gene matrix file or download demo files from http://app.cgu.edu.tw/circlnc/  

5. Perform differential expresion analysis  
``` ./deseq2.r  ```  

6. Select the gene you intrerested, then perform co-expression analysis  
```./correlation.r ```

7. Perfrom gene enrichment, RBP and miRNA sponge analysis  
  * lncRNA  
      + RBP ```./triple_network_lncRNA_RBP_2step.R```  
      + miRNA sponge ```./triple_network_lncRNA_sponge_2step.R```  
  * circRNA  
      + RBP ```./triple_network_circRNA_RBP_2step.R```  
      + miRNA sponge ```./triple_network_circRNA_sponge_2step.R```  

8. Scatter plot  
```./scatterplot.r```  

9. Co-expressed gene heatmap  
```./heatmap.R```


