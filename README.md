# My_homemade_functions
I will deposit all my homemade functions that will generate progressively

YT2 :
This function aim to plot and save mFuzz plot and gene list using R mFuzz pakage.
mFuzz allow to find cluster though time eerie analysis.

I used it to find patterns in sorted list of cell lines

      # In more details YT2_Mfuzz:
      -   creates a working output directory including the dfname !  
      -   This function allow to import a df by its name.\
      -   lunch the mFuzz package depending on parameters we choose (=\> test with different c (number of clusters) and m (fuzzification parameter) parameters \#\#\#\#)\
      -   plot and save the PDF of cluster in the folder of our choise with the auto_generated name\
      -   it extract the geneset per each cluster\
      -   save the table with the same name in the same folder.\
          ?acore = storing all the score of belonging to each cluster\


# My_Function_to_source.R 
It contains all the functions I created to use in R in single-cell-RNAseq especially for Seurat object manipulation and plotting.
  # make_unique_png
  
  # LoopForPlotingHighlights
  
  # FolderCreationR
  
  # AddMarkerGenesToSeurat2_FAM
  
  # AddSignatureAUCscoreSeuratObject
  
  # PlotAUCscoreByQuantiles
