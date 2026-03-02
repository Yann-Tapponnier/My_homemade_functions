######## CREATING A LOOP FOR PLOTING THE HIGHLIGHTS BY ORIG.IDENT ############
LoopForPlotingHighlights <-  function (provide_object_names,
                                       Reductions = "umap",
                                       output_path, 
                                       Ncol = 2, 
                                       Width = 2000, 
                                       Height = 2000,
                                       metadata_to_highlight= "orig.ident",
                                       MultiplePlotOnSingleFile = TRUE,
                                       TitleSize = 20,
                                       GifGeneration = FALSE){
  
  
  # Package
  require(Seurat)
  require(patchwork)
  
  # Multiple Object loop
  for (y in 1:length(provide_object_names)){
    obj_name <- provide_object_names[y]        # This is just the string name, like "join_Harmony"
    Obj <- get(obj_name)          
    
    # Preparation for each object
    P <-  list()
    i <- 0
    CellToHiglight <- names(table(Obj@meta.data[[paste0(metadata_to_highlight)]]))
    
    
      # Loop
      for (name in CellToHiglight) {
        print(name)
        i <- i+1 
        
        P[[i]] <-  DimPlot(Obj, reduction = Reductions, cells.highlight = rownames(Obj@meta.data[Obj@meta.data[[metadata_to_highlight]] == name, ]))+
          ggtitle(print(name))+
          theme(
            plot.title = element_text(
              hjust = 0.5,
              size = TitleSize
            ))
        print(i)
      }
      
      
      # Checking/Creating the folder
      if (!dir.exists(output_path)) {
        dir.create(output_path, recursive = TRUE)
        message("Created folder: ", output_path)
      } else {
        message("Folder already exists: ", output_path)
      }
      
      if (MultiplePlotOnSingleFile){
      # Printing in Single Plot
      png(filename = paste0(output_path, "/Highlights_",obj_name,"_" ,metadata_to_highlight ,"_", Reductions, ".png") , width = Width, height = Height )
        print(wrap_plots(plots = P, ncol = Ncol ))
      dev.off() 
    
    } else {
      # Checking/Creating the folder Highlights
      output_path2 <- paste0(output_path,"/highlights")
        if (!dir.exists(output_path2)) {
          dir.create(output_path2, recursive = TRUE)
          message("Created folder: ", output_path2)
        } else {
          message("Folder already exists: ", output_path2)
        }
            
          
          i <- 0
          # Looping for plot
          for (name in CellToHiglight) {
            i <- i+1 
        
          # Printing in multiple Plots
          png(filename = paste0(output_path2, "/Highlights_",obj_name,"_" ,metadata_to_highlight ,"_",name,"_", Reductions, ".png") , width = Width, height = Height )
           print(P[[i]])
          dev.off()
          }
            if (GifGeneration){
              # Generating GIF with magick
              require(magick)

                imgs <- list.files(output_path2, full.names = TRUE, pattern = ".png$")
                # Remmettre dans l'ordre des jours
                imgs <- imgs[order(as.numeric(sub(".*_D([0-9]+)_.*", "\\1", imgs)))]
                
                img_list <- lapply(imgs, image_read)
                ## join the images together
                img_joined <- image_join(img_list)
                ## animate at 2 frames per second
                img_animated <- image_animate(img_joined, fps = 2)
                ## save to disk
                image_write(image = img_animated,
                            path = paste0(output_path2, "/Highlights_",metadata_to_highlight,"_animated.gif"))
            }
               
    }
  }
}
  


                                  ##### Exemple of use
                                  # LoopForPlotingHighlights(provide_object_names = "int_obj1",
                                  #                          metadata_to_highlight = "orig.ident",
                                  #                          output_path = "report/4_hiPAR_further_analysis/",
                                  #                          Width = 1600,
                                  #                          Height = 1200,
                                  #                          MultiplePlotOnSingleFile = FALSE,
                                  #                          TitleSize = 40,
                                  #                          GifGeneration = FALSE)
              


  ############################# Directory arborescence creation #######################
# Liste des dossiers à créer
FolderCreationR <-  function(){
  folders <- c(
  "data",
  "figures",
  "html_export",
  "output_seurat",
  "report",
  "signature",
  paste0("report/",tools::file_path_sans_ext(
    basename(rstudioapi::getActiveDocumentContext()$path)
  )),
  paste0("signature/",tools::file_path_sans_ext(
    basename(rstudioapi::getActiveDocumentContext()$path)
  ))
)

    # Création des dossiers s’ils n’existent pas déjà
    for (f in folders) {
      if (!dir.exists(f)) {
        dir.create(f, recursive = TRUE)
        message(paste("✅ Dossier créé :", f))
      } else {
        message(paste("ℹ️ Dossier déjà existant :", f))
      }
    }
}




###### MODIFIIYING THE FUNCTION OF EMILE TO REALLY USE FIND ALL MARKERS  ######
###### AddMarkers to seurat with FINDALLMARKER for REAL THIS TIME !!!
# remove the NA cluster temporarly
# Store the info for cerebro
# BE CAREFULL it is different WILCOXAUC and FIND ALL MARKERS /////      1 is sensitive  the other is specific : 1 pay more weight on foldchange / the other on % of cell expressing it and NOT in OUTclusters.
AddMarkerGenesToSeurat2_FAM <- function(seuratObj, 
                                   group_by, 
                                   filterOutNAcells=T ,
                                   P_val_adj = 0.05,
                                   pct_in_min=0, pct_out_min=0, 
                                   logFCposMin=0, logFCnegMax=-0){
  require(presto)
  require(dplyr)
  seuratObj_saved <- seuratObj
  
  # This is made to filter out temporalry the NA cluster cells (to not take it as a CLUSTER)
  if(filterOutNAcells){
    seuratObj <- subset(seuratObj, 
                        cells =  rownames(which(seuratObj@meta.data[,group_by, drop=F]!='NA', arr.ind = T))) # Remoiving false String NA
    seuratObj <- subset(seuratObj, 
                        cells =  rownames(which(!is.na(seuratObj@meta.data[,group_by, drop=F]), arr.ind = T))) # Remoiving actual NA
  }
  Idents(seuratObj) <- group_by
  seuratObj <- subset(seuratObj, downsample =100 )                         
  
  markers <- FindAllMarkers(seuratObj,
                            group.by = group_by, 
                            seurat_assay = 'RNA') %>% 
    dplyr::filter(p_val_adj < P_val_adj & (pct.1 >= pct_in_min | pct.2 >= pct_out_min ) & (avg_log2FC >= logFCposMin | avg_log2FC <= logFCnegMax)) %>%
    dplyr::select(-1) %>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = T) # .by_group = T     make descending order by group !!
  
  markers <- markers %>% dplyr::select (gene, cluster, avg_log2FC, p_val_adj, pct.1, pct.2) # rearanging the column order AND !!group_by, !! because group_by is a dynamical value
  colnames(markers) <- c("gene", group_by, "avg_log2FC", "p_val_adj", "pct.in", "pct.out")
  
  markers <- group_split(markers) 
  markers_df <- list()
  markers_df_pos <- list()
  markers_df_neg <- list()
  for (i in 1:length(markers)){ 
    markers_df[[i]] <- as.data.frame(markers[[i]])
    markers_df[[i]] <- markers_df[[i]] %>% mutate_at(.vars = group_by, .funs = as.factor)
    markers_df_pos[[i]] <- markers_df[[i]] %>% filter(avg_log2FC > 0)
    markers_df_neg[[i]] <- markers_df[[i]] %>% filter(avg_log2FC < 0)
  }
  
  #SAVE the output table in MISC to cerebro to read it !! 
  seuratObj_saved@misc$marker_genes$cerebro_seurat[[paste0(group_by,"_FAM")]] <- do.call(what = rbind, markers_df)
  seuratObj_saved@misc$marker_genes$cerebro_seurat[[paste0(group_by,"_FAM", '_upreg_only')]] <- do.call(what = rbind, markers_df_pos)
  seuratObj_saved@misc$marker_genes$cerebro_seurat[[paste0(group_by,"_FAM", '_downreg_only')]] <- do.call(what = rbind, markers_df_neg)
  return(seuratObj_saved)
}


########################################### make_unique_png ############################################
# I want to create a function to automatically check the if the name exist, otherwise it and _1 _2 incrementally.

make_unique_png <- function(path, filename) {
  full_path <- file.path(path, filename)
  
  if (!file.exists(full_path)) {
    return(full_path)
  }
  
  base <- sub("\\.png$", "", filename)
  i <- 1
  
  repeat {
    new_name <- paste0(base, "_", i, ".png")
    new_path <- file.path(path, new_name)
    if (!file.exists(new_path)) {
      return(new_path)
    }
    i <- i + 1
  }
}



###############################  Automation to plot AUC score cutted by quantiles ############################### 
PlotAUCscoreByQuantiles <- function(Seurat_Obj , #should be a Seurat Object
                                    Probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.99),
                                    Sig_names , #should be a character vector
                                    OutputFolderPath = getwd(),
                                    Obj_name = NULL,
                                    Width = 600,
                                    Height = 400,
                                    Res = 100) {
  
  require(patchwork)
  ###### Retreiving Dynamical name under of the real Seurat_obj ONLY when ploting directly OTHERWISE it uses the obj_name already caming from the function above
  if (is.null(Obj_name)) {
    Obj_name <- deparse(substitute(Seurat_Obj))
  }
  
  ##### Ploting with threshold
  for (sig in Sig_names) {
    # Define the quantiles and keep probs as names
    thresholds <- quantile(FetchData(Seurat_Obj, vars = sig)[, 1], probs = Probs)
    names(thresholds) <- Probs  # assign the probs as names
    
    message("Plotting all thresholds for ", sig)
    plots_list <- list()  # store plots for this signature
    
    # Safely extract the metadata column as numeric
    sig_vector <- FetchData(Seurat_Obj, vars = sig)[, 1]
    
    
    for (threshold_name in names(thresholds)) {
      threshold_value <- thresholds[threshold_name]  # the numeric value
      
      # Create a temporary highlight column
      highlight_col <- paste0(sig, "_highlight_", threshold_name)
      Seurat_Obj[[highlight_col]] <- ifelse(sig_vector > threshold_value, sig_vector, 0)
      
      # Create the FeaturePlot
      p <- FeaturePlot(Seurat_Obj, features = highlight_col) +
        ggtitle(paste0(sig, " (quantile ", threshold_name, ")")) +
        NoLegend()
      
      plots_list[[threshold_name]] <- p
    }
    
    # Combine all plots with patchwork
    ncol_layout <- ceiling(sqrt(length(thresholds)))
    combined_plot <- wrap_plots(plots_list, ncol = ncol_layout)
    
    # Save as a single PNG
    png(filename = paste0(OutputFolderPath,Obj_name,"_", sig, ".png"),
        width = ncol_layout*Width , height = ncol_layout*Height, res = Res)
    print(combined_plot)
    dev.off()
  }
}



# #### Example of use
# PlotAUCscoreByQuantiles(Seurat_Obj = int_obj1,
#                         Sig_names = names (List_Makers_upper),
#                         Probs = 0.75,
#                         OutputFolderPath = "report/4_hiPAR_further_analysis/",
#                         Width = 1600,
#                         Height = 1200,
#                         Res = 150)








##########################################  Automatically calculate and save the AUCscore of genesets ##########################################
AddSignatureAUCscoreSeuratObject <-  function (
    Seurat_Obj , # should be seurat obj
    Genesets , # should be named list of marker genes
    Output_path= paste0(getwd(),"/"), # need for last /
    Plot = F,
    Probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.99),
    verbose = TRUE
){
###### Info ##### 
  
  # Obj1 is the Seurat object itself
  if (!inherits(Seurat_Obj, "Seurat")) {
    stop("`Seurat_Obj` must be a Seurat object")
  }
  
  if (!is.list(Genesets)) {
    stop("`Genesets` should be named list of marker genes")
  }
  if (!grepl("/$", Output_path)) { 
    stop("Output_path should end with '/'. Please add a trailing slash.") 
  }
  
  if (verbose) print("So far the input seems OK")
  
  #### library
  if (verbose) print( "Loading libraries")
    require (Seurat)
    require (tidyverse)
    require (AUCell)
  if (verbose) print("Libraries loaded")
  
  ###### Retreiving Dynamical name under of the real Seurat_obj
  Obj_name <- deparse(substitute(Seurat_Obj))
  
  ###### 
  Counts <- GetAssayData(Seurat_Obj, assay = "RNA", layer= "counts") #AUCell is not afected if using COUNTS or DATA !!
  cells_rankings <-AUCell_buildRankings(Counts, nCores = 12, plotStats = TRUE)
  
  # Calculate enrichment scores
  if (verbose) print("Starting AUCscore actual computation")
  cells_AUC <- AUCell_run(exprMat = Counts,
                          geneSets = Genesets, 
                          aucMaxRank=nrow(cells_rankings)*0.05)
  
  
  
  # Plot histogram
    #Creating a unique name
  if (!exists("make_unique_png", mode = "function")) {
    source("/Users/administrateur/Desktop/Bio_info/My_codes/My_Functions_to_source.R")
  }
  
  # I need to remove the final "/" to not have 2 but keep it for the plot function
    Output_path_temp <- sub("/+$", "", Output_path) 
    
    png_unique_name <- make_unique_png(
        path = Output_path_temp,
        filename = paste0(Obj_name, "_AUCscore_histograms.png")
      )
    
  # Combine all plots with patchwork
  ncol_layout <- ceiling(sqrt(length(Genesets)))
  if (verbose) print("Ploting the AUCscore disctributions")
  
  png(filename = png_unique_name , width = ncol_layout*600, height = ncol_layout*450 )
    par(mfrow = c(ncol_layout, ncol_layout)) # to save it in xXx format
        for (i in seq(1:length(Genesets))){
         name <- names(Genesets)[i]
    
           hist(getAUC(cells_AUC)[name, ], main = paste0("AUC for ",name))
    
         }
  dev.off()
  par(mfrow = c(1, 1)) # to revert toward normal plotting  
  if (verbose) print("Ploting the AUCscore disctributions is over and it worked well ;-)")
  
  ### Check
  
  auc_matrix <- as.data.frame(t(getAUC(cells_AUC)))  # Transpose to match Seurat: cells = rows
  answer <- all(rownames(auc_matrix) %in% colnames(Seurat_Obj))  # Should be TRUE
  print(paste0( "testing is all rownames of AUC_MATRIX are the colnames of ", Obj_name,' and the answer is : ', answer))
  
  # Intergrating the result in Metadata
  Seurat_Obj <- AddMetaData(Seurat_Obj, metadata = auc_matrix)
  

  
    if (verbose) print("I you wanted to plot, it is ploting the AUCscores in the UMAP")
    ####### Automatically plot all the newly calculated   
  if (verbose) print("I have the impression the it will be legend... wait for it... ") 
    if (Plot == T) {
      
      PlotAUCscoreByQuantiles(Seurat_Obj = Seurat_Obj , #should be a character vector
                              Probs = Probs,
                              Sig_names = names(Genesets), #should be a character vector
                              OutputFolderPath = Output_path,
                              Obj_name = Obj_name) # need the Output_path containing the last "/" 
      
      if (verbose) print("Dary, Legendary !!")
      if (verbose) print("By the way, it is over, enjoy !")
    } else { if (verbose) print("but as you did not, it is over,
                                have a good one !")
     }
  
  # return the modified object to modify Seurat obj outide
  return(Seurat_Obj)  
  
}

                                            
                                              #### Exemle of use
                                              # ds_objf1 <- AddSignatureAUCscoreSeuratObject (Seurat_Obj = ds_objf1,
                                              #                                               Genesets = Genesets = list(deadend_up_FM_top100_FC = deadend_up_FM_top100_FC,
                                              #                                                                          deadend_up_wc_top100_FC = deadend_up_wc_top100_FC,
                                              #                                                                          deadend_up_wc_top100_AUC = deadend_up_wc_top100_AUC),
                                              #                                               Output_path = "report/8_RT_in_vivo_WO_ANSES_ds_Signatures_clean/" ,
                                              #                                               Plot = TRUE,
                                              #                                               verbose = T)





######################## Exporting TABLE of markers ########################
ExportingTablesOfMarkers <- function (markers, # Markers should be a vector of names
                                      output_path = getwd()
                                      ){
  for (marker in markers){
    write.table(get(marker), paste0(output_path,marker,".tsv"), sep = "\t", row.names = F)
  }
}




############################################################################## 
################################### plotly ################################### 
##############################################################################
PlotlyForLassoSelection  <- function(Seurat_obj,
                                     Reduction = "umap",
                                     Color.by = "orig.ident",
                                     Threshold = 0,
                                     Colors = "viridis") {
  
        ## 🔍 QC — vérifications des arguments
        if (!inherits(Seurat_obj, "Seurat")) {
          stop("❌ Seurat_obj must be a Seurat object")
        }
        
        if (!is.character(Color.by)) {
          stop("❌ Color.by must be a single character string (metadata name)")
        }
        
        if (!Color.by %in% colnames(Seurat_obj@meta.data)) {
          stop(paste0(
            "❌ Color.by = '", Color.by,
            "' not found in Seurat object metadata"
          ))
        }
        
        if (!is.numeric(Threshold) || length(Threshold) != 1 ||
            Threshold < 0 || Threshold > 1) {
          stop("❌ Threshold must be a single numeric value between 0 and 1")
        }
        # #0 Libraris and QC
        require(plotly)
        require(Seurat)
        require(shiny)
        
        
        
        # #1 Conversion of Seurat Object
        umap_df <- as.data.frame(Embeddings(Seurat_obj, reduction = Reduction))
        umap_df$cell <- colnames(Seurat_obj)
        if (Threshold==0){
          umap_df$cluster <- (Seurat_obj[[Color.by]]) %>% unlist()
        }else{
          umap_df$cluster <- (Seurat_obj[[Color.by]]) %>% unlist() 
          thr <- quantile(umap_df$cluster, Threshold, na.rm = TRUE)
          umap_df$cluster[umap_df$cluster <= thr] <- 0
        }
        
        # #2 — Créer l’interface Shiny (UI)  
        ui <- fillPage(
          h3("UMAP interactive (lasso)"), 
          plotlyOutput("umap_plot", height = "85vh"), # Height adjusted for fillPage of the umap plot. 85% for umap 15% for selected UMIs
          actionButton("save", "Sauvegarder la sélection"),
          verbatimTextOutput("info")
        )
        
        
        #3 — Créer le serveur Shiny
        server <- function(input, output, session) {
          ## 🔁 Sélection réactive unique
          selected_cells <- reactive({
            sel <- event_data("plotly_selected", source = "umap_source")
            if (is.null(sel)) return(NULL)
            unique(sel$key)
          })
          
          ## 1️⃣ Plot UMAP
          output$umap_plot <- renderPlotly({
            plot_ly(
              data = umap_df,
              x = ~umap_1,
              y = ~umap_2,
              key = ~cell,
              type = "scatter",
              mode = "markers",
              color = ~cluster,
              colors = Colors,
              marker = list(size = 5),
              source = "umap_source"
            ) %>%
              layout(dragmode = "lasso")
          })
          
          
          ## 2️⃣ Print des cellules sélectionnées
          output$info <- renderPrint({
            cells <- selected_cells()
            if (is.null(cells)) return("Aucune sélection")
            cells
          })
          
          
          ## 3️⃣ Nombre total de cellules sélectionnées
          output$n_cells <- renderText({
            cells <- selected_cells()
            if (is.null(cells)) return("0 cellule")
            paste(length(cells), "cellules sélectionnées")
            message("UMI récupérés pour ", length(cells), " cellules")
          })
          
          observeEvent(input$save, {
            assign(
              "cells_selected",
              selected_cells(),
              envir = .GlobalEnv
            )
          })
        }
        
        #4 — Lancer l’application Shiny
        shinyApp(ui = ui, server = server)
      }



              # #### Exemple of launch :
              # PlotlyForLassoSelection(
              #   Seurat_obj = int_obj1,
              #   Reduction = "umap",
              #   Color.by = "NK_Cl3_up_top100_FC",
              #   Threshold = 0.9)
              # 







############################################################################## 
################################### hdWGCNA ################################## 
##############################################################################
# Automation of ploting the regulon scores for a TF of interest in the UMAP with FeaturePlot and saving it in the right folder with the name of the TF and the day of interest (D7 here)

PlotingRegulon_Scores <- function(Seurat_obj, 
                                  Output_path = getwd(),
                                  Sufixe_name = "",
                                  Cur_tfs){
  
  for (Cur_tf in Cur_tfs){
    # select a TF of interest
    cur_tf <- Cur_tf
    
    # add the regulon scores to the Seurat metadata
    hd1$pos_regulon_score <- pos_regulon_scores[,cur_tf]
    hd1$neg_regulon_score <- neg_regulon_scores[,cur_tf]
    
    # plot using FeaturePlot
    p1 <- FeaturePlot(hd1, feature=cur_tf) + umap_theme()
    p2 <- FeaturePlot(hd1, feature='pos_regulon_score', cols=c('lightgrey', 'red')) + umap_theme()
    p3 <- FeaturePlot(hd1, feature='neg_regulon_score', cols=c('lightgrey', 'seagreen')) + umap_theme()
    
    png(filename = paste0(Output_path, cur_tf, Sufixe_name, ".png"), width = 2400, height = 900, res = 150)
    print(p1 | p2 | p3)
    dev.off()
    
    print(p1 | p2 | p3)
    
  }
}


# PlotingRegulon_Scores(Seurat_obj = hd1,
#                           Output_path = "report/7_hiPAR_hdWGCNA_D7_Transcription_Factor/hdWGCNA_Regulons_Scores/",
#                           Sufixe_name = "_D7",
#                           Cur_tf = "POU5F1")
# 









