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
      Ncol
      if (MultiplePlotOnSingleFile){
      # Printing in Single Plot
      png(filename = paste0(output_path, "/Highlights_",obj_name,"_" ,metadata_to_highlight ,"_", Reductions, ".png") , width = Width, height = Height )
        print(wrap_plots(plots = P, ncol =  ))
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
                imgs <- imgs[order(as.numeric(sub(".*_[A-Z]([0-9]+)_.*", "\\1", imgs)))]
                
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
  


                                  # ##### Exemple of use
                                  # LoopForPlotingHighlights (provide_object_names = vector,
                                  #                                        Reductions = "umap",
                                  #                                        output_path, 
                                  #                                        Ncol = 2, 
                                  #                                        Width = 2000, 
                                  #                                        Height = 2000,
                                  #                                        metadata_to_highlight= "orig.ident",
                                  #                                        MultiplePlotOnSingleFile = TRUE,
                                  #                                        TitleSize = 20,
                                  #                                        GifGeneration = FALSE)
              


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
PlotAUCscoreByQuantiles <- function(Seurat_Obj,
                                    Probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.99),
                                    Sig_names,
                                    OutputFolderPath = getwd(),
                                    Obj_name = NULL,
                                    Width = 600,
                                    Height = 400,
                                    Res = 100,
                                    Wrapped_plot = FALSE,
                                    Suffixe = NULL) { # Should be the Sig_names general name when wrapped_plot is used ex "Scheibinger_signatures"
  
  require(patchwork)
  if (!dir.exists(OutputFolderPath)) dir.create(OutputFolderPath, recursive = TRUE)
  
  if (is.null(Obj_name)) {
    Obj_name <- deparse(substitute(Seurat_Obj))
  }
  
  # ── WRAPPED : on accumule tous les plots de toutes les signatures ──
  if (Wrapped_plot) {
    all_plots <- list()
  }
  
  for (sig in Sig_names) {
    
    thresholds <- quantile(FetchData(Seurat_Obj, vars = sig)[, 1], probs = Probs)
    names(thresholds) <- Probs
    message("Plotting all thresholds for ", sig)
    
    sig_vector <- FetchData(Seurat_Obj, vars = sig)[, 1]
    plots_list <- list()
    
    for (threshold_name in names(thresholds)) {
      threshold_value <- thresholds[threshold_name]
      
      highlight_col <- paste0(sig, "_highlight_", threshold_name)
      Seurat_Obj[[highlight_col]] <- ifelse(sig_vector > threshold_value, sig_vector, 0)
      
      p <- FeaturePlot(Seurat_Obj, features = highlight_col) +
        ggtitle(paste0(sig, " (quantile ", threshold_name, ")")) +
        NoLegend()
      
      plots_list[[threshold_name]] <- p
    }
    
    if (Wrapped_plot) {
      # ✅ Accumule les plots de cette signature dans la liste globale
      all_plots <- c(all_plots, plots_list)
      
    } else {
      # ✅ Sauvegarde un PNG par signature
      ncol_layout <- ceiling(sqrt(length(thresholds)))
      combined_plot <- wrap_plots(plots_list, ncol = ncol_layout)
      
      png(filename = paste0(OutputFolderPath, Obj_name, "_", sig, ".png"),
          width = ncol_layout * Width,
          height = ncol_layout * Height,
          res = Res)
      print(combined_plot)
      dev.off()
    }
  }
  
  # ✅ Wrapped : sauvegarde un seul PNG avec toutes les signatures
  if (Wrapped_plot) {
    ncol_layout <- ceiling(sqrt(length(all_plots)))
    combined_plot <- wrap_plots(all_plots, ncol = ncol_layout)
    
    png(filename = paste0(OutputFolderPath, Obj_name, "_all_signatures_",Suffixe,".png"),
        width = ncol_layout * Width,
        height = ncol_layout * Height,
        res = Res)
    print(combined_plot)
    dev.off()
    
    message("Saved wrapped plot: ", Obj_name, "_all_signatures.png")
  }
}



# #### Example of use
# PlotAUCscoreByQuantiles(Seurat_Obj = int_obj1,
#                         Sig_names = names (List_Makers_upper),












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


######################## Exporting TABLE of markers from Seurat object ########################
ExportingTablesOfMarkersFromSeuratObj <- function(Seurat_obj,
                                                  Vec_markers = NULL,
                                                  output_path = getwd()) {
  
  # ✅ Crée le dossier si il n'existe pas (récursif)
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
    message("Created output directory: ", output_path)
  } else {
    message("Output directory already exists: ", output_path)
  }
  
  if (is.null(Vec_markers)) {
    Markers_lists <- names(Seurat_obj@misc$marker_genes$cerebro_seurat)[
      !grepl("only", names(Seurat_obj@misc$marker_genes$cerebro_seurat))
    ]
  } else {
    Markers_lists <- Vec_markers
  }
  
  for (Markers in Markers_lists) {
    df <- Seurat_obj@misc$marker_genes$cerebro_seurat[[Markers]]
    write.table(df, paste0(output_path, "/", Markers, ".tsv"), sep = "\t", row.names = FALSE)
  }
}
                              # Exemple of use :
                              # ExportingTablesOfMarkersFromSeuratObj (Seurat_obj = int_obj2,
                              #                                        Vec_markers = NULL, # if NULL it export every marker list in Seurat_obj@misc$marker_genes$cerebro_seurat
                              #                                        output_path = "signature/5_NBEI1_integration_clustering_DOX_PlusE14")  # --> Can create folders



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
  
  ####### QCs ############
  if (!inherits(Seurat_obj, "Seurat")) {
    stop("Seurat_obj is not a Seurat object.")
  }
  message("✓ Seurat object detected")
  
  
  ########## Validation of the presence of the regulon scores in the metadata ############
  # TFs present in both matrices
  common_tfs <- intersect(colnames(pos_regulon_scores),
                          colnames(neg_regulon_scores))
  
  # Keep only requested TFs that are valid
  sub_cur_tfs <- intersect(Cur_tfs, common_tfs)
  
            # Report missing TFs
            removed_tfs <- setdiff(Cur_tfs, sub_cur_tfs)
            
            if (length(removed_tfs) > 0) {
              message("⚠ Removed TFs not present in BOTH regulon matrices: ",
                      paste(removed_tfs, collapse = ", "))
            }
            
            if (length(sub_cur_tfs) == 0) {
              stop("None of the requested TFs are present in both regulon matrices.")
            }
  
  message("✓ Number of TFs retained for plotting: ", length(sub_cur_tfs))
   
  # retreiving the name of the Seurat object for naming the output files         
  seurat_obj_name <- deparse(substitute(Seurat_obj))
            
  ###### Ploting ###########
  for (Cur_tf in sub_cur_tfs){
    # select a TF of interest
    cur_tf <- Cur_tf
    paste0(Output_path, cur_tf, Sufixe_name, ".png")
    # 
    # add the regulon scores to the Seurat metadata
    Seurat_obj$pos_regulon_score <- pos_regulon_scores[,cur_tf]
    Seurat_obj$neg_regulon_score <- neg_regulon_scores[,cur_tf]
    
    # plot using FeaturePlot
    p1 <- FeaturePlot(Seurat_obj, feature=cur_tf) + umap_theme()
    p2 <- FeaturePlot(Seurat_obj, feature='pos_regulon_score', cols=c('lightgrey', 'red')) + umap_theme()
    p3 <- FeaturePlot(Seurat_obj, feature='neg_regulon_score', cols=c('lightgrey', 'seagreen')) + umap_theme()
    
    png(filename = paste0(Output_path,seurat_obj_name,"_", cur_tf, Sufixe_name, ".png"), width = 2400, height = 900, res = 150)
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





############################################################################################################################################################
# =============================================================================
# Leiden clustering — sweep de résolutions + UMAP
# =============================================================================

# ── Paramètres ────────────────────────────────────────────────────────────────
FindNeighborsLeidenMultiRes <- function (Seurat_obj,
                                         K = 20 , # voisins pour le KNN / SNN graph
                                         RESOLUTIONS = 0.5,    # résolutions à tester
                                         DIMS = 1:30, 
                                         OutputPath = "Report/") {   # dimensions PCA utilisées
  # Name of the object
  Name_Obj <- deparse(substitute(Seurat_obj))
  
  # ── 1. KNN / SNN graph ───────────────────────────────────────────────────────
  
  Seurat_obj <- FindNeighbors(Seurat_obj, dims = DIMS, k.param = K, verbose = FALSE)
  
  # ── 2. Leiden à chaque résolution ────────────────────────────────────────────
  
  # FindClusters avec algorithm = 4 → Leiden natif dans Seurat
  for (res in RESOLUTIONS) {
    Seurat_obj <- FindClusters(
      Seurat_obj,
      algorithm        = 4,       # 4 = Leiden
      resolution       = res,
      random.seed      = 42,
      verbose          = FALSE
    )
    # Renommer la colonne pour garder chaque résolution
    col_name <- paste0("leiden_", res)
    Seurat_obj[[col_name]] <- Seurat_obj$seurat_clusters
  }
  
  # ── 3. UMAP — un panel par résolution ────────────────────────────────────────
  plots <- lapply(RESOLUTIONS, function(res) {
    col_name  <- paste0("leiden_", res)
    n_clusters <- nlevels(Seurat_obj[[col_name]][, 1])
    
    DimPlot(Seurat_obj, group.by = col_name, label = TRUE, label.size = 3,
            repel = TRUE, pt.size = 0.3) +
      labs(title = paste0("res = ", res, "  (", n_clusters, " clusters)")) +
      theme_classic(base_size = 10) +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold", size = 10))
  })
  
  p_final <- wrap_plots(plots, ncol = ceiling(sqrt(length(RESOLUTIONS)))) +
    plot_annotation(title    = paste0("Leiden — sweep de résolutions  (K = ", K, ")"),
                    theme    = theme(plot.title = element_text(face = "bold", size = 13)))
  
  ggsave(paste0(OutputPath,Name_Obj,"_leiden_umap_sweep.pdf"), p_final, width = 5 * ceiling(sqrt(length(RESOLUTIONS))),
         height = 5 * ceiling(length(RESOLUTIONS) / ceiling(sqrt(length(RESOLUTIONS)))))
  
  message("✓ Fichier sauvegardé : leiden_umap_sweep.pdf")
  return(Seurat_obj)
}



              # exemple :
              #   FindNeighborsLeidenMultiRes (Seurat_obj,
              #                                K = 20 , # voisins pour le KNN / SNN graph
              #                                RESOLUTIONS = seq(0.1,1,0.1),    # résolutions à tester
              #                                DIMS = 1:30, 
              #                                OutputPath = "Report/3_RA_integration_clustering") 
              # 


#################### Basic function to convert human to mouse gene names ########################
# BiocManager::install("orthogene")

# method = "homologene"  # ✅ local, ultra-rapide, coverage correct
# method = "gprofiler"   # requête réseau mais meilleure coverage
# method = "babelgene"   # autre base locale, bonne alternative

# Basic function to convert human to mouse gene names
convertHumanToMouseGeneList  <- function(x) {
  require("orthogene")
  
  human_genes <- orthogene::convert_orthologs(
    gene_df        = x,
    input_species  = "human",
    output_species = "mouse",
    method         = "homologene"
  )
  
  humanx <- rownames(human_genes)
  print(head(humanx))
  return(humanx)
}

# Basic function to convert mouse to human gene names
convertMouseToHumanGeneList <- function(x) {
  require("orthogene")
  
  human_genes <- orthogene::convert_orthologs(
    gene_df        = x,
    input_species  = "mouse",
    output_species = "human",
    method         = "homologene"
  )
  
  humanx <- rownames(human_genes)
  print(head(humanx))
  return(humanx)
}

                                        # Exemple Humangenes <- convertMouseToHumanGeneList(MouseGenes)








#####################  PlotingRepartitionByCluster
PlotingRepartitionByCluster <- function(
    Seurat_Obj,
    Clustering,
    outputpath = ".",
    plot_bar = TRUE,
    plot_box = FALSE,
    plot_stacked = TRUE
){
  
  library(ggplot2)
  library(patchwork)
  Obj_name <- deparse(substitute(Seurat_Obj))
  for (Clusters in Clustering){
    
    barplot_list <- list()
    boxplot_list <- list()
    
    ## =========================
    ## TABLE: orig.ident × cluster
    ## =========================
    prop_table <- prop.table(
      table(
        Seurat_Obj$orig.ident,
        Seurat_Obj@meta.data[[Clusters]]
      ),
      margin = 2
    )
    
    prop_df <- as.data.frame(prop_table)
    colnames(prop_df) <- c("orig.ident", "cluster", "proportion")
    
    ## =========================================================
    ## BAR + BOX PER CLUSTER (same as before)
    ## =========================================================
    for (cl in unique(prop_df$cluster)){
      
      tmp <- subset(prop_df, cluster == cl)
      
      ## ---------------- BARPLOT ----------------
      if (plot_bar){
        
        p_bar <- ggplot(tmp, aes(x = orig.ident, y = proportion)) +
          geom_bar(stat = "identity") +
          ggtitle(paste0(Clusters, " - Cluster ", cl)) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        barplot_list[[paste0("cluster_", cl)]] <- p_bar
      }
      
      ## ---------------- BOXPLOT ----------------
      if (plot_box){
        
        p_box <- ggplot(tmp, aes(x = orig.ident, y = proportion)) +
          geom_boxplot() +
          ggtitle(paste0(Clusters, " - Cluster ", cl)) +
          theme_bw()
        
        boxplot_list[[paste0("cluster_", cl)]] <- p_box
      }
    }
    
    ## =========================================================
    ## STACKED BARPLOT (REPLACES HISTOGRAM)
    ## =========================================================
    if (plot_stacked){
      
      ## IMPORTANT: normalize by orig.ident for stacked composition
      prop_table2 <- prop.table(
        table(
          Seurat_Obj$orig.ident,
          Seurat_Obj@meta.data[[Clusters]]
        ),
        margin = 1
      )
      
      stacked_df <- as.data.frame(prop_table2)
      colnames(stacked_df) <- c("orig.ident", "cluster", "proportion")
      
      p_stacked <- ggplot(
        stacked_df,
        aes(
          x = orig.ident,
          y = proportion,
          fill = cluster
        )
      ) +
        geom_bar(stat = "identity") +
        ggtitle(paste0("Cluster composition per orig.ident - ", Clusters)) +
        ylab("Relative frequency") +
        xlab("orig.ident") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      outfile_stacked <- file.path(
        outputpath,
        paste0(Obj_name,"_Repartition_", Clusters, "_STACKED.png")
      )
      
      png(outfile_stacked, width = 1200, height = 800, res = 150)
      print(p_stacked)
      dev.off()
    }
    
    ## =========================================================
    ## SAVE BARPLOTS
    ## =========================================================
    if (plot_bar && length(barplot_list) > 0){
      
      nplots <- length(barplot_list)
      ncol <- ceiling(sqrt(nplots))
      
      png(
        file.path(outputpath, paste0(Obj_name,"_Repartition_", Clusters, "_BARPLOTS.png")),
        width = 400 * ncol,
        height = 400 * ceiling(nplots / ncol),
        res = 100
      )
      
      print(wrap_plots(barplot_list, ncol = ncol))
      dev.off()
    }
    
    ## =========================================================
    ## SAVE BOXPLOTS
    ## =========================================================
    if (plot_box && length(boxplot_list) > 0){
      
      nplots <- length(boxplot_list)
      ncol <- ceiling(sqrt(nplots))
      
      png(
        file.path(outputpath, paste0(Obj_name,"_Repartition_", Clusters, "_BOXPLOTS.png")),
        width = 400 * ncol,
        height = 400 * ceiling(nplots / ncol),
        res = 100
      )
      
      print(wrap_plots(boxplot_list, ncol = ncol))
      dev.off()
    }
  }
}
#                         Probs = 0.75,
#                         OutputFolderPath = "report/4_hiPAR_further_analysis/",
#                         Width = 1600,
#                         Height = 1200,
#                         Res = 150)





                  ########## AJOUT DE COULEUR SPECIFIC À RETRAVAILLER ET INTÉGRÉ
                        # if (plot_stacked){
                        #   
                        #   ## IMPORTANT: normalize by orig.ident for stacked composition
                        #   prop_table2 <- prop.table(
                        #     table(
                        #       Seurat_Obj$orig.ident,
                        #       Seurat_Obj@meta.data[[Clusters]]
                        #     ),
                        #     margin = 1
                        #   )
                        #   
                        #   stacked_df <- as.data.frame(prop_table2)
                        #   colnames(stacked_df) <- c("orig.ident", "cluster", "proportion")
                        #   
                        #   ## --- Palette personnalisée ---
                        #   # named vector : les noms DOIVENT correspondre aux valeurs de "cluster"
                        #   # (donc aux levels de Seurat_Obj@meta.data[[Clusters]])
                        #   if (!is.null(cluster_colors)) {
                        #     
                        #     # on vérifie que tous les clusters présents ont une couleur définie
                        #     missing_clusters <- setdiff(unique(stacked_df$cluster), names(cluster_colors))
                        #     if (length(missing_clusters) > 0) {
                        #       warning(paste0(
                        #         "Pas de couleur définie pour : ",
                        #         paste(missing_clusters, collapse = ", "),
                        #         " -> attribution de gris par défaut."
                        #       ))
                        #       extra_colors <- setNames(
                        #         rep("grey80", length(missing_clusters)),
                        #         missing_clusters
                        #       )
                        #       cluster_colors <- c(cluster_colors, extra_colors)
                        #     }
                        #     
                        #     # on force l'ordre des facteurs selon l'ordre du named vector (optionnel mais pratique)
                        #     stacked_df$cluster <- factor(stacked_df$cluster, levels = names(cluster_colors))
                        #   }
                        #   
                        #   p_stacked <- ggplot(
                        #     stacked_df,
                        #     aes(
                        #       x = orig.ident,
                        #       y = proportion,
                        #       fill = cluster
                        #     )
                        #   ) +
                        #     geom_bar(stat = "identity") +
                        #     ggtitle(paste0("Cluster composition per orig.ident - ", Clusters)) +
                        #     ylab("Relative frequency") +
                        #     xlab("orig.ident") +
                        #     theme_bw() +
                        #     theme(axis.text.x = element_text(angle = 45, hjust = 1))
                        #   
                        #   if (!is.null(cluster_colors)) {
                        #     p_stacked <- p_stacked + scale_fill_manual(values = cluster_colors)
                        #   }
                        #   
                        #   outfile_stacked <- file.path(
                        #     outputpath,
                        #     paste0(Obj_name,"_Repartition_", Clusters, "_STACKED.png")
                        #   )
                        #   
                        #   png(outfile_stacked, width = 1200, height = 800, res = 150)
                        #   print(p_stacked)
                        #   dev.off()
                        # }


