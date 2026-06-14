# ============================================
> # STEP 2 - Load and clean data
> # ============================================
> 
> data2 <- read_tsv("/Users/ivancarrier/Downloads/plotCSMD1.txt") %>%
+   rename(
+     ATM_CNA    = `ATM: Putative copy-number alterations from DNAcopy.`,
+     MYC_expr = `MYC: mRNA expression (Illumina HT-12 v3 microarray)`
+   ) %>%
+   filter(!is.na(ATM_CNA), !is.na(MYC_expr))
Rows: 1866 Columns: 4                                                                                                                                                                                              
── Column specification ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
Delimiter: "\t"
chr (3): Sample Id, ATM: Putative copy-number alterations from DNAcopy., Copy Number Alterations
dbl (1): CSMD1: mRNA expression (Illumina HT-12 v3 microarray)

ℹ Use `spec()` to retrieve the full column specification for this data.
ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
Error in `rename()`:
! Can't rename columns that don't exist.
✖ Column `MYC: mRNA expression (Illumina HT-12 v3 microarray)` doesn't exist.
Run `rlang::last_trace()` to see where the error occurred.
> 
> 
> # ============================================
> # STEP 3 - Set factor levels and colors
> # ============================================
> 
> data2$ATM_CNA <- factor(data2$ATM_CNA,
+                          levels = c("Deep Deletion",
+                                     "Shallow Deletion",
+                                     "Diploid",
+                                     "Gain",
+                                     "Amplification"))
> 
> cna_colors <- c(
+   "Deep Deletion"    = "#0000FF",
+   "Shallow Deletion" = "#00BCD4",
+   "Diploid"          = "#808080",
+   "Gain"             = "#FF69B4",
+   "Amplification"    = "#FF0000"
+ )
> 
> 
> # ============================================
> # STEP 4 - Define comparisons
> # ============================================
> 
> comparisons <- list(
+   c("Diploid", "Deep Deletion"),
+   c("Diploid", "Shallow Deletion"),
+   c("Diploid", "Gain"),
+   c("Diploid", "Amplification")
+ )
> 
> 
> # ============================================
> # STEP 5 - Build plot
> # ============================================
> 
> p2 <- ggplot(data2, aes(x = ATM_CNA, y = MYC_expr, color = ATM_CNA)) +
+ 
+   geom_boxplot(
+     outlier.shape = NA,
+     fill = "lightgray",
+     color = "black",
+     width = 0.5
+   ) +
+ 
+   geom_jitter(
+     shape = 1,
+     size = 1.5,
+     width = 0.2,
+     alpha = 0.7
+   ) +
+ 
+   scale_color_manual(values = cna_colors) +
+ 
+   stat_compare_means(
+     comparisons = comparisons,
+     method = "wilcox.test",
+     label = "p.format"
+   ) +
+ 
+   labs(
+     x = "ATM: Putative copy-number alterations from DNAcopy",
+     y = "MYC: mRNA expression (Illumina HT-12 v3 microarray)",
+     title = "ATM Copy Number vs MYC mRNA Expression",
+     subtitle = paste0("METABRIC — ", nrow(data2), " samples")
+   ) +
+ 
+   theme_classic() +
+   theme(
+     legend.position = "none",
+     axis.text.x = element_text(angle = 45, hjust = 1),
+     axis.title = element_text(size = 10),
+     plot.title = element_text(size = 12, face = "bold")
+   )
> 
> 
> # ============================================
> # STEP 6 - View and save
> # ============================================
> 
> print(p2)
> 
> ggsave(
+   "MYC_vs_ATM_CNA_pvalues1.png",
+   plot = p2,
+   width = 7,
+   height = 6,
+   dpi = 300
+ )
> 
> ggsave(
+   "MYC_vs_ATM_CNA_pvalues1.pdf",
+   plot = p2,
+   width = 7,
+   height = 6
+ )
> 
> 
> # ============================================
> # STEP 7 - Get exact p-values
> # ============================================
> 
> MYC_diploid       <- data2 %>% filter(ATM_CNA == "Diploid")          %>% pull(MYC_expr)
> MYC_shallow       <- data2 %>% filter(ATM_CNA == "Shallow Deletion")  %>% pull(MYC_expr)
> MYC_deep          <- data2 %>% filter(ATM_CNA == "Deep Deletion")     %>% pull(MYC_expr)
> MYC_gain          <- data2 %>% filter(ATM_CNA == "Gain")              %>% pull(MYC_expr)
> MYC_amplification <- data2 %>% filter(ATM_CNA == "Amplification")     %>% pull(MYC_expr)
> 
