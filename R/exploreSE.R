library(shiny)
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)
library(plotly)
library(DT)
library(dplyr)
library(shinyWidgets)
library(colourpicker)
library(ggrepel)

# UI
ui <- fluidPage(
  titlePanel("📊 RNA-seq SummarizedExperiment Explorer"),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Data Input"),
      fileInput("se_file", "Upload SummarizedExperiment (.rds)",
                accept = ".rds"),
      checkboxInput("use_demo", "Use Demo Data", value = TRUE),
      hr(),

      conditionalPanel(
        condition = "output.data_loaded",
        h4("Analysis Options"),
        numericInput("top_genes", "Top variable genes for PCA:",
                     value = 500, min = 100, max = 5000, step = 100),
        hr(),

        h5("DE Analysis"),
        conditionalPanel(
          condition = "output.has_precomputed_de",
          selectInput("de_comparison", "Select Comparison:", choices = NULL)
        ),
      ),
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "main_tabs",
        tabPanel("Overview",
                 h3("Dataset Summary"),
                 verbatimTextOutput("data_summary"),
                 hr(),
                 h4("Sample Metadata"),
                 DTOutput("metadata_table")
        ),

        tabPanel("PCA",
                 h3("Principal Component Analysis"),
                 fluidRow(
                   column(3,
                          selectInput("color_var_1", "Color/Group by:", choices = NULL),
                   )
                 ),
                 plotlyOutput("pca_plot", height = "600px"),
                 hr(),
                 verbatimTextOutput("pca_variance")
        ),

        tabPanel("Gene Expression",
                 h3("Gene Expression Plot"),
                 fluidRow(
                   column(3,
                          checkboxGroupInput("groups_to_show", "Include levels:")
                   ),
                   column(3,
                          pickerInput("gene_id", "Select Gene:",
                                      choices = NULL,
                                      options = list(
                                        `live-search` = TRUE,
                                        `live-search-placeholder` = "Search genes...",
                                        size = 10
                                      )),
                          selectInput("plot_type", "Plot Type:",
                                      choices = c("Boxplot" = "box", "Violin" = "violin"))),
                   column(3,
                          selectInput("color_var_2", "Color/Group by:", choices = NULL),
                   )
                 ),
                 plotlyOutput("expr_plot", height = "500px"),
                 hr(),
                 h4("Expression Values"),
                 DTOutput("expr_table"),
                 downloadButton("download_expr", "Download Expression Table")
        ),

        tabPanel("DE Results",
                 h3("Differential Expression Results"),
                 uiOutput("de_status_message"),
                 fluidRow(

                 ),
                 conditionalPanel(
                   condition = "!output.has_precomputed_de",
                   actionButton("run_de", "Run Basic DE Analysis",
                                class = "btn-primary")
                 ),
                 hr(),
                 DTOutput("de_table"),
                 downloadButton("download_de", "Download DE Table")
        ),

        tabPanel("Volcano Plot",
                 h3("Volcano Plot"),
                 fluidRow(
                   column(3,
                          numericInput("padj_cutoff_volcano", "Adjusted p-value cutoff:",
                                       value = 0.05, min = 0, max = 1, step = 0.01)
                   ),
                   column(3,
                          numericInput("lfc_cutoff_volcano", "Log2 Fold Change cutoff:",
                                       value = 1, min = 0, max = 10, step = 0.5)
                   ),
                   column(3,
                          checkboxInput("label_top", "Label top genes", value = TRUE),
                          numericInput("n_labels", "Number to label:",
                                       value = 10, min = 0, max = 50, step = 5)
                   ),
                   column(3,
                          colourInput("up_col_1", "Color for Upregulated", "#d62728"),
                          colourInput("dn_col_1", "Color for Downregulated", "#1f77b4"),
                          colourInput("high_col_1", "Color for Highlights", "#FFD700"))
                 ),
                 fluidRow(
                   column(12,
                          textAreaInput("highlight_genes",
                                        "Highlight specific genes (one per line or comma-separated):",
                                        value = "",
                                        placeholder = "GENE1, GENE2, GENE3\nor\nGENE1\nGENE2\nGENE3",
                                        rows = 3,
                                        width = "100%"),
                          helpText("Enter gene IDs or gene names to highlight in yellow on the plot.")
                   )
                 ),
                 hr(),
                 plotlyOutput("volcano_plot", height = "700px"),
                 hr(),
                 h4("Summary Statistics"),
                 verbatimTextOutput("volcano_summary")
        ),

        tabPanel("Enrichment Results",
                 h3("Enrichment Results"),
                 fluidRow(
                   column(3,
                          numericInput("padj_cutoff_enrichment", "Adjusted p-value cutoff:",
                                       value = 0.05, min = 0, max = 1, step = 0.01)
                   ),
                   column(3,
                          numericInput("n_terms_enrichment", "Number of top terms to show:",
                                       value = 10, min = 3, max = 20, step = 1)
                   ),
                   column(3,
                          colourInput("up_col_2", "Color for Upregulated", "#d62728"),
                          colourInput("dn_col_2", "Color for Downregulated", "#1f77b4"),
                          colourInput("high_col_2", "Color for Highlights", "#FFD700"))
                 ),
                 uiOutput("fe_status_message"),
                 hr(),
                 uiOutput("enrichment_plots", width = "700px")
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {

  # Increase upload size limit to 500MB
  options(shiny.maxRequestSize = 500*1024^2)

  # Reactive values
  rv <- reactiveValues(
    se = NULL,
    vst_data = NULL,
    pca_result = NULL,
    de_results = NULL,
    up_col = "#d62728",
    dn_col = "#1f77b4",
    highlight_col = "#FFD700",
    color_var = NULL
  )

  # syncing selections ------

  observe({
    # for up
    updateSelectInput(session, "up_col_1", selected = rv$up_col)
    updateSelectInput(session, "up_col_2", selected = rv$up_col)
    # down
    updateSelectInput(session, "dn_col_1", selected = rv$dn_col)
    updateSelectInput(session, "dn_col_2", selected = rv$dn_col)
    # highlight
    updateSelectInput(session, "high_col_1", selected = rv$highlight_col)
    updateSelectInput(session, "high_col_2", selected = rv$highlight_col)
  })
  # up
  observeEvent(input$up_col_1, {
    rv$up_col <- input$up_col_1
  })
  observeEvent(input$up_col_2, {
    rv$up_col <- input$up_col_2
  })
  #down
  observeEvent(input$dn_col_1, {
    rv$dn_col <- input$dn_col_1
  })
  observeEvent(input$dn_col_2, {
    rv$dn_col <- input$dn_col_2
  })
  #highlight
  observeEvent(input$high_col_1, {
    rv$highlight_col <- input$high_col_1
  })
  observeEvent(input$high_col_2, {
    rv$highlight_col <- input$high_col_2
  })

  observe({
    # color var
    updateSelectInput(session, "color_var_1", selected = rv$color_var)
    updateSelectInput(session, "color_var_2", selected = rv$color_var)

    if(!is.null(rv$se)){
      updateCheckboxGroupInput(session, "groups_to_show",
                               choices = levels(colData(rv$se)[[rv$color_var]]),
                               selected = levels(colData(rv$se)[[rv$color_var]]))
    }
  })

  # Load demo data -----
  create_demo_data <- function() {
    set.seed(42)
    n_genes <- 1000
    n_samples <- 12

    # Simulate count data
    counts <- matrix(
      rnbinom(n_genes * n_samples, mu = 100, size = 1/0.5),
      nrow = n_genes
    )
    rownames(counts) <- paste0("GENE", seq_len(n_genes))
    colnames(counts) <- paste0("Sample", seq_len(n_samples))

    # Add some differential expression
    de_genes <- 1:100
    counts[de_genes, 1:6] <- counts[de_genes, 1:6] * 3

    # Sample metadata
    colData <- DataFrame(
      condition = factor(rep(c("Control", "Treatment"), each = 6)),
      batch = factor(rep(c("A", "B"), times = 6)),
      replicate = rep(1:6, times = 2)
    )
    rownames(colData) <- colnames(counts)

    # Gene metadata
    rowData <- DataFrame(
      gene_id = rownames(counts),
      gene_name = paste0("Gene_", seq_len(n_genes))
    )

    SummarizedExperiment(
      assays = list(counts = counts),
      colData = colData,
      rowData = rowData
    )
  }

  # Load data ---------
  observe({
    if (input$use_demo) {
      rv$se <- create_demo_data()
    }
  })

  observeEvent(input$se_file, {
    req(input$se_file)
    tryCatch({
      rv$se <- readRDS(input$se_file$datapath)
      updateCheckboxInput(session, "use_demo", value = FALSE)
      showNotification("Data loaded successfully!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error loading file:", e$message), type = "error")
    })
  })

  # precomputed results ---------
  has_precomputed_de <- reactive({
    req(rv$se)

    if(!is.null(rv$se) && is(rv$se, "DeeDeeExperiment") && !is.null(DeeDeeExperiment::getDEANames(rv$se))){
      TRUE
    }else{
      !is.null(metadata(rv$se)$de_results) &&
        is.list(metadata(rv$se)$de_results) &&
        length(metadata(rv$se)$de_results) > 0
    }
  })

  has_precomputed_fe <- reactive({
    req(rv$se)
    if(!is.null(rv$se) && is(rv$se, "DeeDeeExperiment") && !is.null(DeeDeeExperiment::getFEANames(rv$se))){
      TRUE
    }else{
      !is.null(metadata(rv$se)$fe_results) &&
        is.list(metadata(rv$se)$fe_results) &&
        length(metadata(rv$se)$fe_results) > 0
    }
  })

  # Get precomputed DE comparisons
  de_comparisons <- reactive({
    req(has_precomputed_de())
    if(is(rv$se, "DeeDeeExperiment")){
      DeeDeeExperiment::getDEANames(rv$se)
    }else{
      names(metadata(rv$se)$de_results)
    }
  })

  fes <- reactive({
    req(has_precomputed_fe())
    if(is(rv$se, "DeeDeeExperiment")){
      return(names(DeeDeeExperiment::getFEAList(rv$se, input[[input$de_comparison]])))
    }
    return(names(metadata(rv$se)$fe_results[[input$de_comparison]]))
  })

  # Update UI inputs when data loads
  observe({
    req(rv$se)

    # Update grouping variable choices
    col_vars <- colnames(colData(rv$se))

    # Set default to "condition" if it exists, otherwise first column
    default_var <- if ("condition" %in% col_vars) {
      "condition"
    } else {
      col_vars[1]
    }

    updateSelectInput(session, "color_var_1",
                      choices = col_vars,
                      selected = default_var)
    updateSelectInput(session, "color_var_2",
                      choices = col_vars,
                      selected = default_var)


    # Update gene choices with searchable picker
    gene_choices <- rownames(rv$se)
    if ("gene_name" %in% colnames(rowData(rv$se))) {
      names(gene_choices) <- rowData(rv$se)$gene_name
    }
    updatePickerInput(session, "gene_id",
                      choices = gene_choices,
                      selected = gene_choices[1])

    # Update DE comparison choices if precomputed results exist
    if (has_precomputed_de()) {
      comparisons <- de_comparisons()
      updateSelectInput(session, "de_comparison",
                        choices = comparisons,
                        selected = comparisons[1])
      showNotification(
        paste("Found", length(comparisons), "precomputed DE comparison(s)"),
        type = "message"
      )
    }
  })

  observeEvent(input$color_var_1, {
    rv$color_var <- input$color_var_1
    updateCheckboxGroupInput(session, "groups_to_show",
                             choices = levels(as.factor(colData(rv$se)[[input$color_var_1]])),
                             selected = levels(as.factor(colData(rv$se)[[input$color_var_1]])))
  })
  observeEvent(input$color_var_2, {
    rv$color_var <- input$color_var_2
    updateCheckboxGroupInput(session, "groups_to_show",
                             choices = levels(as.factor(colData(rv$se)[[input$color_var_2]])),
                             selected = levels(as.factor(colData(rv$se)[[input$color_var_2]])))
  })


  # VST transformation -------
  vst_data <- reactive({
    req(rv$se)

    if (is.null(rv$vst_data)) {
      withProgress(message = "Transforming data...", {
        dds <- DESeqDataSet(rv$se, design = ~ 1)
        rv$vst_data <- assay(vst(dds, blind = TRUE))
      })
    }
    rv$vst_data
  })

  # PCA calculation
  pca_data <- reactive({
    req(vst_data(), rv$color_var, input$top_genes)

    vst_mat <- vst_data()

    # Select top variable genes
    rv_genes <- rowVars(vst_mat)
    select_genes <- order(rv_genes, decreasing = TRUE)[1:min(input$top_genes, nrow(vst_mat))]

    # Run PCA
    pca <- prcomp(t(vst_mat[select_genes, ]), scale. = FALSE)
    rv$pca_result <- pca

    # Create data frame for plotting
    pca_df <- data.frame(
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      sample = colnames(rv$se),
      group = colData(rv$se)[[rv$color_var]]
    )

    # Calculate variance explained
    var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

    list(data = pca_df, var = var_explained)
  })

  # Run DE analysis ----------
  observeEvent(input$run_de, {
    req(rv$se, rv$color_var)

    withProgress(message = "Running DESeq2...", {
      tryCatch({
        # Create DESeq2 object
        col_var <- rv$color_var
        design_formula <- as.formula(paste("~", col_var))
        dds <- DESeqDataSet(rv$se, design = design_formula)

        # Filter low counts
        keep <- rowSums(counts(dds)) >= 10
        dds <- dds[keep, ]

        # Run DESeq2
        dds <- DESeq(dds)

        # Get results
        res <- results(dds)
        rv$de_results <- as.data.frame(res) %>%
          tibble::rownames_to_column("gene_id") %>%
          arrange(padj)

        showNotification("DE analysis complete!", type = "message")
      }, error = function(e) {
        showNotification(paste("Error in DE analysis:", e$message), type = "error")
      })
    })
  })

  # current DE results ---------
  current_de_results <- reactive({
    if (has_precomputed_de() && !is.null(input$de_comparison)) {
      # Load precomputed results
      if(is(rv$se, "DeeDeeExperiment")){
        de_res <- DeeDeeExperiment::getDEA(rv$se, dea_name = input$de_comparison)
        colnames(de_res) <- c("log2FoldChange", "pvalue", "padj")
      }else{
        de_res <- metadata(rv$se)$de_results[[input$de_comparison]]
      }


      # Ensure it's a data frame
      if (is(de_res, "DESeqResults")) {
        de_res <- as.data.frame(de_res) %>%
          tibble::rownames_to_column("gene_id")
      } else if (!"gene_id" %in% colnames(de_res) && !is.null(rownames(de_res))) {
        de_res <- de_res %>% as.data.frame() %>%
          tibble::rownames_to_column("gene_id")
      }
      return(de_res)
    } else if (!is.null(rv$de_results)) {
      # Use computed results
      return(rv$de_results)
    } else {
      return(NULL)
    }
  })

  # current FE results ------
  current_fe_results <- reactive({
    if (has_precomputed_fe()) {
      # Load precomputed results
      if(is(rv$se, "DeeDeeExperiment")){
        fe_res <- DeeDeeExperiment::getFEAList(rv$se, dea_name = input$de_comparison)
      }else{
        fe_res <- metadata(rv$se)$fe_results[[input$de_comparison]]
      }
      return(fe_res)
    } else if (!is.null(rv$fe_results)) {
      # Use computed results
      return(rv$fe_results)
    } else {
      return(NULL)
    }
  })

  # Outputs ---------------
  output$data_loaded <- reactive({
    !is.null(rv$se)
  })
  outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)

  output$has_precomputed_de <- reactive({
    has_precomputed_de()
  })
  output$has_precomputed_fe <- reactive({
    has_precomputed_fe()
  })
  outputOptions(output, "has_precomputed_de", suspendWhenHidden = FALSE)
  outputOptions(output, "has_precomputed_fe", suspendWhenHidden = FALSE)

  output$de_status_message <- renderUI({
    if (has_precomputed_de()) {
      comparisons <- de_comparisons()
      tagList(
        p(
          icon("check-circle", class = "text-success"),
          strong(paste("Found", length(comparisons), "precomputed DE comparison(s):")),
          br(),
          paste(comparisons, collapse = ", ")
        )
      )
    } else {
      p(
        icon("info-circle"),
        "No precomputed DE results found. Run basic DE analysis or upload data with results in metadata(se)$de_results."
      )
    }
  })
  output$fe_status_message <- renderUI({
    if (has_precomputed_fe()) {
      comparisons <- names(current_fe_results())
      tagList(
        p(
          icon("check-circle", class = "text-success"),
          strong(paste("Found", length(comparisons), "precomputed functional enrichment(s):")),
          br(),
          paste(comparisons, collapse = ", ")
        )
      )
    } else {
      p(
        icon("info-circle"),
        "No precomputed FE results found. Run basic FE analysis or upload data with results in metadata(se)$fe_results."
      )
    }
  })

  output$data_summary <- renderPrint({
    req(rv$se)
    cat("SummarizedExperiment Object\n")
    cat("===========================\n\n")
    cat("Dimensions:\n")
    cat("  Genes:", nrow(rv$se), "\n")
    cat("  Samples:", ncol(rv$se), "\n\n")
    cat("Assays:", paste(names(assays(rv$se)), collapse = ", "), "\n\n")
    cat("Sample Metadata Columns:\n")
    cat(" ", paste(colnames(colData(rv$se)), collapse = ", "), "\n\n")
    if (ncol(rowData(rv$se)) > 0) {
      cat("Gene Metadata Columns:\n")
      cat(" ", paste(colnames(rowData(rv$se)), collapse = ", "), "\n")
    }
  })

  output$metadata_table <- renderDT({
    req(rv$se)
    datatable(
      as.data.frame(colData(rv$se)),
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = TRUE
    )
  })

  output$pca_plot <- renderPlotly({
    req(pca_data())

    pca_info <- pca_data()
    df <- pca_info$data
    var <- pca_info$var

    p <- ggplot(df, aes(x = PC1, y = PC2, color = group, text = sample)) +
      geom_point(size = 4, alpha = 0.8) +
      labs(
        x = paste0("PC1 (", var[1], "%)"),
        y = paste0("PC2 (", var[2], "%)"),
        color = rv$color_var
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "right")

    ggplotly(p, tooltip = c("text", "group")) %>%
      layout(hovermode = "closest")
  })

  output$pca_variance <- renderPrint({
    req(rv$pca_result)
    var_explained <- round(100 * rv$pca_result$sdev^2 / sum(rv$pca_result$sdev^2), 2)
    cat("Variance Explained by PCs:\n")
    for (i in 1:min(10, length(var_explained))) {
      cat(sprintf("  PC%d: %.2f%%\n", i, var_explained[i]))
    }
  })
  ## expression plot-------
  output$expr_plot <- renderPlotly({
    req(rv$se, input$gene_id, rv$color_var)

    gene <- input$gene_id
    counts_data <- assay(rv$se, "counts")[gene, ]

    plot_df <- data.frame(
      s_a_m_p_l_e = colnames(rv$se),
      expression = counts_data,
      group = forcats::as_factor(colData(rv$se)[[rv$color_var]])
    ) %>%
      dplyr::filter(group %in% input$groups_to_show)

    gene_label <- if ("gene_name" %in% colnames(rowData(rv$se))) {
      rowData(rv$se)[gene, "gene_name"]
    } else {
      gene
    }

    p <- ggplot(plot_df, aes(x = group, y = expression, fill = group, text = s_a_m_p_l_e)) +
      labs(
        title = paste("Expression:", gene_label),
        x = rv$color_var,
        y = "Normalized Counts"
      ) +
      scale_x_discrete(labels = \(x) stringr::str_wrap(stringr::str_replace_all(x, "_", " "), 10))+
      theme_minimal(base_size = 14) +
      theme(legend.position = "none")

    if (input$plot_type == "box") {
      p <- p + geom_boxplot(alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 2)
    } else {
      p <- p + geom_violin(alpha = 0.7) +
        geom_jitter(width = 0.1, alpha = 0.5, size = 2)
    }

    ggplotly(p)
  })

  output$expr_table <- renderDT({
    req(rv$se, input$gene_id, input$groups_to_show)

    gene <- input$gene_id
    expr_data <- data.frame(
      Sample = colnames(rv$se),
      colData(rv$se),
      Gene = input$gene_id,
      Count = assay(rv$se, "counts")[gene, ],
      g_r_o_u_p = colData(rv$se)[[rv$color_var]]
    )%>%
      dplyr::filter(g_r_o_u_p %in% input$groups_to_show) %>%
      dplyr::select(-g_r_o_u_p)

    datatable(
      expr_data,
      options = list(pageLength = 12, scrollX = TRUE),
      rownames = FALSE
    ) %>%
      formatRound("Count", digits = 0)
  })

  output$de_table <- renderDT({
    de_data <- current_de_results()
    de_data <- dplyr::select(de_data, any_of(c("gene_id", "baseMean", "log2FoldChange", "pvalue", "padj"))) %>%
      mutate(across(any_of(c("baseMean", "log2FoldChange")), \(x) round(x, 2)),
             across(any_of(c("pvalue", "padj")), \(x) round(x, 4)))
    req(de_data)

    datatable(
      de_data,
      options = list(pageLength = 25, scrollX = TRUE),
      rownames = FALSE,
      filter = "top"
    )
  })

  output$download_expr <- downloadHandler(
    filename = function() {
      comparison_name <- stringr::str_replace(input$gene_id, "[^a-zA-Z0-9_-]", "_")

      paste0(comparison_name, "_expression_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(rv$se, input$gene_id, input$groups_to_show)

      gene <- input$gene_id
      expr_data <- data.frame(
        Sample = colnames(rv$se),
        colData(rv$se),
        Gene = input$gene_id,
        Count = assay(rv$se, "counts")[gene, ],
        g_r_o_u_p = colData(rv$se)[[rv$color_var]]
      )%>%
        dplyr::filter(g_r_o_u_p %in% input$groups_to_show) %>%
        dplyr::select(-g_r_o_u_p)
      readr::write_excel_csv2(expr_data, file)
    }
  )

  output$download_de <- downloadHandler(
    filename = function() {
      comparison_name <- if (has_precomputed_de() && !is.null(input$de_comparison)) {
        gsub("[^a-zA-Z0-9_-]", "_", input$de_comparison)
      } else {
        "DE_results"
      }
      paste0(comparison_name, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      de_data <- current_de_results()
      req(de_data)
      readr::write_excel_csv2(de_data, file)
    }
  )

  ## Volcano plot ----
  output$volcano_plot <- renderPlotly({
    de_data <- current_de_results()
    req(de_data, input$padj_cutoff_volcano, input$lfc_cutoff_volcano, rv$up_col, rv$dn_col, rv$highlight_col)
    colors_acute <- c("Up" = rv$up_col, "Down" = rv$dn_col, "NS" = "grey70", "Highlighted" = rv$highlight_col)

    # Prepare data
    volcano_df <- de_data %>%
      dplyr::filter(!is.na(padj) & !is.na(log2FoldChange) & !is.infinite(log2FoldChange)) %>%
      mutate(
        neg_log10_padj = -log10(padj),
        sig = case_when(
          padj < input$padj_cutoff_volcano & log2FoldChange > input$lfc_cutoff_volcano ~ "Up",
          padj < input$padj_cutoff_volcano & log2FoldChange < -input$lfc_cutoff_volcano ~ "Down",
          TRUE ~ "NS"
        )
      )

    # Get gene names for hover
    if ("gene_name" %in% colnames(rowData(rv$se))) {
      gene_lookup <- setNames(rowData(rv$se)$gene_name, rownames(rv$se))
      volcano_df$gene_name <- gene_lookup[volcano_df$gene_id]
    } else {
      volcano_df$gene_name <- volcano_df$gene_id
    }

    # Parse highlighted genes
    highlight_list <- c()
    if (!is.null(input$highlight_genes) && nchar(trimws(input$highlight_genes)) > 0) {
      # Split by newlines and commas, trim whitespace
      highlight_list <- input$highlight_genes %>%
        strsplit("[\n,]") %>%
        unlist() %>%
        trimws() %>%
        .[nchar(.) > 0]
    }

    # Mark highlighted genes
    volcano_df <- volcano_df %>%
      mutate(
        highlighted = gene_id %in% highlight_list | gene_name %in% highlight_list,
        display_category = case_when(
          highlighted ~ "Highlighted",
          sig == "Up" ~ "Up",
          sig == "Down" ~ "Down",
          TRUE ~ "NS"
        )
      )

    # Identify top genes to label (excluding highlighted genes since they'll be labeled anyway)
    genes_to_label <- data.frame()
    if (input$label_top && input$n_labels > 0) {
      top_genes <- volcano_df %>%
        dplyr::filter(sig != "NS" & !highlighted) %>%
        arrange(padj) %>%
        head(input$n_labels)
      genes_to_label <- rbind(genes_to_label, top_genes)
    }

    # Always label highlighted genes
    if (length(highlight_list) > 0) {
      highlighted_genes <- volcano_df %>%
        dplyr::filter(highlighted)
      genes_to_label <- rbind(genes_to_label, highlighted_genes)
    }

    # Color scheme - highlighted genes in bright yellow/gold

    # Reorder so highlighted genes are plotted on top
    volcano_df <- volcano_df %>%
      arrange(highlighted)

    # Create plot
    p <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj),
                                color = display_category, text = paste("log2 FC:", round(log2FoldChange, 2), "\nGene:", gene_name, "\nadjusted p-value:", format(padj, digits=4)))) +
      geom_point(aes(size = highlighted, alpha = ifelse(highlighted, 1, 0.6))) +
      scale_size_manual(values = c("TRUE" = 3, "FALSE" = 1.5), guide = "none") +
      scale_alpha_identity() +
      scale_color_manual(values = colors_acute,
                         name = "Category",
                         breaks = c("Up", "Down", "Highlighted", "NS"),
                         labels = c("Up-regulated", "Down-regulated", "Highlighted", "Not significant")) +
      geom_vline(xintercept = c(-input$lfc_cutoff_volcano, input$lfc_cutoff_volcano),
                 linetype = "dashed", color = "grey30", linewidth = 0.5) +
      geom_hline(yintercept = -log10(input$padj_cutoff_volcano),
                 linetype = "dashed", color = "grey30", linewidth = 0.5) +
      labs(
        x = "Log2 Fold Change",
        y = "-Log10 Adjusted P-value",
        title = if (has_precomputed_de() && !is.null(input$de_comparison)) {
          paste("Volcano Plot:", input$de_comparison)
        } else {
          "Volcano Plot"
        }
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "right",
        panel.grid.minor = element_blank()
      )

    # Add labels
    if (nrow(top_genes) > 0) {
      p <- p + geom_text(
        data = top_genes,
        aes(label = gene_name),
        size = 5,
        show.legend = FALSE
      )
    }

    ggplotly(p, tooltip = c("text")) %>%
      layout(hovermode = "closest") %>%
      style(textposition = "right")
  })
  output$volcano_summary <- renderPrint({
    de_data <- current_de_results()
    req(de_data, input$padj_cutoff_volcano, input$lfc_cutoff_volcano, rv$up_col, rv$dn_col, rv$highlight_col)

    colors_acute <- c("Up" = rv$up_col, "Down" = rv$dn_col, "NS" = "grey70", "Highlighted" = rv$highlight_col)
    summary_df <- de_data %>%
      dplyr::filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
      summarise(
        total_genes = n(),
        sig_genes = sum(padj < input$padj_cutoff_volcano),
        up_regulated = sum(padj < input$padj_cutoff_volcano & log2FoldChange > input$lfc_cutoff_volcano),
        down_regulated = sum(padj < input$padj_cutoff_volcano & log2FoldChange < -input$lfc_cutoff_volcano),
        not_significant = sum(padj >= input$padj_cutoff_volcano |
                                (abs(log2FoldChange) < input$lfc_cutoff_volcano & padj < input$padj_cutoff_volcano))
      )

    cat("Differential Expression Summary\n")
    cat("================================\n\n")
    cat(sprintf("Total genes tested: %d\n", summary_df$total_genes))
    cat(sprintf("Significant (padj < %.3f): %d (%.1f%%)\n",
                input$padj_cutoff_volcano,
                summary_df$sig_genes,
                100 * summary_df$sig_genes / summary_df$total_genes))
    cat(sprintf("\nUp-regulated (LFC > %.2f): %d\n",
                input$lfc_cutoff_volcano, summary_df$up_regulated))
    cat(sprintf("Down-regulated (LFC < -%.2f): %d\n",
                input$lfc_cutoff_volcano, summary_df$down_regulated))
    cat(sprintf("Not significant: %d\n", summary_df$not_significant))
  })
  ## enrichment plots -------
  output$enrichment_plots <- renderUI({
    current_fes <- current_fe_results()
    colors_acute <- c("Up" = rv$up_col, "Down" = rv$dn_col, "NS" = "grey70", "Highlighted" = rv$highlight_col)
    req(current_fes, input$padj_cutoff_enrichment, input$n_terms_enrichment, colors_acute)


    plots <- purrr::imap(current_fes, function(enrich, name) {

      id <- paste0("e_plot_", name)
      DIR <- ifelse(stringr::str_detect(name, "up") | stringr::str_detect(name, "UP"),
                    "Up", "Down")
      plotlyOutput(outputId = id, width = "700px")

      output[[id]] <- renderPlotly({
        fe_data <- enrich
        if(stringr::str_detect(name, "[Gg][Oo]")){
          go_df <- fe_data %>%
            dplyr::filter(p.adjust < input$padj_cutoff_enrichment) %>%
            dplyr::slice_max(FoldEnrichment, n = input$n_terms_enrichment) %>%
            dplyr::mutate(Description = forcats::fct_reorder(Description, FoldEnrichment),
                          dir = DIR)
          p2 <- ggplot(go_df, aes(FoldEnrichment,
                                  Description, text = stringr::str_wrap(stringr::str_replace_all(geneID, "\\/", ", "),
                                                                        width = 60)))+
            geom_col(aes(fill = dir))+
            scale_y_discrete(labels = \(x) stringr::str_wrap(x, width = 60))+
            scale_fill_manual(values = colors_acute)+
            theme_light(base_size = 14)+
            labs(x = "Fold Enrichment over Background",
                 title = if (DIR == "Up") {
                   paste("Top", input$n_terms_enrichment, "upregulated GO Terms")
                 }else {
                   paste("Top", input$n_terms_enrichment, "downregulated GO Terms")
                 })+
            theme(axis.title.y = element_blank(), legend.position = "none")

          ggplotly(p2, tooltip = c("text")) %>%
            layout(hovermode = "closest") %>%
            style(textposition = "right")
        }else if(stringr::str_detect(name, "gsea")){
          gsea_df <- fe_data %>%
            dplyr::filter(p.adjust < input$padj_cutoff_enrichment) %>%
            group_by(sign(NES)) %>%
            dplyr::slice_max(abs(NES), n = round(input$n_terms_enrichment / 2)) %>%
            ungroup() %>%
            dplyr::mutate(Description = stringr::str_remove(Description, "HALLMARK_") %>%
                            stringr::str_remove("REACTOME") %>%
                            stringr::str_replace_all("_", " ") %>%
                            forcats::fct_reorder(NES),
                          dir = ifelse(NES > 0, "Up", "Down"))

          p2 <- ggplot(gsea_df, aes(NES, Description))+
            geom_col(aes(fill = dir))+
            scale_y_discrete(labels = \(x) stringr::str_wrap(x, width = 60))+
            scale_fill_manual(values = colors_acute)+
            theme_light(base_size = 14)+
            labs(x = "NES (Normalized Enrichment Score)",
                 title = if (stringr::str_detect(name, "HALLMARK")) {
                   paste("Top", input$n_terms_enrichment, "enriched hallmark gene sets")
                 }else {
                   paste("Top", input$n_terms_enrichment, "enriched Reactome gene sets")
                 })+
            theme(axis.title.y = element_blank(), legend.position = "none")

          ggplotly(p2, tooltip = c("text")) %>%
            layout(hovermode = "closest") %>%
            style(textposition = "right")
        }
      })

    }

    )

  })
}

# Run app
exploreSE <- function(){
  shinyApp(ui, server)
}
