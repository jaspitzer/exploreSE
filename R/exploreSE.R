# UI
ui <- shiny::fluidPage(
  shiny::titlePanel("RNA-seq SummarizedExperiment Explorer"),

  shiny::sidebarLayout(
    shiny::sidebarPanel(
      width = 3,
      shiny::h4("Data Input"),
      shiny::fileInput(
        "se_file",
        "Upload SummarizedExperiment (.rds)",
        accept = ".rds"
      ),
      shiny::checkboxInput("use_demo", "Use Demo Data", value = TRUE),
      shiny::hr(),

      shiny::conditionalPanel(
        condition = "output.data_loaded",
        shiny::h4("Analysis Options"),

        shiny::h5("DE Analysis"),
        shiny::conditionalPanel(
          condition = "output.has_precomputed_de",
          shiny::selectInput(
            "de_comparison",
            "Select Comparison:",
            choices = NULL
          )
        ),
      ),
    ),

    shiny::mainPanel(
      width = 9,
      shiny::tabsetPanel(
        id = "main_tabs",
        shiny::tabPanel(
          "Overview",
          shiny::h3("Dataset Summary"),
          shiny::verbatimTextOutput("data_summary"),
          shiny::hr(),
          shiny::h4("Sample Metadata"),
          DT::DTOutput("metadata_table")
        ),

        shiny::tabPanel(
          "PCA",
          shiny::h3("Principal Component Analysis"),
          shiny::fluidRow(
            shiny::column(
              3,
              shiny::selectInput(
                "color_var_1",
                "Color/Group by:",
                choices = NULL
              ),
              shiny::numericInput(
                "top_genes",
                "Top variable genes for PCA:",
                value = 500,
                min = 100,
                max = 5000,
                step = 100
              ),
              shiny::hr(),
            )
          ),
          plotly::plotlyOutput("pca_plot", height = "600px"),
          shiny::hr(),
          shiny::verbatimTextOutput("pca_variance")
        ),

        shiny::tabPanel(
          "Gene Expression",
          shiny::h3("Gene Expression Plot"),
          shiny::fluidRow(
            shiny::column(
              3,
              shiny::checkboxGroupInput("groups_to_show", "Include levels:")
            ),
            shiny::column(
              3,
              shinyWidgets::pickerInput(
                "gene_id",
                "Select Gene:",
                choices = NULL,
                options = list(
                  `live-search` = TRUE,
                  `live-search-placeholder` = "Search genes...",
                  size = 10
                )
              ),
              shiny::selectInput(
                "plot_type",
                "Plot Type:",
                choices = c("Boxplot" = "box", "Violin" = "violin")
              )
            ),
            shiny::column(
              3,
              shiny::selectInput(
                "color_var_2",
                "Color/Group by:",
                choices = NULL
              ),
            )
          ),
          plotly::plotlyOutput("expr_plot", height = "500px"),
          shiny::hr(),
          shiny::h4("Expression Values"),
          DT::DTOutput("expr_table"),
          shiny::downloadButton("download_expr", "Download Expression Table")
        ),

        shiny::tabPanel(
          "DE Results",
          shiny::h3("Differential Expression Results"),
          shiny::uiOutput("de_status_message"),
          shiny::fluidRow(),
          shiny::conditionalPanel(
            condition = "!output.has_precomputed_de",
            shiny::actionButton(
              "run_de",
              "Run Basic DE Analysis",
              class = "btn-primary"
            )
          ),
          shiny::hr(),
          shiny::plotOutput("de_plot", height = "700px"),
          shiny::hr(),
          DT::DTOutput("de_table"),
          shiny::downloadButton("download_de", "Download DE Table")
        ),

        shiny::tabPanel(
          "Volcano Plot",
          shiny::h3("Volcano Plot"),
          shiny::fluidRow(
            shiny::column(
              3,
              shiny::numericInput(
                "padj_cutoff_volcano",
                "Adjusted p-value cutoff:",
                value = 0.05,
                min = 0,
                max = 1,
                step = 0.01
              )
            ),
            shiny::column(
              3,
              shiny::numericInput(
                "lfc_cutoff_volcano",
                "Log2 Fold Change cutoff:",
                value = 1,
                min = 0,
                max = 10,
                step = 0.5
              )
            ),
            shiny::column(
              3,
              shiny::checkboxInput(
                "label_top",
                "Label top genes",
                value = TRUE
              ),
              shiny::numericInput(
                "n_labels",
                "Number to label:",
                value = 10,
                min = 0,
                max = 50,
                step = 5
              )
            ),
            shiny::column(
              3,
              colourpicker::colourInput(
                "up_col_1",
                "Color for Upregulated",
                "#d62728"
              ),
              colourpicker::colourInput(
                "dn_col_1",
                "Color for Downregulated",
                "#1f77b4"
              ),
              colourpicker::colourInput(
                "high_col_1",
                "Color for Highlights",
                "#FFD700"
              )
            )
          ),
          shiny::fluidRow(
            shiny::column(
              12,
              shiny::textAreaInput(
                "highlight_genes",
                "Highlight specific genes (one per line or comma-separated):",
                value = "",
                placeholder = "GENE1, GENE2, GENE3\nor\nGENE1\nGENE2\nGENE3",
                rows = 3,
                width = "100%"
              ),
              shiny::helpText(
                "Enter gene IDs or gene names to highlight in yellow on the plot."
              )
            )
          ),
          shiny::hr(),
          plotly::plotlyOutput("volcano_plot", height = "700px"),
          shiny::hr(),
          shiny::h4("Summary Statistics"),
          shiny::verbatimTextOutput("volcano_summary")
        ),

        shiny::tabPanel(
          "Enrichment Results",
          shiny::h3("Enrichment Results"),
          shiny::fluidRow(
            shiny::column(
              3,
              shiny::numericInput(
                "padj_cutoff_enrichment",
                "Adjusted p-value cutoff:",
                value = 0.05,
                min = 0,
                max = 1,
                step = 0.01
              )
            ),
            shiny::column(
              3,
              shiny::numericInput(
                "n_terms_enrichment",
                "Number of top terms to show:",
                value = 10,
                min = 3,
                max = 20,
                step = 1
              )
            ),
            shiny::column(
              3,
              colourpicker::colourInput(
                "up_col_2",
                "Color for Upregulated",
                "#d62728"
              ),
              colourpicker::colourInput(
                "dn_col_2",
                "Color for Downregulated",
                "#1f77b4"
              ),
              colourpicker::colourInput(
                "high_col_2",
                "Color for Highlights",
                "#FFD700"
              )
            )
          ),
          shiny::uiOutput("fe_status_message"),
          shiny::hr(),
          shiny::uiOutput("enrichment_plots", width = "700px")
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  # Increase upload size limit to 500MB
  options(shiny.maxRequestSize = 500 * 1024^2)

  # Reactive values
  rv <- shiny::reactiveValues(
    se = NULL,
    vst_data = NULL,
    pca_result = NULL,
    de_results = NULL,
    up_col = "#d62728",
    dn_col = "#1f77b4",
    highlight_col = "#FFD700",
    color_var = NULL
  )

  has_precomputed_de <- shiny::reactive({
    shiny::req(rv$se)
    .check_precomputed_de(rv$se)
  })

  has_precomputed_fe <- shiny::reactive({
    shiny::req(rv$se)
    .check_precomputed_fe(rv$se)
  })

  .create_dir_color_observers(input, session, rv)
  .create_interest_color_observers(input, session, rv)
  .observe_demo_data(input, session, rv)
  .observe_inital_obj(input, session, rv)
  .observe_load_file(input, session, rv)

  # precomputed results ---------

  # Get precomputed DE comparisons
  de_comparisons <- shiny::reactive({
    shiny::req(has_precomputed_de())
    .de_results_names(rv$se)
  })

  fes <- shiny::reactive({
    shiny::req(has_precomputed_fe())
    .fe_results_names(rv$se, input$de_comparison)
  })

  # Update UI inputs when data loads

  # observes the SE object and triggers a reload of the the UI, dropdown choices, etc

  .observe_se_load(input, session, rv)

  # VST transformation -------
  vst_data <- shiny::reactive({
    shiny::req(rv$se)

    if (is.null(rv$vst_data)) {
      shiny::withProgress(message = "Transforming data...", {
        dds <- DESeq2::DESeqDataSet(rv$se, design = ~1)
        rv$vst_data <- SummarizedExperiment::assay(DESeq2::vst(
          dds,
          blind = TRUE
        ))
      })
    }
    rv$vst_data
  })

  # PCA calculation
  pca_data <- shiny::reactive({
    shiny::req(vst_data(), rv$color_var, input$top_genes)

    vst_mat <- vst_data()

    # Select top variable genes
    rv_genes <- matrixStats::rowVars(vst_mat)
    select_genes <- order(rv_genes, decreasing = TRUE)[
      seq_len(min(input$top_genes, nrow(vst_mat)))
    ]

    # Run PCA
    pca <- stats::prcomp(t(vst_mat[select_genes, ]), scale. = FALSE)
    rv$pca_result <- pca

    # Create data frame for plotting
    pca_df <- data.frame(
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      sample = colnames(rv$se),
      group = SummarizedExperiment::colData(rv$se)[[rv$color_var]]
    )

    # Calculate variance explained
    var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

    list(data = pca_df, var = var_explained)
  })

  # Run DE analysis ----------
  shiny::observeEvent(input$run_de, {
    shiny::req(rv$se, rv$color_var)

    shiny::withProgress(message = "Running DESeq2...", {
      tryCatch(
        {
          # Create DESeq2 object
          col_var <- rv$color_var
          design_formula <- stats::as.formula(paste("~", col_var))
          dds <- DESeq2::DESeqDataSet(rv$se, design = design_formula)

          # Filter low counts
          keep <- rowSums(BiocGenerics::counts(dds)) >= 10
          dds <- dds[keep, ]

          # Run DESeq2
          dds <- DESeq2::DESeq(dds)

          # Get results
          res <- DESeq2::results(dds)
          rv$de_results <- as.data.frame(res) %>%
            tibble::rownames_to_column("gene_id") %>%
            dplyr::arrange(padj)

          shiny::showNotification("DE analysis complete!", type = "message")
        },
        error = function(e) {
          shiny::showNotification(
            paste("Error in DE analysis:", e$message),
            type = "error"
          )
        }
      )
    })
  })

  # current DE results ---------
  current_de_results <- shiny::reactive({
    if (has_precomputed_de() && !is.null(input$de_comparison)) {
      # Load precomputed results

      de_res <- .de_result(rv$se, input$de_comparison) %>%
        tibble::rownames_to_column("gene_id")
      return(de_res)
    } else if (!is.null(rv$de_results)) {
      # Use computed results
      return(rv$de_results)
    } else {
      return(NULL)
    }
  })

  # current FE results ------
  current_fe_results <- shiny::reactive({
    if (has_precomputed_fe()) {
      # Load precomputed results
      if (methods::is(rv$se, "DeeDeeExperiment")) {
        fe_res <- .fe_result(rv$se, input$de_comparison)
      } else {
        fe_res <- S4Vectors::metadata(rv$se)$fe_results[[input$de_comparison]]
      }
      return(fe_res)
    } else {
      return(NULL)
    }
  })

  # Outputs ---------------
  output$data_loaded <- shiny::reactive({
    !is.null(rv$se)
  })
  shiny::outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)

  output$has_precomputed_de <- shiny::reactive({
    has_precomputed_de()
  })
  output$has_precomputed_fe <- shiny::reactive({
    has_precomputed_fe()
  })
  shiny::outputOptions(output, "has_precomputed_de", suspendWhenHidden = FALSE)
  shiny::outputOptions(output, "has_precomputed_fe", suspendWhenHidden = FALSE)

  output$de_status_message <- shiny::renderUI({
    if (has_precomputed_de()) {
      comparisons <- de_comparisons()
      htmltools::tagList(
        htmltools::p(
          shiny::icon("check-circle", class = "text-success"),
          htmltools::strong(paste(
            "Found",
            length(comparisons),
            "precomputed DE comparison(s):"
          )),
          htmltools::br(),
          paste(comparisons, collapse = ", ")
        )
      )
    } else {
      htmltools::p(
        shiny::icon("info-circle"),
        "No precomputed DE results found. Run basic DE analysis or upload data with results in metadata(se)$de_results."
      )
    }
  })
  output$fe_status_message <- shiny::renderUI({
    if (has_precomputed_fe()) {
      comparisons <- names(current_fe_results())
      htmltools::tagList(
        htmltools::p(
          shiny::icon("check-circle", class = "text-success"),
          htmltools::strong(paste(
            "Found",
            length(comparisons),
            "precomputed functional enrichment(s):"
          )),
          htmltools::br(),
          paste(comparisons, collapse = ", ")
        )
      )
    } else {
      htmltools::p(
        shiny::icon("info-circle"),
        "No precomputed FE results found. Run basic FE analysis or upload data with results in metadata(se)$fe_results."
      )
    }
  })

  output$data_summary <- shiny::renderPrint({
    shiny::req(rv$se)
    cat("SummarizedExperiment Object\n")
    cat("===========================\n\n")
    cat("Dimensions:\n")
    cat("  Genes:", nrow(rv$se), "\n")
    cat("  Samples:", ncol(rv$se), "\n\n")
    cat(
      "Assays:",
      paste(names(SummarizedExperiment::assays(rv$se)), collapse = ", "),
      "\n\n"
    )
    cat("Sample Metadata Columns:\n")
    cat(
      " ",
      paste(colnames(SummarizedExperiment::colData(rv$se)), collapse = ", "),
      "\n\n"
    )
    if (ncol(SummarizedExperiment::rowData(rv$se)) > 0) {
      cat("Gene Metadata Columns:\n")
      cat(
        " ",
        paste(colnames(SummarizedExperiment::rowData(rv$se)), collapse = ", "),
        "\n"
      )
    }
  })

  output$metadata_table <- DT::renderDT({
    shiny::req(rv$se)
    DT::datatable(
      as.data.frame(SummarizedExperiment::colData(rv$se)),
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = TRUE
    )
  })

  output$pca_plot <- plotly::renderPlotly({
    shiny::req(pca_data())

    pca_info <- pca_data()
    df <- pca_info$data
    var <- pca_info$var

    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = PC1, y = PC2, color = group, text = sample)
    ) +
      ggplot2::geom_point(size = 4, alpha = 0.8) +
      ggplot2::labs(
        x = paste0("PC1 (", var[1], "%)"),
        y = paste0("PC2 (", var[2], "%)"),
        color = rv$color_var
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(legend.position = "right")

    plotly::ggplotly(p, tooltip = c("text", "group")) %>%
      plotly::layout(hovermode = "closest")
  })

  output$pca_variance <- shiny::renderPrint({
    shiny::req(rv$pca_result)
    var_explained <- round(
      100 * rv$pca_result$sdev^2 / sum(rv$pca_result$sdev^2),
      2
    )
    cat("Variance Explained by PCs:\n")
    for (i in seq_len(min(10, length(var_explained)))) {
      cat(sprintf("  PC%d: %.2f%%\n", i, var_explained[i]))
    }
  })
  ## expression plot-------
  output$expr_plot <- plotly::renderPlotly({
    shiny::req(rv$se, input$gene_id, rv$color_var)

    gene <- input$gene_id
    counts_data <- SummarizedExperiment::assay(rv$se, "counts")[gene, ]

    plot_df <- data.frame(
      s_a_m_p_l_e = colnames(rv$se),
      expression = counts_data,
      group = forcats::as_factor(SummarizedExperiment::colData(rv$se)[[
        rv$color_var
      ]])
    ) %>%
      dplyr::filter(group %in% input$groups_to_show)

    gene_label <- if (
      "gene_name" %in% colnames(SummarizedExperiment::rowData(rv$se))
    ) {
      SummarizedExperiment::rowData(rv$se)[gene, "gene_name"]
    } else {
      gene
    }

    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = group, y = expression, fill = group, text = s_a_m_p_l_e)
    ) +
      ggplot2::labs(
        title = paste("Expression:", gene_label),
        x = rv$color_var,
        y = "Normalized Counts"
      ) +
      ggplot2::scale_x_discrete(labels = \(x) {
        stringr::str_wrap(stringr::str_replace_all(x, "_", " "), 10)
      }) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(legend.position = "none")

    if (input$plot_type == "box") {
      p <- p +
        ggplot2::geom_boxplot(alpha = 0.7) +
        ggplot2::geom_jitter(width = 0.2, alpha = 0.5, size = 2)
    } else {
      p <- p +
        ggplot2::geom_violin(alpha = 0.7) +
        ggplot2::geom_jitter(width = 0.1, alpha = 0.5, size = 2)
    }

    plotly::ggplotly(p)
  })

  output$expr_table <- DT::renderDT({
    shiny::req(rv$se, input$gene_id, input$groups_to_show)

    gene <- input$gene_id
    expr_data <- data.frame(
      Sample = colnames(rv$se),
      SummarizedExperiment::colData(rv$se),
      Gene = input$gene_id,
      Count = SummarizedExperiment::assay(rv$se, "counts")[gene, ],
      g_r_o_u_p = SummarizedExperiment::colData(rv$se)[[rv$color_var]]
    ) %>%
      dplyr::filter(g_r_o_u_p %in% input$groups_to_show) %>%
      dplyr::select(-g_r_o_u_p)

    DT::datatable(
      expr_data,
      options = list(pageLength = 12, scrollX = TRUE),
      rownames = FALSE
    ) %>%
      DT::formatRound("Count", digits = 0)
  })

  output$de_table <- DT::renderDT({
    de_data <- current_de_results()
    de_data <- dplyr::select(
      de_data,
      tidyselect::any_of(c(
        "gene_id",
        "baseMean",
        "log2FoldChange",
        "pvalue",
        "padj"
      ))
    ) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::any_of(c("baseMean", "log2FoldChange")),
          \(x) round(x, 2)
        ),
        dplyr::across(tidyselect::any_of(c("pvalue", "padj")), \(x) round(x, 4))
      )
    shiny::req(de_data)

    DT::datatable(
      de_data,
      options = list(pageLength = 25, scrollX = TRUE),
      rownames = FALSE,
      filter = "top"
    )
  })

  output$de_plot <- shiny::renderPlot({
    # shiny::req(de_data, input$padj_cutoff_volcano, input$lfc_cutoff_volcano, rv$up_col, rv$dn_col, rv$highlight_col)
    colors_acute <- c(
      "Up" = rv$up_col,
      "Down" = rv$dn_col,
      "NS" = "grey70",
      "Highlighted" = rv$highlight_col
    )
    de_data <- current_de_results()

    .plot_des(
      de_data,
      input$de_comparison,
      input$padj_cutoff_volcano,
      input$lfc_cutoff_volcano,
      colors_acute
    )
  })

  output$download_expr <- shiny::downloadHandler(
    filename = function() {
      comparison_name <- stringr::str_replace(
        input$gene_id,
        "[^a-zA-Z0-9_-]",
        "_"
      )

      paste0(comparison_name, "_expression_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      shiny::req(rv$se, input$gene_id, input$groups_to_show)

      gene <- input$gene_id
      expr_data <- data.frame(
        Sample = colnames(rv$se),
        SummarizedExperiment::colData(rv$se),
        Gene = input$gene_id,
        Count = SummarizedExperiment::assay(rv$se, "counts")[gene, ],
        g_r_o_u_p = SummarizedExperiment::colData(rv$se)[[rv$color_var]]
      ) %>%
        dplyr::filter(g_r_o_u_p %in% input$groups_to_show) %>%
        dplyr::select(-g_r_o_u_p)
      readr::write_excel_csv2(expr_data, file)
    }
  )

  output$download_de <- shiny::downloadHandler(
    filename = function() {
      comparison_name <- if (
        has_precomputed_de() && !is.null(input$de_comparison)
      ) {
        gsub("[^a-zA-Z0-9_-]", "_", input$de_comparison)
      } else {
        "DE_results"
      }
      paste0(comparison_name, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      de_data <- current_de_results()
      shiny::req(de_data)
      readr::write_excel_csv2(de_data, file)
    }
  )

  ## Volcano plot ----
  output$volcano_plot <- plotly::renderPlotly({
    de_data <- current_de_results()
    shiny::req(
      de_data,
      input$padj_cutoff_volcano,
      input$lfc_cutoff_volcano,
      rv$up_col,
      rv$dn_col,
      rv$highlight_col
    )
    colors_acute <- c(
      "Up" = rv$up_col,
      "Down" = rv$dn_col,
      "NS" = "grey70",
      "Highlighted" = rv$highlight_col
    )

    highlight_vec <- c()
    if (
      !is.null(input$highlight_genes) &&
        nchar(trimws(input$highlight_genes)) > 0
    ) {
      # Split by newlines and commas, trim whitespace
      highlight_vec <- input$highlight_genes %>%
        strsplit("[\n,]") %>%
        unlist() %>%
        trimws() %>%
        .[nchar(.) > 0]
    }

    .plot_volcano(
      RES = de_data,
      NAME = input$de_comparison,
      padj_CO = input$padj_cutoff_volcano,
      fc_CO = input$lfc_cutoff_volcano,
      highlights = highlight_vec,
      COLS = colors_acute,
      LABEL_TOP = input$label_top,
      TOPN = input$n_labels
    )
  })
  output$volcano_summary <- shiny::renderPrint({
    de_data <- current_de_results()
    shiny::req(
      de_data,
      input$padj_cutoff_volcano,
      input$lfc_cutoff_volcano,
      rv$up_col,
      rv$dn_col,
      rv$highlight_col
    )

    colors_acute <- c(
      "Up" = rv$up_col,
      "Down" = rv$dn_col,
      "NS" = "grey70",
      "Highlighted" = rv$highlight_col
    )
    summary_df <- de_data %>%
      dplyr::filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
      dplyr::summarise(
        total_genes = dplyr::n(),
        sig_genes = sum(padj < input$padj_cutoff_volcano),
        up_regulated = sum(
          padj < input$padj_cutoff_volcano &
            log2FoldChange > input$lfc_cutoff_volcano
        ),
        down_regulated = sum(
          padj < input$padj_cutoff_volcano &
            log2FoldChange < -input$lfc_cutoff_volcano
        ),
        not_significant = sum(
          padj >= input$padj_cutoff_volcano |
            (abs(log2FoldChange) < input$lfc_cutoff_volcano &
              padj < input$padj_cutoff_volcano)
        )
      )

    cat("Differential Expression Summary\n")
    cat("================================\n\n")
    cat(sprintf("Total genes tested: %d\n", summary_df$total_genes))
    cat(sprintf(
      "Significant (padj < %.3f): %d (%.1f%%)\n",
      input$padj_cutoff_volcano,
      summary_df$sig_genes,
      100 * summary_df$sig_genes / summary_df$total_genes
    ))
    cat(sprintf(
      "\nUp-regulated (LFC > %.2f): %d\n",
      input$lfc_cutoff_volcano,
      summary_df$up_regulated
    ))
    cat(sprintf(
      "Down-regulated (LFC < -%.2f): %d\n",
      input$lfc_cutoff_volcano,
      summary_df$down_regulated
    ))
    cat(sprintf("Not significant: %d\n", summary_df$not_significant))
  })
  ## enrichment plots -------
  output$enrichment_plots <- shiny::renderUI({
    current_fes <- current_fe_results()
    colors_acute <- c(
      "Up" = rv$up_col,
      "Down" = rv$dn_col,
      "NS" = "grey70",
      "Highlighted" = rv$highlight_col
    )
    shiny::req(
      current_fes,
      input$padj_cutoff_enrichment,
      input$n_terms_enrichment,
      colors_acute
    )

    plots <- purrr::imap(current_fes, function(enrich, name) {
      id <- paste0("e_plot_", name)

      output[[id]] <- plotly::renderPlotly({
        .plot_fe(
          FE = enrich,
          NAME = name,
          padj_CO = input$padj_cutoff_enrichment,
          N_terms = input$n_terms_enrichment,
          COLS = colors_acute
        )
      })
    })
  })
}

# Run app
#' exploreSE
#'
#'
#' @description
#' this runs the explorer app
#'
#' @param file a string leading to a file you want loaded as a default
#' @param object a summarizedExperiment object (or DeeDeeExperiment) to load in
#'
#' @returns an application
#' @export
#'
#' @examples
#' app <- exploreSE()
#' if(interactive()){
#' shiny::runApp(app, port = 1234)
#'}

exploreSE <- function(file = NULL, object = NULL) {
  shiny::shinyOptions(set_file_name = file, set_object = object)
  shiny::shinyApp(ui, server)
}
