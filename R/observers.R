.create_dir_color_observers <- function(input, session, rv) {
  shiny::observe({
    # for up
    shiny::updateSelectInput(session, "up_col_1", selected = rv$up_col)
    shiny::updateSelectInput(session, "up_col_2", selected = rv$up_col)
    # down
    shiny::updateSelectInput(session, "dn_col_1", selected = rv$dn_col)
    shiny::updateSelectInput(session, "dn_col_2", selected = rv$dn_col)
    # highlight
    shiny::updateSelectInput(session, "high_col_1", selected = rv$highlight_col)
    shiny::updateSelectInput(session, "high_col_2", selected = rv$highlight_col)
  })

  shiny::observeEvent(input$up_col_1, {
    rv$up_col <- input$up_col_1
  })
  shiny::observeEvent(input$up_col_2, {
    rv$up_col <- input$up_col_2
  })
  #down
  shiny::observeEvent(input$dn_col_1, {
    rv$dn_col <- input$dn_col_1
  })
  shiny::observeEvent(input$dn_col_2, {
    rv$dn_col <- input$dn_col_2
  })
  #highlight
  shiny::observeEvent(input$high_col_1, {
    rv$highlight_col <- input$high_col_1
  })
  shiny::observeEvent(input$high_col_2, {
    rv$highlight_col <- input$high_col_2
  })
}

.create_interest_color_observers <- function(input, session, rv) {
  shiny::observe({
    # color var
    shiny::updateSelectInput(session, "color_var_1", selected = rv$color_var)
    shiny::updateSelectInput(session, "color_var_2", selected = rv$color_var)

    if (!is.null(rv$se)) {
      shiny::updateCheckboxGroupInput(
        session,
        "groups_to_show",
        choices = levels(SummarizedExperiment::colData(rv$se)[[rv$color_var]]),
        selected = levels(SummarizedExperiment::colData(rv$se)[[rv$color_var]])
      )
    }
  })

  shiny::observeEvent(input$color_var_1, {
    rv$color_var <- input$color_var_1
  })
  shiny::observeEvent(input$color_var_2, {
    rv$color_var <- input$color_var_2
  })

  # shiny::observeEvent(rv$color_var, {
  #   if (!is.null(rv$se)) {
  #     shiny::updateCheckboxGroupInput(
  #       session,
  #       "groups_to_show",
  #       choices = levels(as.factor(SummarizedExperiment::colData(rv$se)[[
  #         rv$color_var
  #       ]])),
  #       selected = levels(as.factor(SummarizedExperiment::colData(rv$se)[[
  #         rv$color_var
  #       ]]))
  #     )
  #   }
  # })
}

.create_rowvar_observer <- function(input, session, rv) {
  shiny::observeEvent(input$row_data_var, {
    shiny::req(input$row_data_var)
    rv$row_var <- input$row_data_var
  })

  shiny::observeEvent(rv$row_var, {
    if (
      !is.null(rv$se) &&
        !is.null(SummarizedExperiment::rowData(rv$se)) &&
        rv$row_var %in% colnames(SummarizedExperiment::rowData(rv$se))
    ) {
      gene_choices <- SummarizedExperiment::rowData(rv$se)[,
        rv$row_var,
        drop = T
      ]
      # if ("gene_name" %in% colnames(SummarizedExperiment::rowData(rv$se))) {
      #   names(gene_choices) <- SummarizedExperiment::rowData(rv$se)$gene_name
      # }

      shinyWidgets::updatePickerInput(
        session,
        "gene_id",
        choices = gene_choices,
        selected = gene_choices[1]
      )
    }
  })
}

.gene_ident_observer <- function(input, session, rv) {}

.observe_demo_data <- function(input, session, rv) {
  shiny::observe({
    if (input$use_demo) {
      rv$se <- .create_demo_data()
    }
  })
}

.observe_inital_obj <- function(input, session, rv) {
  set_file <- shiny::getShinyOption("set_file_name")
  set_object <- shiny::getShinyOption("set_object")

  shiny::observe({
    shiny::req(set_file)
    rv$se <- readRDS(set_file)
    shiny::updateCheckboxInput(session, "use_demo", value = FALSE)
    shiny::showNotification("Data loaded successfully!", type = "message")
  })

  shiny::observe({
    shiny::req(set_object)
    rv$se <- set_object
    shiny::updateCheckboxInput(session, "use_demo", value = FALSE)
    shiny::showNotification("Data loaded successfully!", type = "message")
  })
}


.observe_se_load <- function(
  input,
  session,
  rv,
  has_precomputed_de,
  de_comparisons
) {
  shiny::observeEvent(rv$se, {
    shiny::req(rv$se)

    # Update grouping variable choices
    col_vars <- colnames(SummarizedExperiment::colData(rv$se))
    row_vars <- colnames(SummarizedExperiment::rowData(rv$se))

    # Set default to "condition" if it exists, otherwise first column
    default_col_var <- if ("condition" %in% col_vars) {
      "condition"
    } else {
      col_vars[1]
    }

    shiny::updateSelectInput(
      session,
      "color_var_1",
      choices = col_vars,
      selected = default_col_var
    )
    shiny::updateSelectInput(
      session,
      "color_var_2",
      choices = col_vars,
      selected = default_col_var
    )

    default_row_var <- if ("gene_name" %in% row_vars) {
      "gene_name"
    } else if ("SYMBOL" %in% stringr::str_to_upper(row_vars)) {
      row_vars[which(stringr::str_to_upper(row_vars) == "SYMBOL")[1]]
    } else {
      row_vars[1]
    }
    shiny::updateSelectInput(
      session,
      "row_data_var",
      choices = row_vars,
      selected = default_row_var
    )

    # Update gene choices with searchable picker
    gene_choices <- SummarizedExperiment::rowData(rv$se)[,
      default_row_var,
      drop = T
    ]
    # if ("gene_name" %in% colnames(SummarizedExperiment::rowData(rv$se))) {
    #   names(gene_choices) <- SummarizedExperiment::rowData(rv$se)$gene_name
    # }

    shinyWidgets::updatePickerInput(
      session,
      "gene_id",
      choices = gene_choices,
      selected = gene_choices[1]
    )

    # Update DE comparison choices if precomputed results exist
    if (has_precomputed_de()) {
      comparisons <- de_comparisons()
      shiny::updateSelectInput(
        session,
        "de_comparison",
        choices = comparisons,
        selected = comparisons[1]
      )
      shiny::showNotification(
        paste("Found", length(comparisons), "precomputed DE comparison(s)"),
        type = "message"
      )
    }
  })
}


.observe_load_file <- function(input, session, rv) {
  shiny::observeEvent(input$se_file, {
    shiny::req(input$se_file)
    tryCatch(
      {
        rv$se <- readRDS(input$se_file$datapath)
        shiny::updateCheckboxInput(session, "use_demo", value = FALSE)
        shiny::showNotification("Data loaded successfully!", type = "message")
      },
      error = function(e) {
        shiny::showNotification(
          paste("Error loading file:", e$message),
          type = "error"
        )
      }
    )
  })
}
