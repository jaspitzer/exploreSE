# # se slot checkers and getters

# .react_de <- function(input, session, rv) {
#   shiny::reactive({
#     shiny::req(rv$se)

#     .check_precomputed_de(rv$se)
#   })
# }

# .react_fe <- function(input, session, rv) {
#   shiny::reactive({
#     shiny::req(rv$se)
#     .check_precomputed_fe(rv$se)
#   })
# }

# .react_de_names <- function(input, session, rv) {
#   shiny::reactive({
#     shiny::req(has_precomputed_de())
#     .de_results_names(rv$se)
#   })
# }

# .react_fe_names <- function(input, session, rv) {
#   shiny::reactive({
#     shiny::req(has_precomputed_fe())
#     .fe_results_names(rv$se, input$de_comparison)
#   })
# }
