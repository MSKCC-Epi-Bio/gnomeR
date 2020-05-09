#' @importFrom rlang .data .env %||% set_names sym syms parse_expr expr exprs call2 :=
#' @importFrom dplyr mutate select n group_by ungroup filter pull case_when
#' if_else full_join left_join distinct bind_rows count coalesce arrange rename
#' rename_at bind_cols mutate_all mutate_at slice desc
#' @importFrom utils tail
#' @keywords internal
"_PACKAGE"

# allowing for the use of the dot when piping
utils::globalVariables(".")

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
