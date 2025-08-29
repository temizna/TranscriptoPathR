# R/util_contrasts.R

#' Sanitize a string for filenames/tags
#'
#' Replace any character that is not A–Z, a–z, 0–9, dot, underscore, or dash
#' with an underscore.
#'
#' @param x Character vector to sanitize.
#' @return Character vector with only [A-Za-z0-9._-].
#' @keywords internal
#' @noRd
.safe_tag <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

#' Build a "test_vs_ref" tag from DE selection reactives
#'
#' Safely extracts ref_level() and test_level() from a list returned by
#' mod_de_server() and returns "Test_vs_Ref" with unsafe characters replaced.
#' Falls back to "contrast" if the levels are not available.
#'
#' @param de_sel A list returned by mod_de_server(), containing
#'   group_var, ref_level, and test_level reactives.
#' @return A length-1 character tag such as "TreatmentB_vs_TreatmentA".
#' @keywords internal
#' @noRd
.contrast_tag_from <- function(de_sel) {
  if (is.null(de_sel)) return("contrast")
  ref  <- try(de_sel$ref_level(),  silent = TRUE)
  test <- try(de_sel$test_level(), silent = TRUE)
  if (inherits(ref, "try-error") || inherits(test, "try-error") ||
      is.null(ref) || is.null(test) || ref == "" || test == "") {
    return("contrast")
  }
  paste0(.safe_tag(test), "_vs_", .safe_tag(ref))
}

#' Comparison bridge for filenames and display labels
#'
#' Wrap the list returned by mod_de_server() and add two convenience
#' reactives:
#' - label(): human-readable label built from label_fmt.
#' - tag(): file-safe tag built as "Test_vs_Ref".
#'
#' @param de_sel A list returned by mod_de_server(), containing
#'   group_var, ref_level, and test_level reactives.
#' @param label_fmt Character scalar. An sprintf-style format with two %s
#'   placeholders (first = test, second = ref). Default: "%s vs %s".
#'
#' @return A list of reactives: group_var, ref_level, test_level, label, tag.
#'
#' @examples
#' \dontrun{
#'   de_sel <- mod_de_server("de", filtered_data_rv, filtered_dds_rv, res_reactive)
#'   cmp <- make_cmp_bridge(de_sel)
#'   cmp$label()  # "TreatmentB vs TreatmentA"
#'   cmp$tag()    # "TreatmentB_vs_TreatmentA"
#' }
#' @importFrom shiny reactive
#' @export
make_cmp_bridge <- function(de_sel, label_fmt = "%s vs %s") {
  stopifnot(!is.null(de_sel))
  list(
    group_var  = de_sel$group_var,
    ref_level  = de_sel$ref_level,
    test_level = de_sel$test_level,
    label = shiny::reactive({
      sprintf(
        label_fmt,
        tryCatch(de_sel$test_level(), error = function(e) "TEST"),
        tryCatch(de_sel$ref_level(),  error = function(e) "REF")
      )
    }),
    tag = shiny::reactive({
      .contrast_tag_from(de_sel)
    })
  )
}

#' Get grouping vector for a set of samples, with comparison-aware fallback
#' @param de_sel list returned by mod_de_server()
#' @param samples_df data.frame of sample metadata (rownames = sample IDs)
#' @param sample_ids character vector of sample IDs to align to
#' @return character vector of groups aligned to sample_ids
#' @keywords internal
de_group_vector <- function(de_sel, samples_df, sample_ids) {
  stopifnot(!is.null(samples_df), !is.null(sample_ids))
  # comparison-driven DE?
  if (!is.null(de_sel$group_var) &&
      identical(tryCatch(de_sel$group_var(), error = function(e) NULL), "cmp_group") &&
      !is.null(de_sel$cmp_factor) && !is.null(de_sel$included_samples)) {
    cf  <- tryCatch(as.character(de_sel$cmp_factor()), error = function(e) NULL)
    ids <- tryCatch(de_sel$included_samples(),        error = function(e) NULL)
    if (!is.null(cf) && !is.null(ids) && length(cf) == length(ids)) {
      m <- setNames(cf, ids)
      out <- m[sample_ids]
      return(unname(out))
    }
  }
  # regular DE (no comparison) -> read from metadata column
  gv <- tryCatch(de_sel$group_var(), error = function(e) NULL)
  if (!is.null(gv) && gv %in% colnames(samples_df)) {
    v <- as.character(samples_df[sample_ids, gv, drop = TRUE])
    return(unname(v))
  }
  rep(NA_character_, length(sample_ids))
}

