# R/theme_ios.R

#' iOS-like clean theme (flat + crisp)
#' @return a bslib theme
#' @export
tpr_theme <- function() {
  bslib::bs_theme(
    version = 5,
    
    # Fonts (no font_face, so nothing to lazily load)
    base_font    = bslib::font_google("Inter"),
    heading_font = bslib::font_google("Inter"),
    code_font    = bslib::font_google("JetBrains Mono"),
    
    # iOS-ish palette
    primary = "#007AFF",  # iOS system blue
    success = "#34C759",
    info    = "#5AC8FA",
    warning = "#FF9500",
    danger  = "#FF3B30",
    
    # Surfaces & text
    "body-bg"    = "#FFFFFF",
    "body-color" = "#111111",
    
    # Shape / spacing
    "border-radius"       = "0.75rem",
    "btn-border-radius"   = "0.60rem",
    "input-border-radius" = "0.55rem",
    "input-btn-padding-y" = ".55rem",
    "input-btn-padding-x" = ".90rem",
    
    # Subtle size tweak (optional: slightly smaller base text)
    "font-size-base" = "0.9375rem",  # 15px on 16px root
    
    # Flatten shadows
    "enable-shadows"   = FALSE,
    "box-shadow"       = "none",
    "btn-box-shadow"   = "none",
    "input-box-shadow" = "none"
  ) |>
    bslib::bs_add_variables(
      # Hairline dividers
      "border-color"            = "#E5E5EA",
      "card-border-color"       = "#E5E5EA",
      "input-border-color"      = "#D1D5DB",
      "input-focus-border-color"= "#007AFF",
      "input-focus-box-shadow"  = "none"
    ) |>
    bslib::bs_add_rules("
      /* Scope to app root */
      .tpr-app .nav-tabs { border-bottom: 0; }
      .tpr-app .nav-tabs .nav-link {
        border: 0; background: transparent; color: #6B7280;
        padding: .55rem .85rem; font-weight: 500;
      }
      .tpr-app .nav-tabs .nav-link.active {
        color: #111111; border-bottom: 2px solid var(--bs-primary);
      }

      /* Sub-tabs inside plot panels get a coordinating teal */
      .tpr-app .tpr-plot-tabs .nav-tabs .nav-link.active {
        border-bottom-color: #14B8A6; /* teal-500 */
      }

      /* Sidebar: white + hairline */
      .tpr-app .bslib-sidebar-layout .sidebar {
        background: #FFFFFF; border-right: 1px solid #E5E5EA; padding-top: .75rem;
        box-shadow: none;
      }

      /* Flat surfaces */
      .tpr-surface {
        background: #FFFFFF; border: 1px solid #E5E5EA; border-radius: 12px;
        box-shadow: none; padding: 12px; margin-bottom: 16px;
      }

      /* Inputs / buttons: crisp focus, no shadows */
      .tpr-app .form-control, .tpr-app .form-select,
      .tpr-app .btn, .tpr-app .card, .tpr-app .dropdown-menu, .tpr-app .modal-content {
        box-shadow: none !important;
      }
      .tpr-app .form-control, .tpr-app .form-select {
        background: #FFFFFF !important; border: 1px solid #D1D5DB !important;
      }
      .tpr-app .form-control:focus, .tpr-app .form-select:focus {
        border-color: #007AFF !important; box-shadow: none !important;
      }
      .tpr-app .btn { border: 1px solid #D1D5DB; }

      /* File input progress flat */
      .tpr-app .progress, .tpr-app .progress-bar { box-shadow: none !important; }

      /* Compact spacing */
      .tpr-app .shiny-input-container { margin-bottom: 0.75rem; }
      .tpr-app .help-block { color: #6B7280; }
    ")
}
