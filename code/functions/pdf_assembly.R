pdf_panel <- function(
  id, x, y, width, label_x = x - 1.5, label_y = y, source = NULL, crop = TRUE
) {
  list(
    id = id,
    x = x,
    y = y,
    width = width,
    label_x = label_x,
    label_y = label_y,
    source = source,
    crop = crop
  )
}

assemble_pdf_figure <- function(
  stem, page, panels, repo, figure_root = file.path(repo, "figures")
) {
  current_work <- tempfile("figure_composite_")
  on.exit(unlink(current_work, recursive = TRUE), add = TRUE)
  cropped_root <- file.path(current_work, "cropped")
  dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(cropped_root, recursive = TRUE, showWarnings = FALSE)

  required_programs <- c(pdfcrop = "pdfcrop", xelatex = "xelatex")
  tex_bins <- c("/Library/TeX/texbin", path.expand("~/.local/TinyTeX/bin/x86_64-linux"))
  tex_bins <- tex_bins[dir.exists(tex_bins)]
  program_paths <- vapply(required_programs, function(program) {
    candidates <- c(unname(Sys.which(program)), file.path(tex_bins, program))
    candidates <- candidates[nzchar(candidates) & file.access(candidates, mode = 1L) == 0L]
    if (length(candidates)) candidates[[1L]] else ""
  }, character(1))
  if (any(!nzchar(program_paths))) {
    stop(
      "Figure assembly requires: ",
      paste(names(program_paths)[!nzchar(program_paths)], collapse = ", "),
      call. = FALSE
    )
  }

  run_checked <- function(command, args) {
    # pdfcrop writes tmp-pdfcrop-*.tex to TMPDIR; point explicitly to current_work
    # to ensure cleanup on exit and prevent cluttering the repository root.
    output <- system2(command, args = args, stdout = TRUE, stderr = TRUE,
                      env = c(paste0("TMPDIR=", shQuote(current_work)),
                              paste0("PATH=", shQuote(paste(
                                c(tex_bins, Sys.getenv("PATH")), collapse = .Platform$path.sep
                              )))))
    status <- attr(output, "status")
    if (!is.null(status) && status != 0L) {
      stop(
        "Command failed: ", command, " ", paste(args, collapse = " "),
        "\n", paste(output, collapse = "\n"),
        call. = FALSE
      )
    }
    invisible(output)
  }

  panel_ids <- vapply(panels, `[[`, character(1), "id")
  if (anyDuplicated(panel_ids)) {
    stop("Duplicate panel ids in ", stem, call. = FALSE)
  }
  source_files <- vapply(
    panels,
    function(current) {
      if (is.null(current$source)) paste0(stem, current$id, ".pdf") else current$source
    },
    character(1)
  )
  source_paths <- setNames(file.path(figure_root, source_files), panel_ids)

  cropped_paths <- file.path(cropped_root, basename(source_paths))
  names(cropped_paths) <- names(source_paths)
  for (id in names(source_paths)) {
    if (isFALSE(panels[[match(id, panel_ids)]]$crop)) {
      if (!file.copy(source_paths[[id]], cropped_paths[[id]], overwrite = TRUE)) {
        stop("Unable to copy panel: ", source_paths[[id]], call. = FALSE)
      }
      next
    }
    run_checked(
      program_paths[["pdfcrop"]],
      c(
        "--margins", shQuote("3 3 3 3"),
        shQuote(source_paths[[id]]), shQuote(cropped_paths[[id]])
      )
    )
  }

  tex_path <- function(path) normalizePath(path, winslash = "/", mustWork = TRUE)
  page_width <- page[[1L]]
  page_height <- page[[2L]]
  nodes <- unlist(lapply(panels, function(current) {
    c(
      sprintf(
        paste0(
          "  \\node[anchor=north west,inner sep=0pt] at (%s,%s) ",
          "{\\includegraphics[width=%smm]{\\detokenize{%s}}};"
        ),
        current$x, current$y, current$width,
        tex_path(cropped_paths[[current$id]])
      ),
      sprintf(
        "  \\panellabel{%s}{%s}{%s}",
        current$label_x, current$label_y, current$id
      )
    )
  }))

  tex_lines <- c(
    "\\documentclass{article}",
    sprintf(
      "\\usepackage[paperwidth=%smm,paperheight=%smm,margin=0mm]{geometry}",
      page_width, page_height
    ),
    "\\usepackage{fontspec}",
    "\\setsansfont{Arial}",
    "\\usepackage{graphicx}",
    "\\usepackage{tikz}",
    "\\pagestyle{empty}",
    "\\setlength{\\parindent}{0pt}",
    "\\newcommand{\\panellabel}[3]{%",
    paste0(
      "  \\node[anchor=north west,inner sep=0pt,font=\\sffamily",
      "\\bfseries\\fontsize{10}{10}\\selectfont] at (#1,#2) {#3};%"
    ),
    "}",
    "\\begin{document}",
    "\\noindent",
    "\\begin{tikzpicture}[x=1mm,y=1mm]",
    sprintf(
      "  \\useasboundingbox (0,0) rectangle (%s,%s);",
      page_width, page_height
    ),
    nodes,
    "\\end{tikzpicture}",
    "\\end{document}"
  )
  tex_file <- file.path(current_work, paste0(stem, ".tex"))
  writeLines(tex_lines, tex_file, useBytes = TRUE)
  run_checked(
    program_paths[["xelatex"]],
    c(
      "-halt-on-error", "-interaction=nonstopmode",
      paste0("-output-directory=", shQuote(current_work)), shQuote(tex_file)
    )
  )

  compiled_pdf <- file.path(current_work, paste0(stem, ".pdf"))
  output_pdf <- file.path(figure_root, paste0(stem, ".pdf"))
  run_checked(
    program_paths[["pdfcrop"]],
    c(
      "--margins", shQuote("6 6 6 6"),
      shQuote(compiled_pdf), shQuote(output_pdf)
    )
  )
  message("Saved assembled ", stem, " to ", output_pdf)
  invisible(output_pdf)
}
