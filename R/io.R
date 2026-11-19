# Internal image I/O and class-preservation helpers for pTFCE
# Supports RNifti::niftiImage, oro.nifti S4 objects, plain R arrays, and file paths.

is_missing_arg <- function(x) {
  is.null(x) || (is.atomic(x) && length(x) == 1L && is.na(x))
}

read_volume_input <- function(x, arg_name = "img", allow_4d = FALSE) {
  kind <- "array"
  template <- NULL

  if (is.character(x)) {
    if (length(x) != 1L || is.na(x) || !nzchar(x)) {
      stop(sprintf("`%s` file path must be a non-empty character scalar.", arg_name), call. = FALSE)
    }
    if (!file.exists(x)) {
      stop(sprintf("`%s` file does not exist: %s", arg_name, x), call. = FALSE)
    }
    nii <- RNifti::readNifti(x)
    kind <- "rnifti"
    template <- nii
    arr <- as.array(nii)
  } else if (inherits(x, "niftiImage")) {
    kind <- "rnifti"
    template <- x
    arr <- as.array(x)
  } else if (isS4(x) && (methods::is(x, "nifti") || methods::is(x, "anlz"))) {
    kind <- "oro"
    template <- x
    arr <- as.array(x@.Data)
  } else if (is.numeric(x) || is.logical(x)) {
    arr <- as.array(x)
    if (!is.null(dim(x))) {
      dim(arr) <- dim(x)
    }
  } else {
    stop(
      sprintf("Unsupported `%s` type (%s). Supply an RNifti image, oro.nifti object, array, or NIfTI file path.",
              arg_name, paste(class(x), collapse = "/")),
      call. = FALSE
    )
  }

  orig_dim <- dim(arr)
  if (is.null(orig_dim) || length(orig_dim) < 2L) {
    stop(sprintf("`%s` must be at least 2-dimensional.", arg_name), call. = FALSE)
  }

  storage.mode(arr) <- "double"

  # Normalize dimensions to (nx, ny, nz, nt) while recording original shape
  if (length(orig_dim) == 2L) {
    dim4 <- c(orig_dim[1L], orig_dim[2L], 1L, 1L)
  } else if (length(orig_dim) == 3L) {
    dim4 <- c(orig_dim[1L], orig_dim[2L], orig_dim[3L], 1L)
  } else if (length(orig_dim) >= 4L) {
    trailing <- prod(orig_dim[-(1:3)])
    if (!allow_4d && trailing > 1L) {
      stop(sprintf("`%s` must be a 2D or 3D spatial image (got dimensions %s).",
                   arg_name, paste(orig_dim, collapse = "x")), call. = FALSE)
    }
    dim4 <- c(orig_dim[1L], orig_dim[2L], orig_dim[3L], as.integer(trailing))
  }

  dim(arr) <- if (allow_4d) dim4 else dim4[1:3]

  list(
    data = arr,
    orig_dim = orig_dim,
    dim3 = dim4[1:3],
    nt = dim4[4L],
    kind = kind,
    template = template
  )
}

wrap_volume_output <- function(arr, meta) {
  # Restore original spatial dimensionality (e.g. 2D, 3D, or 3D-with-trailing-1)
  out_dim <- if (length(meta$orig_dim) == 2L) meta$orig_dim else meta$dim3
  if (length(meta$orig_dim) == 4L && meta$orig_dim[4L] == 1L) {
    out_dim <- meta$orig_dim
  }
  dim(arr) <- out_dim

  if (identical(meta$kind, "rnifti")) {
    return(RNifti::updateNifti(arr, template = meta$template))
  }

  if (identical(meta$kind, "oro")) {
    if (!requireNamespace("oro.nifti", quietly = TRUE)) {
      return(arr)
    }
    img_tpl <- meta$template
    snames <- methods::slotNames(img_tpl)
    snames <- setdiff(snames, c(".Data", "dim_", "cal_max", "cal_min", "scl_slope", "scl_inter"))
    out <- oro.nifti::nifti(img = arr, dim = dim(arr))
    class(out) <- class(img_tpl)
    for (islot in snames) {
      methods::slot(out, islot) <- methods::slot(img_tpl, islot)
    }
    oro.nifti::datatype(out) <- 16L # FLOAT32
    oro.nifti::bitpix(out) <- 32L   # FLOAT32
    out@scl_slope <- 0
    out@scl_inter <- 0
    out <- oro.nifti::zero_trans(out)
    r <- range(arr[is.finite(arr)], na.rm = TRUE)
    if (length(r) == 2L && all(is.finite(r))) {
      out@cal_min <- r[1L]
      out@cal_max <- r[2L]
    }
    return(out)
  }

  arr
}
