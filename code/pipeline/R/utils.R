`%+%` <- paste0
`%_%` <- paste
nl    <- function (ll, n = 1) ll %+% Reduce(paste0, rep('\n', n))

die <- function (msg, status = 1) {
  cat(msg)
  quit(save = 'no', status = status)
}

silent_lib_loader <- function (...) {
  libs <- c(...)
  for (l in libs) {
    suppressMessages(library(l, quietly = TRUE, character.only = TRUE,
                             warn.conflicts = FALSE))
  }
}

parse_arg <- function (arguments, arg_name, fun) {
  argument <- getElement(arguments, arg_name)
  tryCatch(fun(argument), error = function (e) {
    die(sprintf('Error while trying to validate "%s" argument:\n\t%s\n',
                arg_name, e$message))
  })
}

parse_arg_pred <- function (arguments, arg_name, pred, error, converter = identity) {
  parse_arg(arguments, arg_name, function (x) {
              res <- converter(x)
              if (pred(res)) {
                res
              } else {
                msg <- error
                if (is.function(error)) {
                  msg <- error(res)
                }
                stop(msg)
              }
  })
}

check_arg_in_list <- function (...) {
  function (aa) !anyNA(pmatch(aa, c(...), duplicates.ok = TRUE))
}

check_if_positive <- function (x) !is.na(x) && x > 0
conv_to_numeric   <- function (x) suppressWarnings(as.numeric(x))
conv_to_integer   <- function (x) suppressWarnings(as.integer(x))

make_main <- function (docstring, libs = c(), main_function = function () {}) {
  main_template <- function (...) {
    do.call(silent_lib_loader, as.list(c('docopt', libs)))
    main_env <- environment(main_function)
    main_env$cli_args <- docopt::docopt(docstring)
    main_function(...)
  }

  function (...) {
    tryCatch(main_template(...), error = function (e) {
      die(e$message)
    })
  }
}

