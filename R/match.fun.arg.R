# Matches up function argument using partial matching
# Observe that if parent call was is a named function this match.fun.arg will
# produces garbage stop message!
# arg (symbol(?)): name of a formal argument or parent function
# multiple.ok (boolean, atomic): are multiple matched ok or not
# msg (character, ?): header of message string in messages issued by this
#   function, otherwise the name of the calling function is used (if unnmaed
#   function this is a mess and hence the option of supplying this string)
match.fun.arg <- function (arg, multiple.ok = FALSE, msg) {

  cl <- sys.call(sys.parent()) # 1st element here big is not names function
  if ( missing(msg) ) {  msg <- deparse(cl[[1]])  }
  formal.args <- formals(sys.function(sys.parent()))
  choices <- eval(formal.args[[deparse(substitute(arg))]])
  i <- pmatch(arg, choices)
  if ( is.na(i) ) {
    stop (
      sprintf(
        "%s: unable to match arg '%s' in '%s'",
        msg,
        arg,
        paste(choices,collapse="', '")
      )
    )
  } else if ( !multiple.ok && length(i)!=1 ) {
    stop (sprintf("%s: illegal multiple matches for arg",msg))
  } else {
    choices[i]
  }

}
