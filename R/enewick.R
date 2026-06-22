# enewick.R -----------------------------------------------------------------
#
# A self-contained extended-Newick / Rich-Newick parser for explicit reticulate
# networks. ape::read.evonet mis-parses nested/stacked (level-k) reticulations
# (leaving phantom unlabeled tips) and fails on hybrid-free trees; this parser
# handles both by construction. See docs/WORKPLAN.md (W1) and docs/DATA.md.
#
# Format: standard Newick where each hybrid node #Hk appears twice -- once as a
# "defining" occurrence carrying its child subtree, `(subtree)#Hk:len::gamma`,
# and once as a leaf placeholder, `#Hk:len::gamma`. The hybrid has two parents:
# the defining-occurrence parent (major/tree edge) and the placeholder parent
# (minor/reticulation edge). Branch fields are `:len:support:gamma`; `::gamma`
# means support is empty.

# Recursive-descent parse into a flat node table. Returns lists indexed by an
# internal preorder node id: parent (0 = root), label, length, gamma. Internal.
parse_enewick <- function(s) {
  s <- sub(";\\s*$", "", trimws(s))
  chars <- strsplit(s, "", fixed = TRUE)[[1]]
  n <- length(chars)
  env <- new.env(parent = emptyenv())
  env$pos <- 1L
  env$parent <- integer(0)
  env$label <- character(0)
  env$length <- numeric(0)
  env$gamma <- numeric(0)

  add_node <- function(parent) {
    env$parent <- c(env$parent, parent)
    env$label  <- c(env$label, NA_character_)
    env$length <- c(env$length, NA_real_)
    env$gamma  <- c(env$gamma, NA_real_)
    length(env$parent)
  }
  peek <- function() if (env$pos <= n) chars[env$pos] else ""
  delim <- c(":", ",", "(", ")")

  read_token <- function() {
    tok <- ""
    while (env$pos <= n && !(chars[env$pos] %in% delim)) {
      tok <- paste0(tok, chars[env$pos]); env$pos <- env$pos + 1L
    }
    tok
  }

  read_annotation <- function(id) {
    lab <- read_token()
    if (nzchar(lab)) env$label[id] <- lab
    fields <- character(0)
    while (peek() == ":") {
      env$pos <- env$pos + 1L            # consume ':'
      fields <- c(fields, read_token())
    }
    if (length(fields) >= 1L && nzchar(fields[1])) {
      env$length[id] <- suppressWarnings(as.numeric(fields[1]))
    }
    # gamma is the last colon field (`:len:support:gamma` or `:len::gamma`)
    if (length(fields) >= 3L && nzchar(fields[length(fields)])) {
      env$gamma[id] <- suppressWarnings(as.numeric(fields[length(fields)]))
    }
  }

  parse_clade <- function(parent) {
    id <- add_node(parent)
    if (peek() == "(") {
      env$pos <- env$pos + 1L             # consume '('
      repeat {
        parse_clade(id)
        ch <- peek()
        if (ch == ",") { env$pos <- env$pos + 1L; next }
        if (ch == ")") { env$pos <- env$pos + 1L; break }
        break                              # malformed; stop this level
      }
    }
    read_annotation(id)
    id
  }

  parse_clade(0L)
  list(parent = env$parent, label = env$label,
       length = env$length, gamma = env$gamma)
}

# Build an evonet (+ gamma table + hybrid-node map) from a parsed node table.
# Merges each hybrid's placeholder occurrence into its defining occurrence as a
# reticulation edge, then renumbers to ape conventions. Internal.
enewick_to_anansi <- function(s) {
  p <- parse_enewick(s)
  P <- p$parent; label <- p$label; len <- p$length; gam <- p$gamma
  nn <- length(P)
  if (nn == 0L) stop("empty parse")
  children <- lapply(seq_len(nn), function(i) which(P == i))

  ishyb <- !is.na(label) & startsWith(label, "#")
  remove <- logical(nn)
  ret_from <- integer(0); ret_to <- integer(0)
  rlab <- character(0); maj <- numeric(0); minr <- numeric(0)
  for (hl in unique(label[ishyb])) {
    ids <- which(label == hl)
    haskid <- vapply(ids, function(i) length(children[[i]]) > 0L, logical(1))
    defs <- ids[haskid]; plhs <- ids[!haskid]
    if (length(defs) != 1L || length(plhs) < 1L) {
      stop("hybrid ", hl, " not in defining + placeholder form")
    }
    D <- defs[1]
    for (Pl in plhs) {
      ret_from <- c(ret_from, P[Pl]); ret_to <- c(ret_to, D)
      rlab <- c(rlab, hl); maj <- c(maj, gam[D]); minr <- c(minr, gam[Pl])
      remove[Pl] <- TRUE
    }
  }

  # Prune nodes that have no sampled (real-taxon) descendants. After merging
  # placeholders away, two such cases arise in nested/stacked reticulations:
  #   * "ghost" donors -- an unlabeled node whose children were all placeholders;
  #   * "leaf" hybrids -- a #Hk defining node whose subtree was all placeholders.
  # Neither is a real taxon, so remove it (iterated, so chains collapse): drop
  # any reticulation it *receives*, and re-point any it *donates* to its parent.
  # Real taxa (labelled, non-"#" tips) are never pruned, so the taxon set is
  # preserved.
  is_real <- !is.na(label) & nzchar(label) & !startsWith(ifelse(is.na(label), "", label), "#")
  repeat {
    rmv <- which(remove)
    childless <- vapply(seq_len(nn), function(i)
      length(setdiff(children[[i]], rmv)) == 0L, logical(1))
    prunable <- which(!remove & childless & P != 0L & !is_real)
    if (!length(prunable)) break
    for (ph in prunable) {
      drop <- ret_to == ph
      if (any(drop)) {
        ret_from <- ret_from[!drop]; ret_to <- ret_to[!drop]
        rlab <- rlab[!drop]; maj <- maj[!drop]; minr <- minr[!drop]
      }
      ret_from[ret_from == ph] <- P[ph]
      remove[ph] <- TRUE
    }
  }

  surv <- which(!remove)
  removed <- which(remove)
  surv_children <- lapply(surv, function(i) setdiff(children[[i]], removed))
  is_tip <- vapply(surv_children, length, integer(1)) == 0L
  tips <- surv[is_tip]
  root <- surv[P[surv] == 0L]
  if (length(root) != 1L) stop("could not identify a unique root")
  internals <- surv[!is_tip]
  ntip <- length(tips); nint <- length(internals)
  if (ntip < 1L) stop("no tips parsed")

  new_id <- integer(nn)
  new_id[tips] <- seq_len(ntip)
  new_id[root] <- ntip + 1L
  rest <- setdiff(internals, root)
  if (length(rest)) new_id[rest] <- (ntip + 2L):(ntip + 1L + length(rest))

  nonroot <- surv[surv != root]
  edge <- cbind(new_id[P[nonroot]], new_id[nonroot])
  storage.mode(edge) <- "integer"
  edge.length <- len[nonroot]

  tip.label <- label[tips]                  # already in new-id (1..ntip) order
  node.label <- rep("", nint)
  for (k in seq_along(rlab)) node.label[new_id[ret_to[k]] - ntip] <- rlab[k]

  ev <- list(edge = edge, Nnode = nint,
             tip.label = tip.label, node.label = node.label)
  if (!any(is.na(edge.length))) ev$edge.length <- edge.length
  ev$reticulation <- if (length(ret_from)) {
    r <- cbind(new_id[ret_from], new_id[ret_to]); storage.mode(r) <- "integer"; r
  } else matrix(integer(0), ncol = 2)
  class(ev) <- c("evonet", "phylo")

  gammas <- if (length(rlab)) {
    data.frame(label = rlab, tree_gamma = maj, ret_gamma = minr,
               stringsAsFactors = FALSE)
  } else {
    data.frame(label = character(), tree_gamma = numeric(),
               ret_gamma = numeric(), stringsAsFactors = FALSE)
  }
  gammas <- gammas[!duplicated(gammas$label), , drop = FALSE]

  hybrid_nodes <- if (length(rlab)) {
    hn <- stats::setNames(new_id[ret_to], rlab); hn[!duplicated(names(hn))]
  } else integer(0)

  issues <- character(0)
  if (any(is.na(edge.length))) {
    issues <- c(issues, "some branch lengths missing; treated as cladogram")
  }

  list(evonet = ev, gammas = gammas, hybrid_nodes = hybrid_nodes, issues = issues)
}
