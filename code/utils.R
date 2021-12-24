

# customized ggplot theme
theme_custom <- ggplot2::theme_bw() +
  ggplot2::theme(axis.text = ggplot2::element_text(size=24),
                 axis.title =ggplot2::element_text(size=24),
                 axis.line = ggplot2::element_blank(),
                 panel.border = ggplot2::element_rect(size=1.5),
                 axis.ticks = ggplot2::element_line(size=1.5),
                 plot.title = ggplot2::element_text(size = 28, hjust =0.5, face = "bold"),
                 legend.text = element_text(size=20), legend.title = element_text(size=20))

theme_full <- ggplot2::theme_bw() + ggplot2::theme(axis.text = ggplot2::element_text(size=15),
                                                   axis.title = ggplot2::element_text(size=16),
                                                   axis.line = ggplot2::element_blank(),
                                                   panel.border = ggplot2::element_rect(size=1.5),
                                                   axis.ticks = ggplot2::element_line(size=1.5),
                                                   plot.title = ggplot2::element_text(size = 16, hjust =0.5, face="bold"),
                                                   panel.grid.major = ggplot2::element_blank(),
                                                   panel.grid.minor = ggplot2::element_blank())

#arctan transformation for Kd values
arcTrans <- function(x,b=2, g = 1) {
  y <- (atan((x+b)*g) + pi/2)/pi
}

#Function to format floats
formatNum <- function(i, limit = 0.01, digits =1, format="e") {
  r <- sapply(i, function(n) {
    if (n < limit) {
      formatC(n, digits = digits, format = format)
    } else {
      format(n, digits = digits)
    }
  })
  return(r)
}

# Function to visualize collapsed target networks
plotTarGroups <- function(ProcessTargetResults, targetsOfInterest = NULL) {

  mapReducedTargets <- ProcessTargetResults$targetCluster

  TarGroups <- lapply(names(mapReducedTargets), function(sourceName) {
    if (length(mapReducedTargets[[sourceName]]) != 0) {
      tibble(source = sourceName, target = mapReducedTargets[[sourceName]])
    }
  }) %>% bind_rows()

  if (!is.null(targetsOfInterest)) {
    TarGroupsSubset <- TarGroups %>% filter_all(any_vars(. %in% targetsOfInterest))
  } else {
    TarGroupsSubset <- TarGroups
  }

  nSource <- length(unique(TarGroupsSubset$source))
  groups <- list(selectedTarget = c(seq(nSource)),
                 mappedTargets = seq(nSource+1, nSource + nrow(TarGroupsSubset)))

  qgraph::qgraph(TarGroupsSubset, groups = groups,
         directed = FALSE, color = ggsci::pal_locuszoom(alpha = 0.7)(2) ,
         vsize = 6,  layout = "spring", label.cex =1, shape="circle",
         label.scale = TRUE, label.color = "black", legend = FALSE)
}


#' Calculate row z-scores of a matrix

mscale <- function(x, center = TRUE, scale = TRUE, censor = NULL, useMad = FALSE){
  if (scale & center) {
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) (y-median(y,na.rm = T))/meanAD(y))
    } else {
      x.scaled <- apply(x, 1, function(y) (y-mean(y,na.rm=T))/sd(y,na.rm = T))
    }
  } else if (center & !scale) {
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) (y-median(y,na.rm=T)))
    } else {
      x.scaled <- apply(x, 1, function(y) (y-mean(y,na.rm=T)))
    }
  } else if (!center & scale) {
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) y/meanAD(y))
    } else {
      x.scaled <- apply(x, 1, function(y) y/sd(y,na.rm = T))
    }
  } else {
    x.scaled <- t(x)
  }

  if (!is.null(censor)) {
    if (length(censor) == 1) {
      x.scaled[x.scaled > censor] <- censor
      x.scaled[x.scaled < -censor] <- -censor
    } else {
      x.scaled[x.scaled > censor[2]] <- censor[2] #higher limit
      x.scaled[x.scaled < censor[1]] <- censor[1] #lower limit
    }
  }
  return(t(as.matrix(x.scaled)))
}


#function for boxplot, given a t-test result table, coef/frequency matrix and cell line annotation object
plotDiffBox <- function(pTab, coefMat, cellAnno, fdrCut = 0.05) {
  #filter
  pTab.sig <- filter(pTab, p.adj <= fdrCut)

  #process genetic background table
  geneBack <- cellAnno
  geneBack <- geneBack[colnames(coefMat),]

  pList <- lapply(seq(nrow(pTab.sig)), function(i) {
    geno <- pTab.sig[i,]$mutName
    target <- pTab.sig[i,]$targetName
    pval <- pTab.sig[i,]$p
    plotTab <- tibble(id = colnames(coefMat),
                      mut = geneBack[, geno],
                      val = coefMat[target,]) %>%
      filter(!is.na(mut))

    if (str_detect(geno, "trisomy12|gain|del")) {
      plotTab <- mutate(plotTab, mut = ifelse(mut %in% c("1",1), "present","absent")) %>%
        mutate(mut = factor(mut, levels = c("absent", "present")))
    } else {
      plotTab <- mutate(plotTab, mut = ifelse(mut %in% c("M","1",1), "Mutated","Unmutated")) %>%
        mutate(mut = factor(mut, levels = c("Unmutated","Mutated")))
    }

    #count cases
    numTab <- group_by(plotTab, mut) %>%
      summarise(n=length(id))

    plotTab <- left_join(plotTab, numTab, by = "mut") %>%
      mutate(mutNum = sprintf("%s\n(n=%s)", mut,n)) %>%
      arrange(mut) %>%
      mutate(mutNum = factor(mutNum, levels = unique(mutNum)))

    if (geno == "FLT3.ITD") {
      genoType = "mutation"
      geno <- "FLT3-ITD"
    } else if (str_detect(geno, "trisomy12|gain|del|IGHV.status")) {
      genoType = ""
      geno <- str_replace(geno, "[.]"," ")
    } else {
      genoType = "mutations"
    }

    titleText <- sprintf("%s ~ %s %s", target, geno, genoType)
    pval <- formatNum(pval, digits = 1, format="e")
    titleText <- bquote(atop(.(titleText), "("~italic("P")~"="~.(pval)~")"))

    ggplot(plotTab, aes(x = mutNum,y = val)) +
      stat_boxplot(geom = "errorbar", width = 0.3) +
      geom_boxplot(outlier.shape = NA, col="black", width=0.4) +
      geom_beeswarm(cex=2, size =2, aes(col = mutNum)) + theme_classic() +
      xlab("") + ylab("Protein dependence") + ggtitle(titleText) + xlab("") +
      scale_color_manual(values = c("#0072B5FF","#BC3C29FF")) +
      theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
            axis.title = element_text(size=18),
            axis.text = element_text(size=18),
            plot.title = element_text(size= 20, face = "bold", hjust = 0.5),
            legend.position = "none")
  })
  names(pList) <- paste0(pTab.sig$targetName,"_", pTab.sig$mutName)
  return(pList)

}

runCamera <- function(exprMat, design, gmtFile, id = NULL,
                      contrast = ncol(design),  method = "camera", pCut = 0.05, direction = "both",
                      ifFDR = FALSE, removePrefix = NULL, plotTitle = "", insideLegend = FALSE,
                      setToHighlight = c(), setMap = NULL) {
  scaleFUN <- function(x) sprintf("%.1f", x)

  #prepare indices
  if (is.null(id)) id <- rownames(exprMat)

  if (is.character(gmtFile)) {
    idx <- limma::ids2indices(loadGSC(gmtFile)$gsc, id)
  } else {
    idx <- limma::ids2indices(gmtFile,id)
  }

  #run camera for fry
  if (method == "camera") {
    res <- limma::camera(exprMat, idx, design, contrast)
  } else if (method == "fry") {
    res <- limma::fry(exprMat, idx, design, contrast)
  }

  #plot enrichment results as bar plot

  plotTab <- res %>% rownames_to_column("Name")

  if (!is.null(removePrefix)) plotTab <- mutate(plotTab, Name = str_remove(Name, removePrefix))
  if(!is.null(setMap)) {
    plotTab <- mutate(plotTab, newName = setMap[match(Name, setMap$setName),]$pathwayName) %>%
      mutate(Name = ifelse(is.na(newName), Name, newName))
  }

  plotTab <- plotTab %>%
    mutate(Direction= factor(Direction, levels =c("Down","Up"))) %>%
    arrange(desc(Direction),desc(PValue)) %>%
    mutate(Name = factor(Name, levels = Name))

  if (direction != "both") {
    plotTab <- filter(plotTab, Direction == direction)
  }

  if (ifFDR) {
    plotTab <- dplyr::filter(plotTab, FDR <= pCut)
  } else {
    plotTab <- dplyr::filter(plotTab, PValue <= pCut)
  }

  colAxis <- ifelse(plotTab$Name %in% setToHighlight, "red", "black")

  #for higlighting sets
  if (nrow(plotTab) == 0) {
    print("No sets passed the criteria")
    return(list(enrichTab = res, enrichPlot = NULL))
  } else {
    p <- ggplot(data = plotTab, aes(x = Name, y = -log10(PValue), fill = Direction)) +
      geom_bar(position = "dodge", stat = "identity", width = 0.5) +
      scale_fill_manual(values = c(Up =  "#BC3C29FF", Down = "#0072B5FF")) +
      coord_flip() + xlab("") +
      scale_y_continuous(labels=scaleFUN) +
      ylab(expression(-log[10]*'('*italic(P)~value*')')) + ggtitle(plotTitle) +
      theme_full +
      theme(axis.text.y = element_text(size= 12, color = colAxis),
            axis.text.x = element_text(size=12),
            plot.title = element_text(size=14, face = "bold"))
    if (insideLegend) {
      p <- p + theme(legend.position = c(0.8,0.2))
    } else {
      p <- p + theme(legend.position = "right")
    }

    return(list(enrichTab = res, enrichPlot = p))
  }
}


# Function to perform a t-test or ANOVA test given a protein dependence matrix and a metadata object
diffImportance <- function(coefMat, Annotation) {
  #process genetic background table
  geneBack <- Annotation
  geneBack <- geneBack[colnames(coefMat),]
  keepCols <- apply(geneBack,2, function(x) length(unique(na.omit(x))) >=2 & all(table(x)>6))
  geneBack <- geneBack[,keepCols]

  pTab <- lapply(rownames(coefMat), function(targetName) {
    lapply(colnames(geneBack), function(mutName) {
      impVec <- coefMat[targetName, ]
      genoVec <- geneBack[, mutName]
      resTab <- data.frame(targetName = targetName, mutName = mutName,
                           stringsAsFactors = FALSE)
      if (length(unique(na.omit(genoVec))) == 2) {
        #binary feature, usting t.test
        res <- t.test(impVec ~ genoVec, var.equal = TRUE, na.action = na.omit)
        resTab$p <- res$p.value
        resTab$FC <- (res$estimate[[2]]-res$estimate[[1]])/abs(res$estimate[[1]])
      } else if (length(unique(na.omit(genoVec))) >=3) {
        #using anova
        res <- anova(lm(impVec ~ genoVec, na.action = na.omit))
        #get the group mean difference
        diffTab <- data.frame(val = impVec, gr = genoVec) %>%
          dplyr::filter(!is.na(gr)) %>% dplyr::group_by(gr) %>%
          dplyr::summarise(meanVal = mean(val))
        resTab$p = res$`Pr(>F)`[1]
        resTab$FC = max(diffTab$meanVal)-min(diffTab$meanVal)
      }
      resTab
    }) %>% dplyr::bind_rows()
  }) %>% dplyr::bind_rows() %>% dplyr::arrange(p) %>% dplyr::mutate(p.adj = p.adjust(p, method = "BH"))
  pTab
}
