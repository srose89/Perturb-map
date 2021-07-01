## 7/1/2021 ##
## PERTURB-map functions ##

## code used for the analysis of PERTURB-map data

# libraries required
require(sf)
require(ggvoronoi)
require(dplyr)

#### scd ####
#' function to use the single cell debarcoder strategy at a certain cutoff
#' to call barcodes for each cell
#' Input:
#' @pdat - procode epitope intensity data
#' @delta.cut - delta cutoff used for debarcoding
#' @norm.method - whether to do percentile based or rescale
#' @out.method - whether to return a vector of procodes or a data frame with procodes and deltas
#' with no thresholding on cutoff
#' Output:
#' @pc - sorted procode tags of debarcoded cells separated by ';'

scd <- function(pdat, 
                delta.cut = .2, 
                out.method = 'vector')
{
  
  # normalize epitope tags
  
  pdat.norm <- rescale.norm(pdat)  
  
  # calculate deltas
  deltas <- calc.delta(pdat.norm)
  # change column names of pdat norm to just the tags themselves
  colnames(pdat.norm) <- gsub("Nucleus: | mean", "", colnames(pdat.norm))
  # for just the vector output
  if(out.method == 'vector'){
    # get barcodes for a certain delta threshold
    pc <- sapply(seq(from = 1, to = nrow(pdat.norm)), function(idx){
      if(deltas[idx] > delta.cut){
        # get channel ranks
        channel.rank <- rank(pdat.norm[idx,], ties.method = 'first')
        # take the top three channel names as the barcode
        bc <- paste(colnames(pdat.norm)[channel.rank > (ncol(pdat.norm) - 3)] %>% 
                      sort(), collapse = ";")
        # TO FOLLOW UP: can check if the barcode is part of library list
        
      }
      # if the delta is not above the cutoff then return NA
      else{
        bc <- NA
      }
      return(bc)
    })
    return(pc)
  }
  # for the data frame output
  if(out.method == 'df'){
    # get barcodes for all delta thresholds
    pc <- sapply(seq(from = 1, to = nrow(pdat.norm)), function(idx){
      # get channel ranks
      channel.rank <- rank(pdat.norm[idx,], ties.method = 'first')
      # take the top three channel names as the barcode
      paste(colnames(pdat.norm)[channel.rank > (ncol(pdat.norm) - 3)] %>% 
              sort(), collapse = ";")
    })
    pc.df <- data.frame(procode = pc, 
                        delta = deltas)
    # return the barcode for each cell along with a delta
    return(pc.df)
    
  }
}


#### rescale.norm ####
#' function to normalize the epitope tag intensity to percentile in the tissue
#' Input:
#' @pdat - epitope tag intensity in a matrix
#' Output:
#' - rescaled epitope tag intensities

rescale.norm <- function(pdat){
  pdat.norm <- apply(pdat, 2, function(x){
    scales::rescale(x, to = c(0,1))
  })
  return(pdat.norm)
}

#### calc.delta ####
#' this is a function to use the single cell debarcoder type functionality
#' in debarcoding perturb map data files
#' Input: 
#' @pdat matrix of normalized epitope tag intensity data from segmentation
#' Output:
#' deltas for each cell based on difference between third highest and
#' 4th highest channel
#' 
#' - 
calc.delta <- function(pdat.norm){
  
  # calculate deltas for everything
  deltas <- sapply(seq(from = 1, to = nrow(pdat.norm), by = 1), function(idx){
    # get channel ranks
    channel.rank <- rank(pdat.norm[idx,], ties.method = 'first')
    
    # return the delta
    pdat.norm[idx,channel.rank == (ncol(pdat.norm) - 2)] - pdat.norm[idx, channel.rank == (ncol(pdat.norm) - 3)]
  }) %>% as.numeric()
  
  return(deltas)
}

#### clusterPCmap ####
## function that takes the processed image data and subclusters
## each deconvoluted pro-code to remove outliers and separate independent tumors
clusterPCmap <- function(image.dat, mp = 10, eps.value = 80, 
                         method = 'dbscan', 
                         dist.cutoff = 200, 
                         out.method = 'df'){
  
  # add an index to image dat so that I can rearrange at the end
  image.dat <- mutate(image.dat, 
                      idx = seq(from = 1, to = nrow(image.dat), by = 1))
  # get the unique tag combinations available in datasets
  tags <- unique(image.dat$tags[!(is.na(image.dat$tags))])
  #print(tags)
  # assign clusters for each cell within each procode
  pc.r <- do.call(dplyr::bind_rows, lapply(tags, function(tag){
    #print(tag)
    image.pc <- dplyr::filter(image.dat, tags == tag)
    # for dbscan clustering
    # take x and y cooridinates of segmented cells
    clusters <- dbscan::dbscan(dplyr::select(image.pc, `Centroid X px`, `Centroid Y px`), 
                               minPts = mp , eps = eps.value)$cluster
    image.pc <- dplyr::bind_cols(image.pc, data.frame(cluster = clusters)) %>%
      tidyr::unite(tag_cluster, tags, cluster, remove = F)
    
    return(image.pc)
    
  }))
  
  ## need to put in a command to reassign clusters that have fewer than 10 points
  ## to the 0 cluster
  # how many cells per cluster?
  cpc <- table(pc.r$tag_cluster)
  rm_clusters <- names(cpc[cpc < 10])
  pc.r <- bind_rows(dplyr::filter(pc.r, !(tag_cluster %in% rm_clusters)), 
                    dplyr::filter(pc.r, tag_cluster %in% rm_clusters) %>%
                      mutate(cluster = 0) %>%
                      mutate(tag_cluster = paste(tags, cluster, sep = "_")))
  
  # bind with na cells
  pc.r <- dplyr::bind_rows(pc.r, dplyr::filter(image.dat, is.na(tags)))
  
  # add a gene cluster column for plotting and analysis
  pc.r <- mutate(pc.r, gene_cluster = paste(Gene, cluster, sep = "_"))
  
  # arrange in the original order and remove the index
  pc.r <- dplyr::arrange(pc.r, idx) %>% 
    dplyr::select(-idx)
  # output based on the method selected
  if(out.method == 'df'){
    return(pc.r)  
    break()
  }
  if(out.method == 'assignments'){
    return(pc.r[,c('tag_cluster', 'gene_cluster')])
    break()
  }
}
  
  
#### pcBorder ####
# function that produces a border around a set of points using the ahull method
#  need to iterate over alpha values to draw borders until a circle is formed
# by connecting in igraph
# input: data frame resulting from the dbscan cluster annotation
# value: will return a list of sp polygon objects for each procode with a defined cluster
## this will also include sizes for the polygons in a named list indexed by polygon or size
  
pcBorder <- function(pcdat, alpha_start = 80){
  
  # list of procodes to draw borders for
  pc.list <- unique(pcdat$tag_cluster)
  # filter out na values and outlier points defined by clustering
  pc.list <- pc.list[!(is.na(pc.list)) & !(grepl("_0", pc.list))]
  
  # loop over each procode to draw borders and get out sp polygon objects
  pc.pol <- lapply(pc.list, function(pc){
    
    # filter for the procode
    pcdat.f <- dplyr::filter(pcdat, tag_cluster == pc)
    
    # loop over function creating border and polygon around a set of points
    # until I have an alpha value that works
    repeat{
      print(paste("drawing borders w/ alpha", alpha_start))
      pol <- try(borderPolygon(pcdat.f, alpha_start))
      
      if(!(inherits(pol, "try-error")))
        break
      
      alpha_start <- alpha_start + 10
    }
    # get the size of the created polygon
    r <- list(polygon = pol, 
              size = st_area(pol))
    return(r)
  }
  )
  names(pc.pol) <- pc.list
  return(pc.pol)
  
}

#### borderPoints ####
# for a list of boundaries I need to get the overlapping points associated

borderPoints <- function(pcdat, 
                         borders, 
                         border_scale = NULL){
  
  # first convert the segmented cell data into an sf object
  cells.sf <- st_as_sf(pcdat[,c('Centroid X px', 'Centroid Y px')],
                       coords = c('Centroid X px', 'Centroid Y px'))
  # if specified modify the border by a certain amount to get points
  # more in the center or more on the outside
  if(!(is.null(border_scale))){
    #borders <- lapply(borders, st_buffer, border_scale, endCapStyple = "FLAT")
    borders <- lapply(borders, function(b){
      b <- scalePoly(b, border_scale)
    })
  }
  # then for each border in the borders list extract the selected points in the 
  # pcdat data frame
  b.points <- lapply(borders, function(b){
    
    pcdat[st_intersects(cells.sf, b$polygon, sparse = F),]
    #as.data.frame() %>%
    #dplyr::select(`Centroid X px` = X, `Centroid Y px` = Y) %>%
    #left_join(., pcdat, by = c('Centroid X px', 'Centroid Y px'))
    
  })
  names(b.points) <- names(borders)
  
  return(b.points)
  
}

#### sizeFilterPoly ####
# function to filter out polygons that do not meet a certain area


sizeFilterPoly <- function(p.l, area_threshold = 100){

  # create a return index to get the polygons that meet the area requirements
  ret_idx <- unlist(lapply(p.l, function(x){
    #print(st_area(x))
    x$size > area_threshold
  }))
  # print out how many tumors are filtered out
  print(paste("kept", sum(ret_idx), "out of", length(p.l), "lesions"))
  # filter polygon list
  p.l.f <- p.l[ret_idx]
  # return filtered polygon list
  return(p.l.f)
}

#### compositionQuant ####
## this will be a function to look at cellular composition
## based on different markers within tumor boundaries
# input: 
# 1 x,y datapoints with columns of phenotype info
# 2 boundary definitions
# thing to append onto the colnames to identify what regions they are

compositionQuant <- function(pcdat, borders, region) {
  # get cells within the boundaries
  pc.cells <- borderPoints(pcdat, borders)
  
  # within each region, what is the number of cells, frequency, and cells per area
  r <- do.call(rbind, 
               lapply(names(pc.cells), function(pc){
                 # get the index of columns to apply over, basically all that are not
                 # the centroid coordinates
                 marker_idx <- !(grepl("Centroid", colnames(pc.cells[[pc]])))
                 # get sums of positive cells for all markers
                 #s <- apply(pc.cells[[pc]][,3:ncol(pc.cells[[pc]])], 2, function(x){
                 s <- apply(pc.cells[[pc]][,marker_idx], 2, function(x){
                   sum(x)
                 })
                 # convert to list so I can bind it together later
                 s.list <- as.list(s)
                 
                 
                 # get frequencies of cells
                 f <- as.list(s / nrow(pc.cells[[pc]]))
                 names(f) <- paste(names(f), "freq")
                 
                 
                 # get cells per area
                 a <- as.list(s / borders[[pc]]$size)
                 names(a) <- paste(names(a), 'cells_per_area')
                 
                 # put all these columns together
                 p <- as.data.frame(cbind(tag_cluster = pc, 
                                          size = borders[[pc]]$size,
                                          as.data.frame(s.list), 
                                          as.data.frame(f), 
                                          as.data.frame(a) 
                 ))
                 return(p)
               }))
  colnames(r) <- c('tag_cluster', paste(colnames(r)[2:ncol(r)], region, sep = "."))
  return(r)
}