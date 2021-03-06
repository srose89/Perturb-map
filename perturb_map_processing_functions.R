## 7/1/2021 ##
## PERTURB-map functions ##

## code used for the analysis of PERTURB-map data

# libraries required
require(sf)
require(ggvoronoi)
require(dplyr)
require(dbscan)


#### scd ####
#' function to use the single cell debarcoder strategy at a certain cutoff
#' to call barcodes for each cell
#' Input:
#' @pdat - procode epitope intensity data
#' @delta.cut - delta cutoff used for debarcoding
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
#' function that takes the processed image data and subclusters
#' each deconvoluted pro-code to remove outliers and separate independent tumors
#' input: 
#' @image.dat - must have a column named tags corresponding to the assigned procode
#' and x,y coordinates of segmented cells. A column with the linked Gene must also be
#' present
clusterPCmap <- function(image.dat, mp = 10, eps.value = 80, 
                         method = 'dbscan', 
                         dist.cutoff = 200, 
                         out.method = 'df'){
  
  # add an index to image dat so that I can rearrange at the end
  image.dat <- mutate(image.dat, 
                      idx = seq(from = 1, to = nrow(image.dat), by = 1))
  # get the unique tag combinations available in datasets
  tags <- unique(image.dat$tags[!(is.na(image.dat$tags))])

  # assign clusters for each cell within each procode
  pc.r <- do.call(dplyr::bind_rows, lapply(tags, function(tag){

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
#' function that produces a border around a set of points using the ahull method
#' need to iterate over alpha values to draw borders until a circle is formed
#' by connecting in igraph
#' input: 
#' @pcdat - data frame resulting from the dbscan cluster annotation
#' output: will return a list of sf polygon objects for each procode with a defined cluster
#' this will also include sizes for the polygons in a named list indexed by polygon or size
  
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

#### borderPolygon ####

# function to draw alpha shape around points and return polygon vector
#' input:
#' @pcdat - is segmented cell data with X, Y coordinates 
#' output is a polygon vector (sf object)

borderPolygon <- function(pcdat, alpha_start){
  # draw the alpha shape around this set of points
  pc.shape <- ashape(cbind(pcdat[,"Centroid X px"], pcdat[,"Centroid Y px"]), 
                     alpha = alpha_start)
  # assemble the points into connected lines
  pc.g <- graph.edgelist(cbind(as.character(pc.shape$edges[, "ind1"]), 
                               as.character(pc.shape$edges[, "ind2"])), directed = FALSE)
  
  # run tests to make sure alpha value used forms the correct shape of graph
  if (!is.connected(pc.g)) {
    stop("Graph not connected")
  }
  if (any(degree(pc.g) != 2)) {
    stop("Graph not circular")
  }
  if (clusters(pc.g)$no > 1) {
    stop("Graph composed of more than one circle")
  }
  
  # I will now just return the xy coords of the polygon for plotting in ggplot
  cutg = pc.g - E(pc.g)[1]
  # find chain end points
  ends = names(which(degree(cutg) == 1))
  path = unlist(get.shortest.paths(cutg, ends[1], ends[2])[[1]])
  # this is an index into the points
  pathX = as.numeric(V(pc.g)[path]$name)
  # join the ends
  pathX = c(pathX, pathX[1])

  
  # now for this I will return an sf polygon object because this will interface
  # with other spatial testing procedures
  # for plotting this can always be converted back to a matrix by using 
  # the st_coordinates function
  pc.poly <- st_polygon(list(pc.shape$x[pathX,]))
  return(pc.poly)
  
}

#### borderPoints ####
#' for a list of boundaries I need to get the overlapping points associated
#' input:
#' @pcdat - segmented cell data with xy coordinates
#' @borders - list of polygon objects
#' @border_scale - number if the border needs to be scaled by an affine transform
#' output: list of pcdat data frame filtered for points within each polygon
borderPoints <- function(pcdat, 
                         borders, 
                         border_scale = NULL){
  
  # first convert the segmented cell data into an sf object
  cells.sf <- st_as_sf(pcdat[,c('Centroid X px', 'Centroid Y px')],
                       coords = c('Centroid X px', 'Centroid Y px'))
  # if specified modify the border by a certain amount to get points
  # more in the center or more on the outside
  if(!(is.null(border_scale))){
    borders <- lapply(borders, function(b){
      b <- scalePoly(b, border_scale)
    })
  }
  # then for each border in the borders list extract the selected points in the 
  # pcdat data frame
  b.points <- lapply(borders, function(b){
    
    pcdat[st_intersects(cells.sf, b$polygon, sparse = F),]

    
  })
  names(b.points) <- names(borders)
  
  return(b.points)
  
}

#### sizeFilterPoly ####
#' function to filter out polygons that do not meet a certain area
#' input:
#' @p.l - list of polygons to filter
#' @area_threshold - size threshold for filtering out polygons
#' output:
#' size filtered polygon list

sizeFilterPoly <- function(p.l, area_threshold = 100){

  # create a return index to get the polygons that meet the area requirements
  ret_idx <- unlist(lapply(p.l, function(x){
    x$size > area_threshold
  }))
  # print out how many tumors are filtered out
  print(paste("kept", sum(ret_idx), "out of", length(p.l), "lesions"))
  # filter polygon list
  p.l.f <- p.l[ret_idx]
  # return filtered polygon list
  return(p.l.f)
}

#### scalePoly ####
#' function to scale a polygon by specific amount through an affine transformation
#' creating a setting where a portion can be subtracted from the transformation to
#' create a polygon with a hole in it
#' input:
#' @p - polygon to scale
#' @pscale - numeric constant to scale by
#' @subtract - if a subtraction of inner area is needed
#' @sub_scale - fraction of area within polygon to be subtracted
scalePoly <- function(p, pscale, 
                      subtract = FALSE, 
                      sub_scale = NULL){
  # extract polygon from list object
  p <- p$polygon
  # normal transform
  if(subtract == F){
    p <- (p - st_centroid(p)) * pscale + st_centroid(p)
    r <- list(polygon = p, 
              size = st_area(p))
    return(r)
  }
  # subtracting center region
  else{
    # check for value of sub_scale
    stopifnot(sub_scale > 0 & sub_scale < 1)
    p <- st_difference(((p - st_centroid(p)) * pscale + st_centroid(p)), 
                       ((p - st_centroid(p)) * sub_scale + st_centroid(p)))
    r <- list(polygon = p, 
              size = st_area(p))
    return(r)
  }
}

#### compositionQuant ####
## this will be a function to look at cellular composition
## based on different markers within tumor boundaries
#' input: 
#' @pcdat - x,y datapoints with columns of phenotype info (binary 1 or 0 for marker positivity)
#' @borders - boundary definitions (polygons)
#' @region - name of region being quantified

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