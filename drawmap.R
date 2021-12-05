drawmap <-
  function(xy = spa, clusters, main = "Clusters in SORad", tcol = "black") {

    plot(
      xy,
      asp = 1,
      type = "n",
      main = main,
      xlab = "Longitude",
      ylab = "Latitude",
      bg = 'transparent'

    )
    # Add the clusters
    k <- length(levels(factor(clusters)))
    for (i in 1:k)
    {
      points(
        xy[clusters == i, 1],
        xy[clusters == i, 2],
        pch = 20,
        cex = 3,
        col = i + 1,
        bg = i + 1
      )
    }
    #text(xy,
    #     row.names(factors$Site),
    #     cex = 0.7,
    #     col = tcol,
    #     font = 2)
    legend(
      "topright",
      paste("Cluster", 1:k),
#      pch = (1:k) + 20,
      pch = 20,
      col = 2:(k + 1),
      pt.bg = 2:(k + 1),
      pt.cex = 2,
      bty = "n"
    )
    
  }
