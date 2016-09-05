library(shiny)
library(ggplot2)
library(RANN)
library(Matrix)
library(RSpectra)

# logifySlider javascript function
JS.logify <-
  "
// function to logify a sliderInput
function logifySlider (sliderId, sci = false) {
if (sci) {
// scientific style
$('#'+sliderId).data('ionRangeSlider').update({
'prettify': function (num) { return ('10<sup>'+num+'</sup>'); }
})
} else {
// regular number style
$('#'+sliderId).data('ionRangeSlider').update({
'prettify': function (num) { return (Math.pow(10, num)); }
})
}
}"

# call logifySlider for each relevant sliderInput
JS.onload <-
  "
// execute upon document loading
$(document).ready(function() {
// wait a few ms to allow other scripts to execute
setTimeout(function() {
// include call for each slider
logifySlider('log_slider', sci = false)
logifySlider('log_slider2', sci = true)
}, 5)})
"

ui <- fluidPage(
  tags$head(tags$script(HTML(JS.logify))),
  tags$head(tags$script(HTML(JS.onload))),
  
  # Application title
  #titlePanel("Shiny Diffusion Maps"),
  
  # Sidebar with slider inputs
  sidebarLayout(
    sidebarPanel(
      #this part of the side panel is for data and parameter selection
      fluidRow(selectInput("data", "Choose data set:", 
                           c("swiss roll","bickley jet (less turbulence)" = "tame bickley", 
                             "bickley jet (more turbulence)" = "wild bickley",
                             "double gyre 60x30", "double gyre 200x100")),
               sliderInput("log_slider2", "Choose value of epsilon (log scale):",
                           min = -3, max = 2, value = -2.0, step = 0.1),
               textOutput("readout2"),
               br(),
               sliderInput("N_k", "Choose number of clusters:",
                           min = 2, max = 20, value = 9)
               ),
      #hr(),
      # this part of the side panel controls what is displayed in plots 1 and 2
      fluidRow(selectInput("toggle", "Plot 1 display:", c("eigenvectors", "eigenvalues", 
                                                            "clustering", "embedding")),
               selectInput("toggle2", "Plot 2 display:", c("eigenvectors", "eigenvalues", 
                                                            "clustering", "embedding"),
                           selected = "embedding"),
               
               sliderInput("k", "Choose eigenfunction to display:", 
                           min = 1, max = 10, value = c(2,3)),
               uiOutput("Tslider"),
               uiOutput("TintervalSlider")
               )
    ),
    
    # Main Panel with Plots 1 and 2
    mainPanel(
      h4("Plot 1"),
      plotOutput("Plot1", height="260px"),
      h4("Plot 2"),
      plotOutput("Plot2", height="260px")
    )
  )
)

server <- function(input, output, session) {
  # set a starting flag, this makes sure that plots are only generated when all
  # reactive conductors are flushed
  values <- reactiveValues(starting = TRUE);
  session$onFlushed(function() {
    values$starting <- FALSE
  })
  # show the actual value of epsilon
  output$readout2 <- reactive({
    paste0("epsilon = ", round(10^input$log_slider2,4))
  })
  # make the time frame slider (with min and max values depending on dataset)
  output$Tslider <- renderUI({
    sliderInput("t", "Choose time slice:", min = 1, max = data()$T, 
                value = floor(data()$T / 2))
  })
  # make the time interval slider (with min and max values depending on dataset)
  output$TintervalSlider <- renderUI({
    sliderInput("t_interval", "Choose time interval:", min = 1, max = data()$T, 
                value = c(1, data()$T))
    # values$starting <- FALSE;
  })
  # retrieve data
  data <- reactive({
    if (input$data == "wild bickley") {
      dataX = read.csv('data/bickley_wild_x.csv', sep = ",", header = FALSE)
      dataY = read.csv('data/bickley_wild_y.csv', sep = ",", header = FALSE)
      r <- 0.2;
    } else if (input$data == "tame bickley") {
      dataX = read.csv('data/bickley_x.csv', sep = ",", header = FALSE)
      dataY = read.csv('data/bickley_y.csv', sep = ",", header = FALSE)
      r <- 0.2;
    }
    else if (input$data == "swiss roll") {
      dataX = read.csv('data/swissroll_x.csv', sep = ",", header = FALSE)
      dataY = read.csv('data/swissroll_y.csv', sep = ",", header = FALSE)
      r <- 0.5;
    } else if (input$data == "double gyre 60x30"){
      dataX = read.csv('data/gyre_x.csv', sep = ",", header = FALSE)
      dataY = read.csv('data/gyre_y.csv', sep = ",", header = FALSE)
      r <- 0.1;
    } else {
      dataX = read.csv('data/gyre_x_fine.csv', sep = ",", header = FALSE)
      dataY = read.csv('data/gyre_y_fine.csv', sep = ",", header = FALSE)
      r <- 0.01;
    }
    return(list(dataX = dataX, dataY = dataY, m = dim(dataX)[1], T = dim(dataX)[2], 
                radius = r))
  })
  # compute list of sparse distance matrices (stored as a list (x, y, v) where x[[j]],
  # y[[j]] and v[[j]] are the vectors that form the jth sparse matrix)
  distances <- reactive({
    x_list = list();
    y_list = list();
    v_list = list();
    T <- data()$T;
    m <- data()$m;
    for (j in 1:T) {
      df <- data.frame(data()$dataX[,j],data()$dataY[,j])
      nearest <- nn2(df, k = 20, searchtype = "radius", radius = data()$radius)
      lv <- sum(nearest$nn.idx > 0)
      x <- rep(0, lv)
      y <- rep(0, lv)
      v <- rep(0, lv)
      icurr <- 1;
      for (i in 1:m){
        idx_there <- (nearest$nn.idx[i,] > 0);
        li <- sum(idx_there);
        x[icurr:(icurr+li-1)] <- i;
        y[icurr:(icurr+li-1)] <- nearest$nn.idx[i,idx_there];
        v[icurr:(icurr+li-1)] <- nearest$nn.dists[i,idx_there];
        icurr = icurr + li;
      }
      x_list[[j]] <- x;
      y_list[[j]] <- y;
      v_list[[j]] <- v;
    }
    return(list(x = x_list, y = y_list, v = v_list))
  })
  # compute diffusion map eigenfunctions
  eigen <- reactive({
    P <- Matrix(0, nrow = data()$m, ncol = data()$m, sparse = TRUE)
    #get epsilon from input slider
    epsilon <- 10^input$log_slider2
    #get window of time slices used for computation from slider t_interval
    Tmin <- as.integer(input$t_interval[1]);
    Tmax <- as.integer(input$t_interval[2]);
    for (j in Tmin:Tmax){
      K <- sparseMatrix(i = distances()$x[[j]], j = distances()$y[[j]], 
                        x = exp(-distances()$v[[j]]^2/epsilon),
                        symmetric = FALSE);
      K <- (1/rowSums(K))*K;
      #K <- t(K) %*% K;
      #K <- (1/rowSums(K))*K;
      P <- P + K;
    }
    P <- 1/(Tmax - Tmin + 1)*P;
    #compute eigenvalues with ARPACK
    return(eigs(P, 20, which = "LM"))
  })
  # retrieves the eigenvector for plotting
  dom_evec <- reactive({
    evecs = Re(eigen()$vectors);
    dom_evec <- evecs[,as.integer(input$k[1])]
    return(sign(max(dom_evec) + min(dom_evec))*dom_evec)
  })
  # k means clustering based on spacetime diffusion coordinates
  cluster_idx <- reactive({
    evecs = Re(eigen()$vectors);
    kcluster <- kmeans(evecs[,2:as.integer(input$N_k)], as.integer(input$N_k), iter.max = 50, nstart = 100);
    return(factor(kcluster$cluster))
  })
  # get time frame for plotting
  plot_frame <- reactive({
    t_disp <- min(data()$T,as.integer(input$t));
    data.frame(x = data()$dataX[,t_disp], y = data()$dataY[,t_disp])
    });
  # generate Plot 1
  output$Plot1 <- renderPlot({
    # delay plotting until everything is flushed
    if (values$starting)
      return(NULL)
    if (input$toggle == "eigenvectors"){
      p1 <- ggplot(plot_frame(), aes(x,y))
      p1 + geom_point(aes(color=dom_evec())) + 
        scale_colour_gradientn(colours=rainbow(5), name = (paste("evec", input$k[1])));
    } else if (input$toggle == "clustering"){
      p1 <- ggplot(plot_frame(), aes(x,y))
      p1 + geom_point(aes(color = cluster_idx())) + 
        scale_color_brewer(palette="Spectral", name = "cluster ID");
    } else if (input$toggle == "eigenvalues"){
      evals <- data.frame(n = as.numeric(1:20), lambda_n = Re(eigen()$values));
      ggplot(evals, aes(n,lambda_n)) + geom_point(size = 5, color = "navy");
    } else {
      evecs = data.frame(x = Re(eigen()$vectors[,input$k[1]]),
                         y = Re(eigen()$vectors[,input$k[2]]));
      evecs$x <- sign(max(evecs$x) + min(evecs$x))*evecs$x;
      evecs$y <- sign(max(evecs$y) + min(evecs$y))*evecs$y;
      p1 <- ggplot(evecs, aes(x,y)) 
      p1 + geom_point(aes(color = cluster_idx())) + 
        scale_color_brewer(palette="Spectral", name = "Cluster ID") + 
        xlab(paste("evec", input$k[1])) + ylab(paste("evec", input$k[2]));
    }
  })
  # generate plot 2 (this is the same code that generated plot 1)                                                                       
  output$Plot2 <- renderPlot({
    # delay plotting until everything is flushed
    if (values$starting)
      return(NULL)
    if (input$toggle2 == "eigenvectors"){
      p2 <- ggplot(plot_frame(), aes(x,y))
      p2 + geom_point(aes(color=dom_evec())) + 
        scale_colour_gradientn(colours=rainbow(5), name = (paste("evec", input$k)));
    } else if (input$toggle2 == "clustering"){
      p2 <- ggplot(plot_frame(), aes(x,y))
      p2 + geom_point(aes(color=cluster_idx())) + 
        scale_color_brewer(palette="Spectral", name = "Cluster ID");
    } else if (input$toggle2 == "eigenvalues"){
      evals <- data.frame(n = as.numeric(1:20), lambda_n = Re(eigen()$values));
      ggplot(evals, aes(n,lambda_n)) + geom_point(size = 5, color = "navy");
    } else {
      evecs = data.frame(x = Re(eigen()$vectors[,input$k[1]]),
                         y = Re(eigen()$vectors[,input$k[2]]));
      evecs$x <- sign(max(evecs$x) + min(evecs$x))*evecs$x;
      evecs$y <- sign(max(evecs$y) + min(evecs$y))*evecs$y;
      p1 <- ggplot(evecs, aes(x,y)) 
      p1 + geom_point(aes(color = cluster_idx())) + 
        scale_color_brewer(palette="Spectral", name = "cluster ID") + 
        xlab(paste("evec", input$k[1])) + ylab(paste("evec", input$k[2]));
    }
  })
  session$onSessionEnded(function() {
    stopApp()
  })
}

shinyApp(ui, server)
