server <- shinyServer(function(input, output, session) {
  
  #####
  #declare global variables
  ####
  nodes <- reactiveVal()
  edges <- reactiveVal()
  modules <- reactiveVal(c())
  
  translate_eids <-  function(teststr, altNames){
    eids <- str_split(teststr, pattern='\\/')[[1]]
    
    tids <- filter(altNames, entrezgene_id %in% eids) %>% pull(hgnc_symbol) %>% str_to_title() %>% paste(., collapse='/')
    
    return(tids)}
  
  
  output$modules <- renderUI({
    checkboxGroupInput(inputId='groupchoices',
               'Choose Modules to Graph',
               choices = modules)
  })
  
  #####
  #Analyze network
  ####
  observeEvent(input$run,{
    req(input$nodefile)
    req(input$edgefile)
    
    #import data
    nodes <<- read_delim(input$nodefile$datapath, delim = '\t')
    edges <<- read_delim(input$edgefile$datapath, delim = '\t')
    
    #format data
    nodes <<- nodes %>% rename('nodeName'='id') %>% mutate(label=id)
    edges <<- edges %>% rename('fromNode'='from', 'toNode'='to', 'weight'='value')
    nodes <<- nodes %>% select_if(~sum(!is.na(.)) > 0)
    edges <<- edges %>% select_if(~sum(!is.na(.)) > 0)
    #specific
    nodes <<- nodes %>% rename('nodeAttr[nodesPresent, ]'='group')
    modules <<- nodes %>% pull(group) %>% unique()
    
    
    #stats of node 
    #print(findFunction('select'))
    network <- graph.edgelist(as.matrix(edges %>% dplyr::select(c('to', 'from'))))
    
    #degree
    degs <- degree(network) %>% 
      data.frame('degree'=.) %>% rownames_to_column('id')
    #betweenness 
    betwn <- betweenness(network) %>% 
      data.frame('betweenness'=.) %>% rownames_to_column('id')
    
    #merge
    nodes <<- nodes %>% left_join(., degs, by ='id') %>% 
                        left_join(., betwn, by ='id') %>% 
                        dplyr::select(label, everything())
    
    #create checkbox choices from modules 
    output$modules <- renderUI({
      checkboxGroupInput(inputId='groupchoices',
                         'Choose Modules to Graph',
                         choices = modules)
    })
    
    #print data tables
    output$nodetable <- DT::renderDataTable({
      DT::datatable(nodes %>% dplyr::select(!label), list(pageLength = 5, scrollX=T))
    })
    output$edgetable <- DT::renderDataTable({
      DT::datatable(edges, list(pageLength = 5, scrollX=T))
    })
    
  })
  
  
  #####
  #Graph network
  ####
  observeEvent(input$graph,{
    
    
    #select / filter
    nnodes <<- nodes %>% filter(group %in% input$groupchoices) %>% mutate(color=group) %>% mutate(value=degree)
    modgenes <- nnodes$id
    eedges <<- edges %>% filter((from %in% modgenes) & (to %in% modgenes))

    #print(nodes)
    #print(edges)
    
    #if using only local degree
    if (input$whichdeg=='Local'){
      lnet <-  graph.edgelist(as.matrix(eedges %>% dplyr::select(c('to', 'from'))))
      ldeg <- degree(lnet) %>% data.frame('local_degree'=.) %>% rownames_to_column('id')
      nnodes <- left_join(nnodes, ldeg, by='id') %>% mutate(value=local_degree)
      nnodes <- nnodes %>% mutate(case_when(is.na(value) ~ 0, .default=nnodes$value))
    }
    
    output$graphnet <- renderVisNetwork({
      visNetwork(nnodes, eedges) %>%
        visIgraphLayout(physics=T) 
    })
  })
  
  #####
  # copy genes from selected modules
  ####
  output$copymods <- renderUI({
    req(input$groupchoices)
    #print('copy')
    #req(input$groupchoices)
    ngenes <- nodes %>% filter(group %in% input$groupchoices) %>% pull(label)
    cliptxt <- print(paste(ngenes, collapse = '\n'))
    rclipButton(inputId='clipbutton', 
                label='Copy module genes to clipboard',
                clipText=cliptxt,
                tooltip='So, you want to copy genes?')
  })

  
  
  
  
  #####
  #Fit powerlaw/Find hubs
  ####
  observeEvent(input$hubs,{
    #filter only if local analysis
    if (input$hublevel=='Local'){
      nnodes <- nodes %>% filter(group %in% input$groupchoices) %>% dplyr::select(!c(degree, betweenness))
      
      #init local network
      nnetwork <- graph.edgelist(as.matrix(edges %>% 
                                           filter(to %in% nnodes$id) %>% 
                                           filter(from %in% nnodes$id) %>% 
                                             dplyr::select(c('to', 'from'))))
      #degree
      degs <- degree(nnetwork) %>% 
        data.frame('degree'=.) %>% 
        rownames_to_column('id')
      degs <- degs %>% 
        mutate(degree=case_when(is.na(degree)~0, .default=degs$degree))
      #betweenness 
      betwn <- betweenness(nnetwork) %>% 
        data.frame('betweenness'=.) %>% 
        rownames_to_column('id')
      betwn <- betwn 
      
      #merge
      nnodes <- nnodes %>% left_join(., degs, by ='id') %>% 
        left_join(., betwn, by ='id') %>% 
        dplyr::select(label, everything()) 
      nnodes <- nnodes %>% 
        mutate(degree=case_when(is.na(degree)~0, .default=nnodes$degree))%>% 
        mutate(betweenness=case_when(is.na(betweenness)~0, .default=nnodes$betweenness))
    } else {nnodes <- nodes}
    
   
    
    #fit power law params
    #print(nnodes)
    pfit <- igraph::fit_power_law(nnodes$degree, start=5, implementation='plfit')
    alpha <- pfit$alpha
    xmin <- pfit$xmin
    
    #plot fit
    ppfit <- displ$new(nnodes %>% filter(degree>0) %>% pull(degree))
    ppfit$setXmin(pfit$xmin)
    ppfit$setPars(pfit$alpha)
    
    output$pfit <- renderPlot({plot(ppfit, xlab='Degree', ylab='Pr(Degree)') 
                              lines(ppfit, col='red')}, width=400, height=400)

        
    
    #find top 20%
    cdf <- dist_cdf(ppfit, lower_tail = T)
    cdf_cutoff <- 0.8
    gdat <- data.frame('Degree'=ppfit$dat, 'CDF'=cdf) %>% 
      mutate(sig=case_when(CDF>cdf_cutoff~T, .default=F))
    lowest_deg <- gdat %>% filter(sig==T) %>% pull(Degree) %>% min()
    
    #plot histogram with top 20%
    output$plhist<- renderPlot({
    ggplot(gdat %>% unique())+
      geom_bar(mapping=aes(x=Degree, y=CDF, fill=sig), stat='identity') +
      scale_fill_manual(values=c('grey','red'))+
      geom_abline(intercept=cdf_cutoff, slope=0, color='red')+
      guides(fill='none')+
      theme_bw()
      
    })
    
    #render table of selected nodes
    output$hubnodes <- DT::renderDataTable({
      DT::datatable((nnodes %>% filter(degree>=lowest_deg) %>% arrange(-degree)), list(pageLength = 10, scrollX=T),rownames = F)
    })
    #set global hubgenes variable
    hubgenes <<- nnodes %>% filter(degree>=lowest_deg) %>% pull(label)
    print(length(hubgenes))
  })
  
  #####
  # copy hub genes identified
  ####
  output$copyhubs <- renderUI({
    req(hubgenes)
    print(hubgenes)
    #req(input$groupchoices)
    cliptxt <- print(paste(hubgenes, collapse = '\n'))
    rclipButton(inputId='clipbutton2', 
                label='Copy hub genes to clipboard',
                clipText=cliptxt,
                tooltip='So, you want to copy genes?')
  })
  
  #####
  # Kegg analysis
  ####
  observeEvent(input$runkegg, {
    
    print('kegg')
    if (input$kegglevel=='Local'){genelist <- strsplit(input$kegg_glist, split='\\s')[[1]]
    } else { genelist <- nodes %>% pull(label) }
    print(genelist)
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    altNames <- tryCatch(getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene_id"),
                      values=genelist,mart= mart,uniqueRows = F),
                      error=function(e) e, warning=function(w){print(w)})
    print(altNames)
    
    results <- clusterProfiler::enrichKEGG(altNames$entrezgene_id, keyType = 'kegg', organism='hsa', pAdjustMethod = 'BH')
    results <- results@result
    results <- results %>% group_by(geneID) %>% mutate(geneID=translate_eids(geneID, altNames)) %>% ungroup()
    
    output$hubtable <- DT::renderDataTable({
      DT::datatable(results, list(pageLength = 10, scrollX=T))
    })
    Sys.sleep(2)
  })
})