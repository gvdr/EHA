#' input/output functions

#' `adjacency_from_paj()` return a weighted adjacency matrix from
#' a Pajek project file
#' it gets `name` of the file (with the explicit extension)
#' and the `base` path (ending with a /) as inputs
adjacency_from_paj <- function(Pajek_project){
  library(network)
  library(igraph)
  library(intergraph)
  
#' We invoke explicitly all the functions environemnt,
#' so to avoid any conflict between the three graph libraries
  Net_paj <- network::read.paj(Pajek_project)
  Net_net <- Net_paj$networks[[1]]
  Web <- intergraph::asIgraph(Net_net)
  igraph::V(Web)$name <- igraph::V(Web)$vertex.names
  weight_name <- igraph::list.edge.attributes(Web)[
    igraph::list.edge.attributes(Web) != "na"
  ]
  weights <- as.vector(igraph::get.edge.attribute(Web,weight_name))
  Web <- igraph::set_edge_attr(Web,"weight", value = weights)
  if("na" %in% igraph::list.edge.attributes(Web)){
    Web <- igraph::remove.edge.attribute(Web,"na")
  }
  Web <- igraph::remove.edge.attribute(Web,weight_name)
  if("na" %in% igraph::list.vertex.attributes(Web)){
    Web <- igraph::remove.vertex.attribute(Web,"na")
  }
  Web <- igraph::remove.vertex.attribute(Web,"vertex.names")
  
  Adjacency_test <- igraph::get.adjacency(Web,attr = "weight",sparse = F)

  return(Adjacency_test)
}


expose_uniqueness <- function(webs_df,Web_name,Transformed,Mode_uniqueness,
                              Arrange = F){
  webs_df %>%
    filter(name_web %in% Web_name) %>%
    select_(Transformed) %>%
    unnest() %>%
    filter(Mode %in% Mode_uniqueness) %>%
    select(-Mode) -> Uniqs
  if(Arrange){
    Uniqs  %>%
      arrange(-MeanDist) %>%
      return()
  } else {
    Uniqs  %>% return()
  }
}

expose_uniqueness_all <- function(webs_df,Transformed,Mode_uniqueness,
                              Arrange = F){
  webs_df %>%
    select_("name_web",Transformed) %>%
    unnest() %>%
    filter(Mode %in% Mode_uniqueness) %>%
    select(-Mode) -> Uniqs
  if(Arrange){
    Uniqs  %>%
      arrange(-MeanDist) %>%
      return()
  } else {
    Uniqs  %>% return()
  }
}

expose_degree <- function(webs_df,Web_name,Degree_measurement,
                          Arrange = F){
  
  data_frame(node_name = webs_df %>%
               filter(name_web %in% Web_name) %>%
               select(adjacency) %>%
               unnest() %>% 
               names()
  ) %>% bind_cols(
    webs_df %>%
      filter(name_web %in% Web_name) %>%
      select_(Degree_measurement) %>%
      unnest()
  ) -> deg
  
  if(Arrange){
    deg  %>%
      arrange_(-Degree_measurement) %>%
      return()
  } else {
    deg  %>% return()
  }
}

expose_degree_all <- function(webs_df,Degree_measurement,
                          Arrange = F){
  
  web_names <- webs_df %$%
    name_web %>%
    unique()
  
  name_nodes <- character(0)
  for(name in web_names){
    webs_df %>%
      filter(name_web %in% name) %>%
      select(adjacency) %>%
      unnest() %>%
      names() %>%
      c(name_nodes,.) -> name_nodes
  }
  
  data_frame(name_web = webs_df %>%
               select(name_web,adjacency) %>%
               unnest() %$% 
               name_web,
             node_name = name_nodes) %>%
  bind_cols(
    webs_df %>%
      select_(Degree_measurement) %>%
      unnest()
  ) -> deg
  
  if(Arrange){
    deg  %>%
      arrange_(-Degree_measurement) %>%
      return()
  } else {
    deg  %>% return()
  }
}