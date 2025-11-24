### _plot_fig.R
library(RColorBrewer)
library(ggdensity)
library(tidyr)

palette_box_0 <- c(
  "#5A8FB7",
  "#FFBA75",
  "#B99AD6",
  "#85D197"

)


palette_line2 <- c(
  "#444444",
  "#999999"

)

palette_line <- c(
 "#377EB8", "#E41A1C","#4DAF4A"
)


palette_box <- c(
  "#E41A1C","#377EB8","#4DAF4A","#984EA3")


# brewer.pal(8, 'Set2')
palette_set2<-c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3")

palette_dark2<-c(
  "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666")

brewer.pal(9, 'Set1')

brewer.pal(9, 'Pastel1')

palette_Pastel1<-c("#FBB4AE","#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F2F2F2")

palette_set1<-c(
"#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00"
,"#FFFF33","#A65628","#F781BF","#999999")



plot_dens2<-function(plot_data,plot_title,target_1,target_2){
  
  # Creates a 2D density contour plot with RdBu color scheme to visualize the relationship 
  # between two continuous variables. This function generates filled contour plots that 
  # show density distributions without axis constraints, making it suitable for data 
  # with varying ranges. The plot uses binned density estimation with 10 levels and 
  # applies a red-blue diverging color palette to highlight density gradients. 
  # The function automatically removes the legend for cleaner visualization.
  
  ggplot(plot_data,aes(x=.data[[target_1]],y= .data[[target_2]]))+
    labs(title=stringr::str_wrap(plot_title,width=40),size=1)+
    # xlim(c(-1,1))+ ylim(c(-1,1))+
    xlab(target_1)+
    ylab(target_2)+
    geom_density_2d_filled(aes(fill = after_stat(level)),bins=10) +
    scale_fill_brewer(palette = "RdBu",direction=-1)  +
    # scale_fill_viridis_d(
    #   name = "Density"
    # ) +
    #geom_density_2d_filled(alpha=0.8)+
    # scale_fill_distiller(name = "Count (%)",
    #                      labels = function(x) paste0(round(x, 1), "%"),palette = 'OrRd',direction=1)+
    # scale_fill_continuous(name = "Count (%)",
    #                      labels = function(x) paste0(round(x, 1), "%"),palette = 'OrRd',direction=1)+
    theme_minimal()+
    theme(legend.position = "none")
}


plot_dens<-function(plot_data,plot_title,target_1,target_2){
  
  # Generates a constrained 2D density contour plot with fixed axis limits from -1 to 1.
  # This specialized version is designed for normalized data within the [-1,1] range,
  # typically used for standardized structural metrics. The function creates filled 
  # contour plots with 10 density bins using the RdBu diverging color scheme, where
  # blue indicates low density and red indicates high density regions. Like plot_dens2,
  # it removes the legend for minimalistic presentation.
  
  ggplot(plot_data,aes(x=.data[[target_1]],y= .data[[target_2]]))+
    labs(title=stringr::str_wrap(plot_title,width=40),size=1)+
    xlim(c(-1,1))+ ylim(c(-1,1))+
    xlab(target_1)+
    ylab(target_2)+
    geom_density_2d_filled(aes(fill = after_stat(level)),bins=10) +
    scale_fill_brewer(palette = "RdBu",direction=-1)  +

    # scale_fill_viridis_d(
    #   name = "Density"
    # ) +
    #geom_density_2d_filled(alpha=0.8)+
    # scale_fill_distiller(name = "Count (%)",
    #                      labels = function(x) paste0(round(x, 1), "%"),palette = 'OrRd',direction=1)+
    # scale_fill_continuous(name = "Count (%)",
    #                      labels = function(x) paste0(round(x, 1), "%"),palette = 'OrRd',direction=1)+
    theme_minimal()+
    theme(legend.position = "none")
}

plot_dens_combn_cf<-function(plot_df,type_var,type_order,topo_used,pdf_file){

  # Produces a combined density plot panel focusing specifically on CF10QS relationships.
  # This function generates multiple 2D density plots comparing CF10QS against 
  # CF10QS.Trend4bi across different categories defined by type_var. It creates a 
  # horizontal panel layout where each category gets two adjacent plots showing the
  # same topological comparison. The output is saved as a PDF with dynamically 
  # calculated dimensions based on the number of categories.
  
  sec_topo<-c("CF10QS.Trend4bi")
  p_list<-list()
  i<-1
  
  for(plot_type in type_order){
    plot_data<-plot_df[which(plot_df[,type_var]==plot_type),]
    
    for(j in 1:length(sec_topo)){
      
      target_1<-"CF10QS" 
      target_2<-sec_topo[j]
      
      plot_title<-paste0(plot_type," (",nrow(plot_data),")")
      
      p_list[[i]]<-plot_dens(plot_data,plot_title,target_1,target_2)
      
      i<-i+1
      
    }
    
  }
  
  p<-ggarrange(plotlist=p_list,ncol=2*length(type_order),nrow=1)
  
  ggsave(pdf_file,p,width=2*2*length(type_order),height=2,limitsize = FALSE)
  
  
}


plot_dens_combn_full<-function(plot_df,type_var,type_order,topo_used,pdf_file){
  
  # Creates comprehensive density plot matrices comparing all combinations of topological 
  # features (LD15QS, CF10QS, and their trend derivatives). The function generates 
  # pairwise comparisons for all possible combinations of topological metrics across
  # multiple categories. It produces a grid layout where rows represent different 
  # categories and columns represent different topological feature pairs, providing
  # a complete overview of multivariate relationships in the dataset.
  
  topo_additional<-paste0(topo_used,c(".Trend4bi"))
  topo_all<-c("LD15QS","CF10QS",topo_additional)
  topo_pair<-combn(topo_all,2)
  
  p_list<-list()
  i<-1
  
  for(plot_type in type_order){
    plot_data<-plot_df[which(plot_df[,type_var]==plot_type),]
    
    for(j in 1:ncol(topo_pair)){
      
      target_1<-topo_pair[1,j]
      target_2<-topo_pair[2,j]
      
      plot_title<-paste0(plot_type," (",nrow(plot_data),")")
      
      p_list[[i]]<-plot_dens(plot_data,plot_title,target_1,target_2)
      
      i<-i+1
      
    }
    
  }
  
  p<-ggarrange(plotlist=p_list,ncol=ncol(topo_pair),nrow=length(type_order))
  
  ggsave(pdf_file,p,width=ncol(topo_pair)*2,height=length(type_order)*2,limitsize = FALSE)
  
  
}



plot_dens_combn<-function(plot_df,type_var,type_order,topo_used,npanels,pdf_prefix){
  # Generates modular density plot sets with flexible panel arrangements. This function
  # creates separate PDF files for each pairwise combination of topological features
  # (LD15QS, CF10QS, CF10QS.Trend4bi), organizing categories into configurable multi-panel
  # layouts. Unlike other density plot functions, it outputs multiple files with 
  # systematic naming based on the feature pairs, allowing for detailed examination
  # of individual relationships across different categorical groupings.
  
  
  # topo_additional<-paste0(topo_used,c(".Trend4bi"))
  topo_all<-c("LD15QS","CF10QS","CF10QS.Trend4bi")
  topo_pair<-combn(topo_all,2)
  
  for(j in 1:ncol(topo_pair)){
    
    
    target_1<-topo_pair[1,j]
    target_2<-topo_pair[2,j]
    
    p_list<-list()
    i<-1
    
    for(plot_type in type_order){
      
      plot_data<-plot_df[which(plot_df[,type_var]==plot_type),]

      
      plot_title<-paste0(plot_type," (",nrow(plot_data),")")
      
      p_list[[i]]<-plot_dens(plot_data,plot_title,target_1,target_2)
      
      i<-i+1
      
    }
    
    
    row_pannels<-ceiling(length(type_order)/npanels)
    # print(row_pannels)
    # p<-ggarrange(plotlist=p_list,ncol=npanels,nrow=row_pannels,common.legend = T)
    p<-ggarrange(plotlist=p_list,ncol=npanels,nrow=row_pannels)
    
    ggsave(paste0(pdf_prefix,"_",target_1,"_",target_2,".pdf"),p,width=npanels*2,height=2*row_pannels,limitsize = FALSE)
    
  }
  

  
  
}



plot_bar_interval_10<-function(df,xvar,yvar,title_label,color_set){

  # Creates sophisticated bar-interval plots that display the relationship between 
  # a continuous x-variable and a continuous y-variable through decile-based aggregation.
  # The function divides the x-variable into 10 equal quantiles, calculates median x-values
  # and mean y-values for each quantile, and displays them as connected points with 
  # error bars representing the interquartile range (25th-75th percentiles). This 
  # visualization effectively shows trends and variability across the x-variable range
  # while providing robust statistical summaries for each segment.
  
    
  df$x<-df[,xvar]
  df$y<-df[,yvar]
  
  df<-df[which(!is.na(df$x)&!is.na(df$y)),]
  
  quantile_labels <- paste0("(", seq(0, 90, by = 10), "%,", seq(10, 100, by = 10), "%]")
  
  df_plot <- df %>%
    mutate(quantile_group = ntile(x, 10)) %>%
    group_by(quantile_group) %>%
    summarise(
      x_median = median(x),  
      y_mean = mean(y),       
      y_q1 = quantile(y, 0.25), # 
      y_q3 = quantile(y, 0.75)  # 
    ) %>%
    ungroup() %>%
    mutate(
      group_label = factor(quantile_group, 
                           levels = 1:10, 
                           labels = quantile_labels)
    )

  plot_title<-paste0(title_label, " ",nrow(df)," sites")

  p<-ggplot(df_plot, aes(x = x_median)) +
    labs(title=plot_title)+
    # ylim(c(y_min,y_max))+
    geom_errorbar(
      aes(ymin = y_q1, ymax = y_q3),
      width = diff(range(df_plot$x_median)) * 0.05,  
      color= color_set[2],
      size = 1
    ) +
    geom_point(aes(y = y_mean), size = 4,color= color_set[1]) +
    geom_line(aes(y = y_mean), linewidth = 1.2,color= color_set[1]) +
    labs(x = xvar, 
         y = yvar) +
    scale_y_continuous(
      expand = expansion(mult = c(0.1, 0.1))
    ) +
    # scale_fill_brewer(palette = "Dark2") +
    scale_fill_manual(values = palette_box) +
    scale_color_manual(values= palette_box) +
    theme_minimal(base_size = 12) +
    theme(
      plot.margin = unit(c(20, 5, 20, 5), "points"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 10),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      plot.caption = element_text(color = "gray50", hjust = 0)
    )  
  p
}


plot_bar_even_10<-function(df,plot_x,plot_y,xmin,xmax){
  
  # Generates proportional distribution line plots across evenly spaced intervals.
  # This function divides the continuous x-variable into 10 equal intervals between
  # specified min and max values, then calculates the percentage distribution of
  # categorical y-variable groups within each interval. The result is a multi-line
  # plot where each line represents a category's proportional change across the
  # x-variable range, providing insights into how category prevalence shifts with
  # changing x-values. Includes automatic x-axis labeling with interval ranges.
  
  df<-df[which(!is.na(df[,plot_x])&!is.na(df[,plot_y])),]
  
  gap<-(xmax-xmin)/10
  
  breaks <-seq(xmin, xmax, length.out = 11) 
  
  quantile_labels <- paste0("(", seq(breaks[1], breaks[10], by = gap), ",", seq(breaks[2], breaks[11], by = gap), "]")

  df$quantile_group <- cut(
    df[,plot_x],
    breaks = breaks,
    include.lowest = TRUE,
    labels = quantile_labels
  )
  
  result <- df %>%
    group_by(.data[[plot_y]]) %>%
    mutate(total_per_label = n()) %>%  
    group_by(.data[[plot_y]], quantile_group) %>%
    summarise(
      count = n(),
      percentage = n() / first(total_per_label) * 100, 
      .groups = "drop"
    ) %>%
    select(-count) %>%
    tidyr::pivot_wider(
      names_from = quantile_group,
      values_from = percentage,
      values_fill = 0
    ) %>%
    ungroup()

  
  plot_data <- result %>%
    pivot_longer(cols = -plot_y, names_to = "quantile", values_to = "percentage")

  color_set<-palette_line
  p<-ggplot(plot_data, aes(x = quantile, fill = .data[[plot_y]])) +
    geom_point(aes(x = quantile,y = percentage, color = .data[[plot_y]]), size = 3) +
    geom_line(aes(x = quantile,y = percentage, group = .data[[plot_y]],color = .data[[plot_y]]), linewidth = 1.2) +
    #scale_color_manual(values=color_set)+ 
    labs(title = "",
         x = plot_x,
         y = paste0("Proportion (%)")) +
    theme_minimal() +
    #scale_fill_brewer(palette = "Dark2") +
    scale_fill_manual(values = color_set) +
    scale_color_manual(values= color_set) +
    scale_x_discrete(limits = quantile_labels) +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

  p
}


plot_bar_even_10_logy<-function(df,plot_x,plot_y,xmin,xmax){
  
  # Creates proportional distribution plots with logarithmic y-axis scaling. This
  # function extends plot_bar_even_10 by applying base-10 logarithmic transformation
  # to the y-axis, which is particularly useful for visualizing proportional data
  # with large dynamic ranges or when emphasizing relative differences at lower
  # percentage values. The log scaling helps reveal patterns that might be obscured
  # in linear-scale plots when dealing with highly skewed proportional distributions.
  
  df<-df[which(!is.na(df[,plot_x])&!is.na(df[,plot_y])),]
  
  gap<-(xmax-xmin)/10
  
  breaks <-seq(xmin, xmax, length.out = 11) 
  
  quantile_labels <- paste0("(", seq(breaks[1], breaks[10], by = gap), ",", seq(breaks[2], breaks[11], by = gap), "]")

  
  df$quantile_group <- cut(
    df[,plot_x],
    breaks = breaks,
    include.lowest = TRUE,
    labels = quantile_labels
  )
  
  result <- df %>%
    group_by(.data[[plot_y]]) %>%
    mutate(total_per_label = n()) %>%  
    group_by(.data[[plot_y]], quantile_group) %>%
    summarise(
      count = n(),
      percentage = n() / first(total_per_label) * 100, 
      .groups = "drop"
    ) %>%
    select(-count) %>%
    tidyr::pivot_wider(
      names_from = quantile_group,
      values_from = percentage,
      values_fill = 0
    ) %>%
    ungroup()
  
  # plot_title<-paste0(title_label, " ",nrow(df)," sites")
  
  plot_data <- result %>%
    pivot_longer(cols = -plot_y, names_to = "quantile", values_to = "percentage")
  
  #color_set<-brewer.pal(4, 'Dark2')
  color_set<-palette_line
  p<-ggplot(plot_data, aes(x = quantile, fill = .data[[plot_y]])) +
    geom_point(aes(x = quantile,y = percentage, color = .data[[plot_y]]), size = 3) +
    geom_line(aes(x = quantile,y = percentage, group = .data[[plot_y]],color = .data[[plot_y]]), linewidth = 1.2) +
    #scale_color_manual(values=color_set)+ 
    labs(title = "",
         x = plot_x,
         y = paste0("Proportion (%)")) +
    theme_minimal() +
    scale_fill_manual(values = color_set) +
    scale_color_manual(values= color_set) +
    scale_x_discrete(limits = quantile_labels) +
    scale_y_log10(breaks = seq(2, 20, by = 3)) +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
  
  p
}

qu_boxplot<-function(plot_df,Qu,QuVar,plotVar,title_label){
  
  # Produces comparative boxplots for extreme quantile groups with statistical testing.
  # This function identifies the top and bottom quantiles (specified by Qu parameter) 
  # of a quantitative variable (QuVar) and creates side-by-side boxplots comparing 
  # the distribution of another variable (plotVar) between these extreme groups.
  # Automatically performs Wilcoxon rank-sum tests and displays significance 
  # annotations between the compared groups. Ideal for identifying differential
  # characteristics between high and low extremes of a distribution.
  
  plot_df<-plot_df[which((!is.na(plot_df[,QuVar]))&(!is.na(plot_df[,plotVar]))),]
  
  low_cut<-quantile(plot_df[,QuVar],Qu)
  high_cut<-quantile(plot_df[,QuVar],1-Qu)
  
  low_label<-paste0(QuVar,"_low(",Qu,")")
  high_label<-paste0(QuVar,"_high(",Qu,")")
  
  plot_df$Group<-"-"
  plot_df$Group[which(plot_df[,QuVar]<=low_cut)]<-low_label
  plot_df$Group[which(plot_df[,QuVar]>=high_cut)]<-high_label
  
  print(table(plot_df$Group))
  
  plot_data<-plot_df[which(plot_df$Group!="-"),]
  
  my_comparisons <- list(c(low_label,high_label))
  
  plot_title<-paste0(title_label, " ",nrow(plot_df)," sites")

  # color_set<-brewer.pal(4, 'Dark2')
  # color_set<-c("#4A4E69", "#9A8C98", "#C9ADA7", "#F2E9E4")
  p<-ggplot(plot_data,aes(x=factor(Group),y= .data[[plotVar]]))+
    labs(title=stringr::str_wrap(plot_title,width=40),size=1)+
    xlab("")+
    ylim(c(min(plot_data[,plotVar])-0.1,max(plot_data[,plotVar])+0.3))+
    geom_boxplot(aes(color=Group),outlier.shape = NA,
                 width = 0.6,
                 alpha = 0.85,
                 lwd = 0.8, 
                 fatten = 1.2 )+
    stat_compare_means(aes(label="p.signif"),comparisons = my_comparisons,method="wilcox.test")+
    # stat_compare_means(aes(label=paste0("p=",..p.format..)),comparisons = my_comparisons,method="wilcox.test")+
    scale_color_manual(values=palette_box)+ 
    scale_fill_manual(values = palette_box) +
    #scale_color_manual(values=color_set)+ 
    #scale_color_manual(palette = "Dark2") +
    #theme_bw() +
    theme_minimal() +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 20))+
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size=7))
  
  p
  
  
  
}



qu_dens<-function(plot_df,Qu,QuVar,x_var,y_var){
  
  # Generates paired 2D density plots for extreme quantile comparison. This function
  # selects the top and bottom quantiles of a specified variable and creates 
  # side-by-side density contour plots showing the bivariate distribution of two
  # other variables within each extreme group. Enables visual comparison of how
  # the relationship between x_var and y_var differs between high and low extremes
  # of the quantile variable, providing insights into conditional distributions
  # at distribution tails.
  
  plot_df<-plot_df[which(!is.na(plot_df[,QuVar])),]
  
  low_cut<-quantile(plot_df[,QuVar],Qu)
  high_cut<-quantile(plot_df[,QuVar],1-Qu)
  
  low_label<-paste0(QuVar,"_low(",Qu,")")
  high_label<-paste0(QuVar,"_high(",Qu,")")
  
  plot_df$Group<-"-"
  plot_df$Group[which(plot_df[,QuVar]<=low_cut)]<-low_label
  plot_df$Group[which(plot_df[,QuVar]>=high_cut)]<-high_label
  
  print(table(plot_df$Group))
  
  plot_data<-plot_df[which(plot_df$Group!="-"),]
  
  p<-list()
  p[[1]]<-plot_dens2(plot_data[which(plot_data$Group==low_label),],low_label,x_var,y_var)
  p[[2]]<-plot_dens2(plot_data[which(plot_data$Group==high_label),],high_label,x_var,y_var)
  
  p_comb<-ggarrange(plotlist=p,nrow=1,ncol=2)
  
  
  p_comb
  
}

