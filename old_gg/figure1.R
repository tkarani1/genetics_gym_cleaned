source('constants.R')

library(tidyverse)
library(scales)

load_stderr_data = function(suffix = '_results_with_std_errors') {
  filenames = list.files('results', full.names = TRUE, pattern = paste0(suffix, '\\.tsv$'))
  stderr_filenames = filenames[!grepl('_grouped', filenames)]
  stderr_data = purrr::map_df(stderr_filenames, readr::read_tsv, id='evaluation')
  if (!('rate_ratio' %in% names(stderr_data))) {
    stderr_data = stderr_data %>% mutate(rate_ratio=NA)
  }
  stderr_data = stderr_data %>%
    mutate(cutoff = 1 - as.numeric(cutoff),
           score_column = gsub('_gene_median', '', score_column),
           score_column = gsub('_percentile', '', score_column),
           evaluation=str_extract(evaluation, paste0('results/(.*)', suffix, '.tsv'), group=1)) %>%
    mutate(score_column=fct_relevel(score_column, names(score_colors_no_percentile)),
           evaluation=fct_relevel(evaluation, names(evaluation_names)),
           result=if_else(is.na(enrichment), rate_ratio, enrichment),
           ylow = result - std_error,
           yhigh = result + std_error,
           cutoff_string=fct_relevel(percent(cutoff), c('1%', '2%', '5%', '10%')))
  return(stderr_data)
}

stderr_data = load_stderr_data()
matched_stderr_data = load_stderr_data('_matched_eval_enr')
pvals_data = load_pairwise_data()

# plot_single_metric_all_cutoffs = function(evaluation_type = 'ddd', top=5,
#                                           manual_score_filter=T, legend=T,
#                                           error_bars=T, data=stderr_data,
#                                           max_point=25) {
#   column = if_else(evaluation_type %in% enrichment_datasets, 'enrichment', 'rate_ratio')
#   column_name = if_else(evaluation_type %in% enrichment_datasets, 'Enrichment', 'Rate ratio')
#   plot_data = data %>%
#     filter(evaluation == evaluation_type) %>%
#     mutate(
#       y := get(!!column),
#       ylow := get(!!column) - std_error,
#       yhigh := get(!!column) + std_error)
#   
#   plot_data = plot_data %>%
#     mutate(yhigh = ifelse(yhigh > max_point, max_point, yhigh))
#   
#   if (manual_score_filter) {
#     plot_data = plot_data %>% filter(score_column %in% scores_to_plot_no_percentile)
#   }
#   if (top > 0) {
#     plot_data %>%
#       filter(cutoff == min(cutoff)) %>%
#       slice_max(y, n=top) %>%
#       pull(score_column) -> metrics
#     plot_data = plot_data %>% filter(score_column %in% metrics)
#   }
#   
#   p = plot_data %>%
#     ggplot + aes(x = cutoff, y = rate_ratio, color = score_column, group = score_column)
#   
#   if (error_bars) {
#     p = p + aes(ymin = ylow, ymax = yhigh) + geom_pointrange(position=position_dodge(width=0.005))
#   } else {
#     p = p + geom_line()
#   }
#   p = p + theme_classic() +
#     scale_x_continuous(labels=percent, name='Top X%', breaks=plot_data %>% pull(cutoff) %>% unique) +
#     scale_color_manual(values=score_colors_no_percentile, labels=score_names_no_percentile, name='Method') +
#     ylab('Rate ratio') +
#     # scale_y_log10() +
#     # geom_segment(aes(xend = cutoff, yend = yhigh),
#     #              arrow = arrow(length = unit(0.2, "cm")),
#     #              data = plot_data %>% filter(yhigh == max_point)) +
#     coord_cartesian(ylim=c(NA, max_point)) 
#   if (legend) {
#     p = p + theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 1)) 
#   } else {
#     p = p + guides(color='none')
#   }
#   return(p)
# }

plot_with_stderrs = function(evaluation_type='ddd', manual_score_filter=T,
                             points=T, cutoff_to_use=0.05, title=T,
                             data=stderr_data) {
  enrichment_datasets = c('fine_mapped', 'genebass', 'schema', 'asc', 'epi25')
  column_name = if_else(evaluation_type %in% enrichment_datasets, 'Enrichment', 'Rate ratio')
  plot_data = data %>% filter(evaluation %in% names(evaluation_names))
  if (manual_score_filter) {
    plot_data = plot_data %>% filter(score_column %in% scores_to_plot_no_percentile)
  }
  if (cutoff_to_use > 0) {
    plot_data = plot_data %>% filter(abs(cutoff - cutoff_to_use)/cutoff < 0.01) %>%
      filter(evaluation == evaluation_type)
  }
  plot_data %>%
    ggplot + aes(y = result) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    geom_hline(yintercept=1, linetype='dashed') +
    theme(legend.position='bottom') +
    coord_cartesian(ylim=c(1, NA)) -> p
  if (points) {
    p = p + aes(color = score_column, ymin = ylow, ymax = yhigh) +
      geom_pointrange() +
      scale_color_manual(values=score_colors_no_percentile, labels=score_names_no_percentile, name='Method')
  } else {
    p = p + aes(fill = score_column) + geom_bar(stat='identity') +
      scale_fill_manual(values=score_colors_no_percentile, labels=score_names_no_percentile, name='Method')
  }
  if (cutoff_to_use == 0) {
    p = p + facet_grid(cols=vars(score_column), rows=vars(evaluation), scales='free',
                       labeller = labeller(evaluation = evaluation_names,
                                           score_column = score_names_no_percentile)) + 
      aes(x = cutoff_string) + xlab('Top X% of variants') +
      guides(color=F) + ylab('Metric (see caption)')
  } else {
    p = p + aes(x = score_column) + 
      scale_x_discrete(name=NULL, labels=score_names_no_percentile) +
      guides(color = guide_legend(nrow = 1)) + ylab(column_name)
    if (title) {
      p = p + labs(title=get_evaluation_name(evaluation_type))
    }
  }
  return(p)
}

figure1 = function(panels, data=stderr_data, width_mult=2.5) {
  all_plots = map(panels, plot_with_stderrs, data=data)

  res = 300
  png('figure1.png', res=res, width=width_mult*length(all_plots)*res, height=3.5*res, units = 'px')
  print(ggpubr::ggarrange(plotlist=all_plots, ncol=length(all_plots), labels='auto',
                          common.legend=T, legend='bottom', align='h'))
  dev.off()
}

figure1(c('ddd', 'schema', 'genebass', 'fine_mapped'))
# figure1(c('schema', 'genebass', 'fine_mapped'), data=matched_stderr_data)
figure1(c('schema', 'fine_mapped'), data=matched_stderr_data, width_mult=3.5)


extended_data_figure2 = function() {
  p1 = plot_with_stderrs(cutoff_to_use = 0, manual_score_filter = F, evaluation_type='genebass')
  
  png('extended_data_figure2.png', res=res, width=180*edf_dpmm, height=170*edf_dpmm, units = 'px')
  print(p1 + theme(strip.background = element_blank(),
                   strip.text.x = element_text(size = 5.5),
                   strip.text.y = element_text(size = 7),
                   axis.text.x = element_text(size = 7)))
  dev.off()
}
extended_data_figure2()


extended_data_figure3 = function() {
  all_pval_plots = pvals_data %>% count(evaluation) %>% pull(evaluation) %>%
    map(function(x) plot_pairwise_pvals(x, title=T))
  
  png('extended_data_figure3.png', res=res, width=180*edf_dpmm, height=170*edf_dpmm, units = 'px')
  print(ggpubr::ggarrange(plotlist=all_pval_plots, ncol=3, nrow=3, labels='auto',
                          common.legend=T, legend='bottom', align='h'))
  dev.off()
}
extended_data_figure3()


extended_data_figure1 = function() {
  roc = read_tsv('results/df_roc_ddd.tsv')
  p = roc %>%
    ggplot + aes(x = fpr, y = tpr, color = score_column) +
    annotate('rect', xmin = 0, xmax = 0.1, ymin = 0, ymax = 1, fill = 'gray90') +
    geom_line() + theme_classic() +
    scale_color_manual(values=score_colors_no_percentile, labels=score_names_no_percentile, name=NULL) +
    theme(legend.position = c(0.99, 0.01), legend.justification = c(1, 0)) +
    scale_x_continuous(labels=percent, 'False Positive Rate') +
    scale_y_continuous(labels=percent, 'True Positive Rate')
  png('extended_data_figure1.png', res=res, width=120*edf_dpmm, height=100*edf_dpmm, units = 'px')
  print(p)
  dev.off()
}
extended_data_figure1()
