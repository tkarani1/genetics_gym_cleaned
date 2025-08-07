source('constants.R')

score_colors_f2 = 
  c('popeve_neg_percentile' = '#007BC3',
    'eve_percentile' = '#54BCD1',
    'MisFit_S_percentile' = '#C70E7B',
    'MisFit_D_percentile' = '#FC6882',
    'score_PAI3D_percentile' = '#F4B95A',
    'esm1_v_neg_percentile' = '#8FDA04',
    'esm_score_neg_percentile' = '#009F3F',
    'am_pathogenicity_percentile' = '#EF7C12',
    'mpc_percentile' = '#B25D91',
    'proteinmpnn_llr_neg_percentile' = '#AF6125',
    'rasp_score_percentile' = '#F4E3C7'
  )


library(tidyverse)
library(scales)

load_grouped_data = function() {
  filenames = list.files('results', full.names = TRUE, pattern = 'with_std_errors\\.tsv$')
  grouped_filenames = filenames[(grepl('_grouped', filenames))]
  grouped_data = purrr::map_df(grouped_filenames, readr::read_tsv, id='evaluation') %>%
    mutate(cutoff = 1 - as.numeric(cutoff),
           evaluation=str_extract(evaluation, 'results/(.*)_results_with_std_errors', group=1)) %>%
    separate_wider_delim(evaluation, '_', names=c('evaluation', NA, 'category'),
                         too_few='align_start', too_many='merge') %>% 
    mutate(score_column=fct_relevel(score_column, names(score_colors_f2)),
           evaluation=fct_relevel(evaluation, names(evaluation_names)),
           category=if_else(is.na(category), 'all', category),
           ymin=enrichment - std_error,
           ymax=enrichment + std_error)
  return(grouped_data)
}
grouped_data = load_grouped_data()

load_contingency_table_data = function(fname) {
  gene_median = read_tsv(paste0('results/', fname)) %>%
    mutate(cutoff = 1 - as.numeric(cutoff), evaluation = 'gnomad') %>%
    mutate(score_column = gsub('_gene_median', '', score_column)) %>%
    select(-all_genes) %>%
    mutate(score_column=fct_relevel(score_column, names(score_colors_f2)),
           evaluation=fct_relevel(evaluation, names(evaluation_names)),
           high_score_neg = high_score_n - high_score_pos,
           total_neg = total_n - total_pos,
           cutoff_string=fct_relevel(percent(cutoff), c('1%', '2%', '5%', '10%')))  %>%
    mutate(result = pmap(list(high_score_pos, high_score_neg, total_pos, total_neg), f_test)) %>%
    unnest(result)
  return(gene_median)
}
ungrouped = load_contingency_table_data('gnomad_ungrouped_contingency_tables.tsv')
gene_median = load_contingency_table_data('gnomad_gene_median_contingency_tables.tsv')

plot_pair = function(evaluation_type = 'clinvar', x_str = 'dominant', y_str = 'recessive',
                     manual_score_filter=T, title=T, ylab='') {
  plot_data = grouped_data %>%
    filter(evaluation == evaluation_type & (category == x_str | category == y_str)) %>%
    select(-c(ymin, ymax))
  if (manual_score_filter) {
    plot_data = plot_data %>% filter(score_column %in% scores_to_plot)
  }
  if (ylab == '') {
    ylab = get_category_name(y_str)
  }
  p = plot_data %>%
    pivot_wider(names_from=category, values_from=c(enrichment, std_error)) %>%
    mutate(
      x := get(!!paste0("enrichment_", x_str)), y = get(!!paste0("enrichment_", y_str)),
      xlow := get(!!paste0("enrichment_", x_str)) - get(!!paste0("std_error_", x_str)),
           xhigh := get(!!paste0("enrichment_", x_str)) + get(!!paste0("std_error_", x_str)),
           ylow := get(!!paste0("enrichment_", y_str)) - get(!!paste0("std_error_", y_str)),
           yhigh := get(!!paste0("enrichment_", y_str)) + get(!!paste0("std_error_", y_str))) %>%
    ggplot + aes(color = score_column, x = x, y = y,
                 xmin = xlow, xmax = xhigh, ymin = ylow, ymax = yhigh) +
    geom_pointrange() + geom_errorbarh(height=0) +
    scale_color_manual(values=score_colors, labels=score_names, name='Method') +
    theme_classic() + 
    xlab(paste0('Enrichment (', get_category_name(x_str), ')')) +
    ylab(paste0('Enrichment (', ylab, ')')) +
    geom_abline(slope=1, intercept=0, linetype='dashed') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    xlim(c(1, NA)) + ylim(c(1, NA))
  if (is.character(title)) {
    p = p + labs(title=get_category_name(title))
  } else {
    if (title) {
      p = p + labs(title=evaluation_names[[evaluation_type]])
    }
  }
  return(p)
}
plot_pair()


plot_single_group = function(evaluation_type, category_name, title=T) {
  p = grouped_data %>%
    filter(evaluation == evaluation_type & category == category_name) %>%
    filter(score_column %in% scores_to_plot) %>%
    ggplot + aes(x = score_column, y = enrichment, color = score_column,
                 ymin = enrichment - std_error, ymax = enrichment + std_error) +
    geom_pointrange() +
    scale_color_manual(values=score_colors, labels=score_names, name='Method') +
    theme_classic() +
    scale_x_discrete(name=NULL, labels=score_names) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    geom_hline(yintercept=1, linetype='dashed') +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       limits=c(0.9, NA), breaks=1:20)
  if (title) {
    p = p + labs(title=evaluation_names[[evaluation_type]])
  }
  return(p)
}

plot_gene_median = function(manual_score_filter = T, plot_data=gene_median,
                            ypoint='estimate', yminpoint='conf.low', ymaxpoint='conf.high',
                            title='') {
  if (manual_score_filter) {
    plot_data = plot_data %>% filter(score_column %in% scores_to_plot)
    print(nrow(plot_data))
  }
  plot_data %>%
    filter(cutoff_string == '5%') %>%
    ggplot + aes(x = score_column, color = score_column) +
    aes_string(y = ypoint, ymin = yminpoint, ymax = ymaxpoint) +
    theme_classic() +
    guides(color = guide_legend(nrow = 1)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ylab('Enrichment') +
    geom_pointrange(size=0.5) -> p
  p = p + scale_color_manual(values=score_colors, labels=score_names, name='Method') +
    scale_x_discrete(name=NULL, labels=score_names)
  if (title != '') {
    p = p + labs(title=title)
  }
  return(p)
}


plot_paired = function(evaluation_type='ddd', manual_score_filter=T,
                                    points=T, cutoff_to_use=0.05, title=T,
                                    data=stderr_data, pair_plot=T) {
  enrichment_datasets = c('fine_mapped', 'genebass', 'schema', 'asc', 'epi25')
  plot_data = data %>%
    filter(abs(cutoff - cutoff_to_use)/cutoff < 0.01) %>%
    filter(evaluation == evaluation_type | evaluation == paste0(evaluation_type, '_gene_median'))
  if (manual_score_filter) {
    plot_data = plot_data %>% filter(score_column %in% scores_to_plot_no_percentile) %>%
      mutate(score_column = fct_drop(score_column))
    score_labels = score_names_no_percentile[scores_to_plot_no_percentile]
  } else {
    score_labels = score_names_no_percentile
  }
  plot_data = plot_data %>%
    mutate(
      median=grepl('median', evaluation),
      xpos=as.numeric(score_column) + if_else(median, 1/8, -1/8)) 
  # column_name = if_else(evaluation_type %in% enrichment_datasets, 'Enrichment', 'Rate ratio')
  # plot_data %>%
  #   ggplot + aes(x = xpos, y = result) + 
  #   aes(color = score_column, ymin = ylow, ymax = yhigh, shape=median) +
  #   geom_pointrange() +
  #   # geom_line(aes(group = score_column)) +
  #   scale_x_continuous(breaks = 1:length(score_labels), labels = score_labels, name=NULL) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  #   geom_hline(yintercept=1, linetype='dashed') +
  #   theme(legend.position='bottom') +
  #   scale_color_manual(values=score_colors_no_percentile, labels=score_names_no_percentile, name='Method') +
  #   scale_shape_manual(values=c('TRUE' = 15, 'FALSE' = 19), name=NULL,
  #                      labels=c('TRUE' = 'Gene score-based', 'FALSE' = 'Overall')) +
  #   guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) + ylab(column_name) +
  #   coord_cartesian(ylim=c(1, NA)) -> p1
  
  plot_data %>%
    group_by(score_column) %>%
    summarize(relative_result = sum(result*!median) - sum(result*median)) %>%
    mutate(median=F) %>%
    union_all(plot_data %>% filter(median) %>% select(score_column, relative_result=result, median)) %>%
    ggplot + aes(x = score_column, y = relative_result) + 
    aes(fill = median, shape=median) +
    geom_bar(stat='identity', position="stack") +
    # geom_line(aes(group = score_column)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    geom_hline(yintercept=1, linetype='dashed') +
    theme(legend.position='bottom') +
    scale_fill_manual(values=c('TRUE' = color_g, 'FALSE' = color_v2g), name=NULL,
                      labels=c('TRUE' = 'Gene score-based', 'FALSE' = 'Overall')) +
    scale_x_discrete(name=NULL, labels=score_names_no_percentile) +
    ylab(column_name) -> p1
  p1
  # plot_data %>%
  #   group_by(score_column) %>%
  #   summarize(relative_rate_ratio = sum(result*median)/sum(result*!median)) %>%
  #   ggplot + aes(x = score_column, y = relative_rate_ratio, fill = score_column) +
  #   geom_bar(stat='identity') +
  #   # geom_line(aes(group = score_column)) +
  #   scale_x_discrete(name=NULL, labels=score_names_no_percentile) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  #   geom_hline(yintercept=1, linetype='dashed') +
  #   theme(legend.position='bottom') +
  #   scale_fill_manual(values=score_colors_no_percentile, labels=score_names_no_percentile, name='Method') +
  #   # coord_cartesian(ylim=c(1, NA)) +
  #   guides(fill = guide_legend(nrow = 2)) + ylab('Relative rate ratio (V->G)') -> p2
  if (title) {
    p1 = p1 + labs(title=get_evaluation_name(evaluation_type))
  }
  # p = ggpubr::ggarrange(p1, p2, nrow=2, common.legend=T, legend = 'bottom')
  # p
  return(p1)
}


test = function() {
  ungrouped %>% 
    transmute(score_column, cutoff, cutoff_string, ungrouped=estimate) %>%
    left_join(gene_median %>% transmute(score_column, cutoff, gene_median=estimate)) %>%
    left_join(grouped_data %>%
                filter(evaluation == 'gnomad' & category == 'all') %>%
                transmute(score_column, cutoff, grouped=enrichment)) %>%
    filter(cutoff_string == "5%") -> test 
  test %>%
    ggplot + aes(x = ungrouped - gene_median, y = grouped) + geom_point()
  test %>%
    ggplot + aes(x = ungrouped, y = grouped) + geom_point()
  
  summary(lm(gene_median ~ ungrouped, data=test))
  summary(lm(grouped ~ ungrouped, data=test))
}

figure2 = function() {
  plots = map(names(evaluation_names), plot_paired)
  all_paired_plots = c(plots[1], plots[4], plots[7], plots[8])
  
  p1 = plot_gene_median(plot_data=ungrouped, title='Overall')
  p2 = plot_gene_median(title='Gene score-based')
  p3 = plot_gene_median(plot_data=grouped_data %>%
                          filter(evaluation == 'gnomad' & category == 'all') %>%
                          mutate(cutoff_string=fct_relevel(percent(cutoff), c('1%', '2%', '5%', '10%'))),
                          ypoint='enrichment', yminpoint='ymin', ymaxpoint='ymax',
                        title='Within-gene')
  all_plots = list(p1, p2, p3)
  top_p = ggpubr::ggarrange(plotlist=all_paired_plots, ncol=length(all_paired_plots), labels='auto',
                            common.legend=T, legend='bottom', align='h')
  bottom_p = ggpubr::ggarrange(plotlist=all_plots, ncol=length(all_plots), labels=c('e', 'f', 'g', 'h'),
                             common.legend=T, legend='bottom', align='h')
  png('figure2.png', res=res, width=2.5*length(all_paired_plots)*res, height=6.5*res, units = 'px')
  print(ggpubr::ggarrange(top_p, bottom_p, nrow=2))
  dev.off()
}
figure2()

all_pair_plots = function(eval_type='gnomad') {
  all_pair_plots = grouped_data %>% filter(category != 'all') %>% count(category) %>% pull(category) %>%
    map(function(x) plot_pair(eval_type, 'all', x, title=x, ylab='gene list'))
  png(paste0('extended_data_figure5_', eval_type, '.png'), res=res, width=9*res, height=9*res, units = 'px')
  print(ggpubr::ggarrange(plotlist=all_pair_plots, ncol=3, nrow=3, labels='auto',
                          common.legend=T, legend='bottom', align='h'))
  dev.off()
}
all_pair_plots()
all_pair_plots('clinvar')

pvals_data_gnomad = load_pairwise_data(T) %>% filter(!is.na(p_value))

extended_data_figure4 = function() {
  p1 = plot_pairwise_pvals('gnomad_ungrouped', data=pvals_data_gnomad, percentile=T, title='Overall')
  p2 = plot_pairwise_pvals('gnomad_gene_median', data=pvals_data_gnomad, percentile=T, title='Gene score-based')
  p3 = plot_pairwise_pvals('gnomad_grouped', data=pvals_data_gnomad, percentile=T, title='Within-gene')
  
  png('extended_data_figure4.png', res=res, width=180*edf_dpmm, height=90*edf_dpmm, units = 'px')
  print(ggpubr::ggarrange(plotlist=list(p1, p2, p3), ncol=3, labels='auto',
                          common.legend=T, legend='bottom', align='h'))
  dev.off()
}
extended_data_figure4()
