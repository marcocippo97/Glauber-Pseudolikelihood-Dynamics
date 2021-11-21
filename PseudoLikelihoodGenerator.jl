module PseudoLikelihoodGenerator

import Statistics: mean, std, cor

import PCA: plot_nat_test_seqs, plot_nat_test_seqs_2, perform_pca, get_projections, plot_in_pc_space  

import DynamicPLM: loadparams, saveparams, plmdcaparam, compute_roc, PlmParam, ReadFasta

export loadparams, saveparams, plmdcaparam, compute_roc, PlmParam, ReadFasta

export plot_nat_test_seqs, plot_nat_test_seqs_2

export mean, std, cor

export dynamic, log_pseudo, dynamic_term, connected_correlation, freq_res, freq_res_2, autocorrelazione, autocorrelazione_dist, blockmean, num_contact

export contrastive, bipartite, mean_dist, get_first, dist, freq_sing_1, freq_sing_2, freq_sing_3, connected_sing_3, test_corr, connected_3, mean_change

include("plmgenerator.jl") 

end 