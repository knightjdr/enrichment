#!/usr/local/bin/Rscript

proteomeFraction = scan('tier1/background-proteome/wilcox_all_fraction.txt')
proteomeLongset = scan('tier1/background-proteome/wilcox_all_longest.txt')
rbpFraction = scan('tier1/background-RBPs/wilcox_all_fraction.txt')
rbpLongest = scan('tier1/background-RBPs/wilcox_all_longest.txt')
tier1Fraction = scan('tier1/background-proteome/wilcox_query_fraction.txt')
tier1Longset = scan('tier1/background-proteome/wilcox_query_longest.txt')
tier1_2Fraction = scan('tier1and2/background-proteome/wilcox_query_fraction.txt')
tier1_2Longset = scan('tier1and2/background-proteome/wilcox_query_longest.txt')
tier1_2_3Fraction = scan('tier1and2and3/background-proteome/wilcox_query_fraction.txt')
tier1_2_3Longset = scan('tier1and2and3/background-proteome/wilcox_query_longest.txt')

# test
wilcox.test(proteomeFraction, proteomeFraction)$p.value

# results
wilcox.test(rbpFraction, proteomeFraction)$p.value
wilcox.test(tier1Fraction, proteomeFraction)$p.value
wilcox.test(tier1_2Fraction, proteomeFraction)$p.value
wilcox.test(tier1_2_3Fraction, proteomeFraction)$p.value
wilcox.test(tier1Fraction, rbpFraction)$p.value
wilcox.test(tier1_2Fraction, rbpFraction)$p.value
wilcox.test(tier1_2_3Fraction, rbpFraction)$p.value

wilcox.test(rbpLongest, proteomeLongest)$p.value
wilcox.test(tier1Longest, proteomeLongest)$p.value
wilcox.test(tier1_2Longest, proteomeLongest)$p.value
wilcox.test(tier1_2_3Longest, proteomeLongest)$p.value
wilcox.test(tier1Longest, rbpLongest)$p.value
wilcox.test(tier1_2Longest, rbpLongest)$p.value
wilcox.test(tier1_2_3Longest, rbpLongest)$p.value