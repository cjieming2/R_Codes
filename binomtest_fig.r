setwd("C:/Users/JM/thesis/mark_work/FIG/allelicSNPs")

x=binom.test(727,55341,0.013171341,"two.sided");
x$p.value
x$conf.int
x=binom.test(917,66560,0.013171341,"two.sided");
x$p.value
x$conf.int
x=binom.test(15803,1209789,0.013171341,"two.sided");
x$p.value
x$conf.int
x=binom.test(590,42099,0.013171341,"two.sided");
x$p.value
x$conf.int
x=binom.test(2251,143024,0.013171341,"two.sided");
x$p.value
x$conf.int
x=binom.test(16984,1322116,0.013171341,"two.sided");
x$p.value
x$conf.int

# AS<non-AS
x=binom.test(727,55341,0.013777043,"less");
x$p.value
x$conf.int
x=binom.test(590,42099,0.015738617,"less");
x$p.value
x$conf.int

