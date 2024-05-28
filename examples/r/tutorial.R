library(reticulate)
library(reshape)
library(ggplot2)

use_python("/usr/local/bin/python3")
rmelts <- import("meltsdynamic")

liquidus = rmelts$MELTSdynamic(1L)

pressure = 500.0
temperature = 1200.0

bulk = c(48.68, 1.01, 17.64, 0.89, 0.03, 7.59, 0.0, 9.10, 0.0, 0.0, 12.45, 2.65, 0.03, 0.08, 0.2, 0.0, 0.0, 0.0, 0.0)

liquidus$engine$setBulkComposition(bulk)
liquidus$engine$pressure = pressure

# Note that as of Nov 18th, 2020 temperature is in Celcius
liquidus$engine$temperature = temperature

liquidus$engine$findLiquidus()
temperature = liquidus$engine$temperature

liquidus$engine$setSystemProperties("Mode", "Fractionate Solids")

ptpath = liquidus$copyAndKeepOutput()

ptpath$engine$calcEquilibriumState(1L, 0L)

while (ptpath$engine$temperature >= 1000) {
  
  ptpath = ptpath$addNodeAfter()
  ptpath$engine$temperature = ptpath$engine$temperature - 3
  
  ptpath$engine$calcEquilibriumState(1L, 1L);
  ptpath$engine$status$message
  
}

liq_comp = ptpath$getListProperty('dispComposition', 'liquid1')
df_test <- as.data.frame(t(liq_comp))

# These are lower case but that could be fixed
colnames(df_test) <- ptpath$endMemberFormulas['bulk']$bulk

temps = ptpath$getListProperty('temperature')
df_test$temps <- temps

# melt is a reshape function - it has nothing to do with MELTS!
df_test_long <- melt(df_test, id=c("temps"))

ggplot(df_test_long) + geom_point(aes(x=temps, y=value, color=variable)) + scale_x_reverse() + scale_y_continuous(limits=c(0,100)) +
                labs(x='T C (system)', y='Wt %', title='Melt composition') + theme_bw() + theme(legend.position = c(0.1, 0.75))

# Can extend the dataframe with rbind e.g:
df_test_long <- rbind(df_test_long, data.frame(temps = df_test$temps, variable='mass', value = ptpath$getListProperty('mass', 'liquid1')))

ggplot(subset(df_test_long, df_test_long$variable %in% c("sio2", "feo", "mgo", "cao", "na2o"))) + 
  geom_point(aes(x=temps, y=value, color=variable)) + scale_x_reverse() + scale_y_continuous(limits=c(0,100)) + 
  scale_color_manual(values=c("#56B4E9", "#E69F00", "#660099", "#999999", "#FFCC99")) + 
  labs(x=expression("T "*~degree*"C (system)"), y='Wt %', title='Melt composition') + theme_bw() + theme(legend.position = c(0.1, 0.75))
