library(bnlearn)
library(MCMCpack)
library(dplyr)
library(tidyr)
library(purrr)
library(mice)
library(magrittr)
library(car)
library(stats)
library(ggplot2)
library(stringr)
library(gtools)
library(igraph)
library(graph)
library(Rgraphviz)
SEM_MICE_R=function(data,original_BN,mech,seed, proportion){
  set.seed(seed+3)   #create missing
  simulated_missdata = ampute(data, prop=proportion, mech=mech)%>% pluck("amp")
  names <- names(simulated_missdata)
  simulated_missdata=simulated_missdata%<>%mutate(across(all_of(names),as.factor))
  set.seed(seed+4)   #use SEM and compare
  check_structure_EM_miss=structural.em(simulated_missdata,maximize="tabu",
                                        max.iter =7,maximize.args = list(score = "bde"),debug = T)
  c1=compare(original_BN ,check_structure_EM_miss)%>% unlist()
  #use MICE 
  imputed_mice_data_miss=mice(simulated_missdata,m=3,maxit=3,printFlag=T,
                              method = 'polyreg', seed = seed+5) %>% complete()
  set.seed(seed+6)
  BN_structure_mice_miss=bnlearn::tabu(imputed_mice_data_miss,score="bde")
  c2=compare(original_BN,BN_structure_mice_miss)%>% unlist()
  cc=rbind(c1,c2) %>% as.data.frame()
  rownames(cc)=c("SEM", "MICE")
  return(cc) }
##################################################
result_R=function(mechanism="MCAR",miss_prop=0.15,data, original_BN){
  z=data.frame()
  for( i in 1:100){
    seed=i
    a=SEM_MICE_R(data=data,original_BN=original_BN,mech=mechanism,
                 seed=seed, proportion=miss_prop)
    z[i,1]=a[1,1]/(a[1,1]+a[1,2])
    z[i,2]=a[1,1]/(a[1,1]+a[1,3])
    z[i,3]=(2*z[i,1]*z[i,2])/(z[i,1]+z[i,2])
    z[i,4]=a[2,1]/(a[2,1]+a[2,2])
    z[i,5]=a[2,1]/(a[2,1]+a[2,3])
    z[i,6]=(2*z[i,4]*z[i,5])/(z[i,4]+z[i,5])  }
  colnames(z)=c("precision","recall","f_m","precision","recall","f_m")
  return(z) }
#########      ASIA      #####
asiadag= model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")
gr=graphviz.plot(asiadag,layout = "dot",shape = "ellipse")
node.attrs = nodeRenderInfo(gr)
node.attrs$col["A"] = "#b3e2cd"
node.attrs$fill["A"] = "#b3e2cd"
node.attrs$col["S"] = "#b3e2cd"
node.attrs$fill["S"] = "#b3e2cd"
node.attrs$col["T"] = "#b3e2cd"
node.attrs$fill["T"] = "#b3e2cd"
node.attrs$col["L"] = "#b3e2cd"
node.attrs$fill["L"] = "#b3e2cd"
node.attrs$col["E"] = "#b3e2cd"
node.attrs$fill["E"] = "#b3e2cd"
node.attrs$col["B"] = "#b3e2cd"
node.attrs$fill["B"] = "#b3e2cd"
node.attrs$col["X"] = "#b3e2cd"
node.attrs$fill["X"] = "#b3e2cd"
node.attrs$col["D"] = "#b3e2cd"
node.attrs$fill["D"] = "#b3e2cd"
nodeRenderInfo(gr) = node.attrs
renderGraph(gr)
#SEM and MICE on ASIA
asia_0.15_MCAR=result_R(mechanism="MCAR",miss_prop=0.15,data=asia,asiadag)
asia_0.3_MCAR=result_R(mechanism="MCAR",miss_prop=0.3,data=asia,asiadag)
asia_0.15_MAR=result_R(mechanism="MAR",miss_prop=0.15,data=asia,asiadag)
asia_0.3_MAR=result_R(mechanism="MAR",miss_prop=0.3,data=asia,asiadag)
asia_0.15_MNAR=result_R(mechanism="MNAR",miss_prop=0.15,data=asia,asiadag)
asia_0.3_MNAR=result_R(mechanism="MNAR",miss_prop=0.3,data=asia,asiadag)

method=c(rep(c("SEM","MICE"),each=100,times=2))
prop=c(rep(c("0.15","0.3"),each=200))
precision_MCAR=data.frame()
precision_MCAR$precision=rbind(asia_0.15_MCAR[,1],asia_0.15_MCAR[,4],
                               asia_0.3_MCAR[,1],asia_0.3_MCAR[,4])
precision_MCAR$method=method
precision_MCAR$prop=prop
recall_MCAR=data.frame()
recall_MCAR$recall=rbind(asia_0.15_MCAR[,2],asia_0.15_MCAR[,5],
                         asia_0.3_MCAR[,2],asia_0.3_MCAR[,5])
recall_MCAR$method=method
recall_MCAR$prop=prop
f_m_MCAR=data.frame()
f_m_MCAR$f_m=rbind(asia_0.15_MCAR[,3],asia_0.15_MCAR[,6],
                   asia_0.3_MCAR[,3],asia_0.3_MCAR[,6])
f_m_MCAR$method=method
f_m_MCAR$prop=prop
precision_MAR=data.frame()
precision_MAR$precision=rbind(asia_0.15_MAR[,1],asia_0.15_MAR[,4],
                              asia_0.3_MAR[,1],asia_0.3_MAR[,4])
precision_MAR$method=method
precision_MAR$prop=prop
recall_MAR=data.frame()
recall_MAR$recall=rbind(asia_0.15_MAR[,2],asia_0.15_MAR[,5],
                        asia_0.3_MAR[,2],asia_0.3_MAR[,5])
recall_MAR$method=method
recall_MAR$prop=prop
f_m_MAR=data.frame()
f_m_MAR$f_m=rbind(asia_0.15_MAR[,3],asia_0.15_MAR[,6],
                  asia_0.3_MAR[,3],asia_0.3_MAR[,6])
f_m_MAR$method=method
f_m_MAR$prop=prop
precision_MNAR=data.frame()
precision_MNAR$precision=rbind(asia_0.15_MNAR[,1],asia_0.15_MNAR[,4],
                               asia_0.3_MNAR[,1],asia_0.3_MNAR[,4])
precision_MNAR$method=method
precision_MNAR$prop=prop
recall_MNAR=data.frame()
recall_MNAR$recall=rbind(asia_0.15_MNAR[,2],asia_0.15_MNAR[,5],
                         asia_0.3_MNAR[,2],asia_0.3_MNAR[,5])
recall_MNAR$method=method
recall_MNAR$prop=prop
f_m_MNAR=data.frame()
f_m_MNAR$f_m=rbind(asia_0.15_MNAR[,3],asia_0.15_MNAR[,6],
                   asia_0.3_MNAR[,3],asia_0.3_MNAR[,6])
f_m_MNAR$method=method
f_m_MNAR$prop=prop
#######   average network  #######
SEM_output_networks1=list(NA,100)
MICE_output_networks1=list(NA,100)
for(i in 1:100){
  seed=i
  set.seed(seed+3)   #create missing
  simulated_missdata= ampute(asia, prop=0.15, mech="MAR")%>% pluck("amp")
  names = names(simulated_missdata)
  simulated_missdata=simulated_missdata%<>%mutate(across(all_of(names),as.factor))
  set.seed(seed+4)   #use SEM and compare
  sem = structural.em(simulated_missdata, maximize ="tabu",return.all =T,
                      max.iter =5,maximize.args = list(score = "bde"),debug = T)
  hamming_sem[i]=sum(asia!=check_structure_EM_miss$imputed)
  SEM_output_networks1[[i]]=sem$dag
  imputed_mice_data_miss=mice(simulated_missdata,m=3,maxit=3,printFlag=T,
                              method = 'polyreg', seed = seed+5) %>% complete()
  set.seed(seed+6)
  mice<- bnlearn::tabu(imputed_mice_data_miss, score = "bde")
  hamming_mice_0.15_mar[i]=sum(asia!=imputed_mice_data_miss)
  MICE_output_networks1[[i]]=mice$arcs  }
nodes = colnames(asia)
# Get strength and relative probability for each arc & direction
arcs_strength_SEM1=custom.strength(SEM_output_networks1,nodes=nodes,cpdag=F)
arcs_strength_MICE1=custom.strength(MICE_output_networks1,nodes=nodes,cpdag=F)
average_network_sem1 <- averaged.network(arcs_strength_SEM1)
gr1=graphviz.plot(average_network_sem1,fontsize=10,layout="neato",shape="ellipse")
average_network_mice1 <- averaged.network(arcs_strength_MICE1)
gr2=graphviz.plot(average_network_mice1,fontsize=10,layout="neato",shape="ellipse")
node.attrs = nodeRenderInfo(gr1)
node.attrs$col["A"] = "pink"
node.attrs$fill["A"] = "pink"
node.attrs$col["S"] = "pink"
node.attrs$fill["S"] = "pink"
node.attrs$col["T"] = "pink"
node.attrs$fill["T"] = "pink"
node.attrs$col["L"] = "pink"
node.attrs$fill["L"] = "pink"
node.attrs$col["E"] = "pink"
node.attrs$fill["E"] = "pink"
node.attrs$col["B"] = "pink"
node.attrs$fill["B"] = "pink"
node.attrs$col["X"] = "pink"
node.attrs$fill["X"] = "pink"
node.attrs$col["D"] = "pink"
node.attrs$fill["D"] = "pink"
nodeRenderInfo(gr1) = node.attrs
renderGraph(gr1)
##########repreat  the same code with mech="MNAR"
##   0.15_MCAR
hamming_sem_0.15_mcar=c()
hamming_mice_0.15_mcar=c() 
for( i in 1:100){
  seed=i
  set.seed(seed+3)   #create missing
  simulated_missdata=ampute(asia, prop=0.15, mech="MCAR")%>%pluck("amp")
  names = names(simulated_missdata)
  simulated_missdata=simulated_missdata%<>%mutate(across(all_of(names),as.factor))
  set.seed(seed+4)   #use SEM and compare
  check_structure_EM_miss=structural.em(simulated_missdata,maximize ="tabu",
                                        return.all =T, max.iter=4,maximize.args =list(score = "bde"),debug = T)
  hamming_sem_0.15_mcar[i]=sum(asia!=check_structure_EM_miss$imputed)
  imputed_mice_data_miss=mice(simulated_missdata, m=3, maxit = 3,printFlag=T,
                              method = 'polyreg', seed = seed+5) %>% complete()
  set.seed(seed+6)
  BN_structure_mice_miss <- bnlearn::tabu(imputed_mice_data_miss,score ="bde")
  hamming_mice_0.15_mcar[i]=sum(asia!=imputed_mice_data_miss)  }
##  0.3_MAR  ####################
hamming_sem_0.3_mar=c()
hamming_mice_0.3_mar=c() 
for( i in 1:100){
  seed=i
  set.seed(seed+3)   #create missing
  simulated_missdata= ampute(asia, prop=0.3, mech="MAR")%>% pluck("amp")
  names <- names(simulated_missdata)
  simulated_missdata=simulated_missdata%<>%mutate(across(all_of(names),as.factor))
  set.seed(seed+4)   #use SEM and compare
  check_structure_EM_miss= structural.em(simulated_missdata, maximize ="tabu"
                                         ,return.all =T,max.iter=4,maximize.args =list(score = "bde"),debug = T)
  hamming_sem_0.3_mar[i]=sum(data!=check_structure_EM_miss$imputed)
  imputed_mice_data_miss=mice(simulated_missdata,m=3,maxit=3,printFlag=T,
                              method = 'polyreg', seed = seed+5) %>% complete()
  set.seed(seed+6)
  BN_structure_mice_miss=bnlearn::tabu(imputed_mice_data_miss,score ="bde")
  hamming_mice_0.3_mar[i]=sum(data!=imputed_mice_data_miss)}
##0.3_MCAR####################
hamming_sem_0.3_mcar=c()
hamming_mice_0.3_mcar=c() 
for( i in 1:100){
  seed=i
  set.seed(seed+3)   #create missing
  simulated_missdata=ampute(asia, prop=0.3, mech="MCAR")%>%pluck("amp")
  names <- names(simulated_missdata)
  simulated_missdata=simulated_missdata%<>%mutate(across(all_of(names),as.factor))
  set.seed(seed+4)   #use SEM and compare
  check_structure_EM_miss=structural.em(simulated_missdata, maximize ="tabu",
                                        return.all =T, max.iter=4,maximize.args =list(score = "bde"),debug = T)
  hamming_sem_0.3_mcar[i]=sum(asia!=check_structure_EM_miss$imputed)
  imputed_mice_data_miss=mice(simulated_missdata,m=3,maxit=3,printFlag=T,
                              method = 'polyreg', seed = seed+5) %>% complete()
  set.seed(seed+6)
  BN_structure_mice_miss <- bnlearn::tabu(imputed_mice_data_miss,score="bde")
  hamming_mice_0.3_mcar[i]=sum(asia!=imputed_mice_data_miss)}
##     0.3_MNAR  ####################
hamming_sem_0.3_mnar=c()
hamming_mice_0.3_mnar=c() 
for( i in 1:100){
  seed=i
  set.seed(seed+3)   #create missing
  simulated_missdata= ampute(asia, prop=0.3, mech="MNAR")%>% pluck("amp")
  names= names(simulated_missdata)
  simulated_missdata=simulated_missdata%<>%mutate(across(all_of(names),as.factor))
  set.seed(seed+4)   #use SEM and compare
  check_structure_EM_miss <- structural.em(simulated_missdata, maximize = "tabu"
                                           ,return.all =T,max.iter=4,maximize.args =list(score = "bde"),debug = T)
  hamming_sem_0.3_mnar[i]=sum(data!=check_structure_EM_miss$imputed)
  imputed_mice_data_miss <- mice(simulated_missdata, m=3, maxit = 3,printFlag=T,
                                 method = 'polyreg', seed = seed+5) %>% complete()
  set.seed(seed+6)
  BN_structure_mice_miss <- bnlearn::tabu(imputed_mice_data_miss, score = "bde")
  hamming_mice_0.3_mnar[i]=sum(data!=imputed_mice_data_miss) }
##                       hamming barplot                   ###############
library(rstatix)
library(ggpubr)
library(gridExtra)
method=c(rep(c("SEM","MICE"),each=100))
ham_0.15_mcar=data.frame()
ham_0.15_mcar$hamming=rbind(hamming_sem_0.15_mcar,hamming_mice_0.15_mcar)
ham_0.15_mcar$method=method                   
stat.test1 <- ham_0.15_mcar %>%
  wilcox_test(hamming~method) %>%
  add_significance("p")
bp_ham_0.15_MCAR <- ggbarplot(ham_0.15_mcar, x = "method", y = "hamming",
                              add="mean_sd",xlab="",color="black",fill="method",palette=c("#00AFBB","yellow"),
                              position=position_dodge(0.8),title="MCAR 0.15")+theme(legend.position="none") 
stat.test1=stat.test1 %>% add_xy_position(fun="mean_sd",x="method",dodge =0.8) 
bp_ham_0.15_MCAR=bp_ham_0.15_MCAR + 
  stat_pvalue_manual(stat.test1,label = "p.signif",tip.length=0.01,hide.ns=T)
#####################################
ham_0.15_mar=data.frame()
ham_0.15_mar$hamming=rbind(hamming_sem_0.15_mar,hamming_mice_0.15_mar)
ham_0.15_mar$method=method 
stat.test2 <- ham_0.15_mar %>%
  wilcox_test(hamming~method) %>%
  add_significance("p")
bp_ham_0.15_MAR <- ggbarplot(ham_0.15_mar, x = "method",y="hamming", 
                             add="mean_sd",xlab ="",color="black",fill="method",palette=c("#00AFBB","yellow"),
                             position=position_dodge(0.8),title ="MAR 0.15")+theme(legend.position="none")
stat.test2=stat.test2 %>% add_xy_position(fun="mean_sd",x="method",dodge=0.8) 
bp_ham_0.15_MAR=bp_ham_0.15_MAR + 
  stat_pvalue_manual(stat.test2,label="p.signif",tip.length=0.01,hide.ns=T)
#####################################
ham_0.15_mnar=data.frame()
ham_0.15_mnar$hamming=rbind(hamming_sem_0.15_mnar,hamming_mice_0.15_mnar)
ham_0.15_mnar$method=method 
stat.test3 <- ham_0.15_mnar %>%
  wilcox_test(hamming~method) %>%
  add_significance("p")
bp_ham_0.15_MNAR <- ggbarplot(ham_0.15_mnar, x="method", y="hamming", 
                              add="mean_sd",color="black",fill="method",palette=c("#00AFBB","yellow"),xlab ="",
                              position=position_dodge(0.8),title ="MNAR 0.15")+theme(legend.position ="none") 
stat.test3=stat.test3 %>% add_xy_position(fun="mean_sd",x="method",dodge=0.8) 
bp_ham_0.15_MNAR=bp_ham_0.15_MNAR + 
  stat_pvalue_manual(stat.test3,label="p.signif",tip.length=0.01,hide.ns=T)
###############################################
ham_0.3_mcar=data.frame()
ham_0.3_mcar$hamming=rbind(hamming_sem_0.3_mcar,hamming_mice_0.3_mcar)
ham_0.3_mcar$method=method 
stat.test4 <- ham_0.3_mcar %>%
  wilcox_test(hamming~method) %>%
  add_significance("p")
bp_ham_0.3_MCAR <- ggbarplot(ham_0.3_mcar, x="method", y="hamming",
                             add="mean_sd",color="black",fill="method",palette=c("#00AFBB","yellow"),xlab ="",
                             position=position_dodge(0.8),title="MCAR 0.3")+theme(legend.position="none")
stat.test4=stat.test4 %>% add_xy_position(fun="mean_sd",x="method",dodge=0.8) 
bp_ham_0.3_MCAR=bp_ham_0.3_MCAR + 
  stat_pvalue_manual(stat.test4,label="p.signif",tip.length=0.01,hide.ns=T)
##############################################
ham_0.3_mar=data.frame()
ham_0.3_mar$hamming=rbind(hamming_sem_0.3_mar,hamming_mice_0.3_mar)
ham_0.3_mar$method=method 
stat.test5 <- ham_0.3_mar %>%
  wilcox_test(hamming~method) %>%
  add_significance("p")
bp_ham_0.3_MAR <- ggbarplot(ham_0.3_mar, x="method",y="hamming",
                            add="mean_sd",color="black",fill="method",palette=c("#00AFBB","yellow"),xlab ="",
                            position = position_dodge(0.8),title="MAR 0.3")+ theme(legend.position = "none") 
stat.test5=stat.test5 %>% add_xy_position(fun="mean_sd",x="method",dodge=0.8) 
bp_ham_0.3_MAR=bp_ham_0.3_MAR + 
  stat_pvalue_manual(stat.test5,label="p.signif",tip.length=0.01,hide.ns=T)
#####################################
ham_0.3_mnar=data.frame()
ham_0.3_mnar$hamming=rbind(hamming_sem_0.3_mnar,hamming_mice_0.3_mnar)
ham_0.3_mnar$method=method 
stat.test6 <- ham_0.3_mnar %>%
  wilcox_test(hamming~method) %>%
  add_significance("p")
bp_ham_0.3_MNAR <- ggbarplot(ham_0.3_mnar, x="method",y="hamming", 
                             add="mean_sd",color="black",fill="method",palette=c("#00AFBB","yellow"),xlab ="",
                             position=position_dodge(0.8),title ="MNAR 0.3")+theme(legend.position="none")
stat.test6=stat.test6 %>% add_xy_position(fun ="mean_sd",x="method",dodge =0.8) 
bp_ham_0.3_MNAR=bp_ham_0.3_MNAR + 
  stat_pvalue_manual(stat.test6,label="p.signif", tip.length=0.01,hide.ns=T)
hamming_plot=grid.arrange(bp_ham_0.15_MCAR,bp_ham_0.3_MCAR,bp_ham_0.15_MAR
                          ,bp_ham_0.3_MAR ,bp_ham_0.15_MNAR,bp_ham_0.3_MNAR,ncol=2)
###############################################
#precision_MCAR 
stat.test1 <- precision_MCAR %>%
  group_by(prop) %>%
  wilcox_test(precision~method) %>%
  add_significance("p")
bp_prec_MCAR=ggbarplot(precision_MCAR,x ="prop",y="precision",add ="mean_sd", 
                       color= "black",fill = "method", palette =c("lightgreen","plum1"), 
                       position=position_dodge(0.8),title = "MCAR")+
  scale_y_continuous(limits = c(0, 1.5), breaks = c(1), labels = c("1"))
stat.test1=stat.test1%>% add_xy_position(fun="mean_sd",x="prop",dodge =0.8) 
bp_prec_MCAR=bp_prec_MCAR+
stat_pvalue_manual(stat.test1,label="p.signif",tip.length=0.01,hide.ns=T)
################################################
#precision_MAR 
stat.test2 <- precision_MAR %>%
  group_by(prop) %>%
  wilcox_test(precision~method) %>%
  add_significance("p")
bp_prec_MAR=ggbarplot(precision_MAR,x ="prop",y="precision",add ="mean_sd", 
                      color= "black",fill = "method", palette =c("lightgreen","plum1"),
                      position = position_dodge(0.8),title = "MAR")+   theme(legend.position = "none")+
  scale_y_continuous(limits = c(0, 1.5), breaks = c(1), labels = c("1"))  
stat.test2=stat.test2 %>% add_xy_position(fun ="mean_sd",x="prop",dodge=0.8) 
bp_prec_MAR=bp_prec_MAR + 
  stat_pvalue_manual(stat.test2,label = "p.signif",tip.length = 0.01,hide.ns=T)
##################################################
#precision_MNAR 
stat.test3 <- precision_MNAR %>%
  group_by(prop) %>%
ilcox_test(precision~method) %>%
  add_significance("p")
bp_prec_MNAR=ggbarplot(precision_MNAR ,x ="prop",y="precision",add ="mean_sd", 
                       color= "black",fill = "method", palette =c("lightgreen","plum1"),
                       position = position_dodge(0.8),title = "MNAR")+   theme(legend.position = "none")+                        
  scale_y_continuous(limits = c(0, 1.5), breaks = c(1), labels = c("1"))  
stat.test3=stat.test3 %>% add_xy_position(fun ="mean_sd", x ="prop",dodge =0.8) 
bp_prec_MNAR=bp_prec_MNAR+
  stat_pvalue_manual(stat.test3,label="p.signif",tip.length=0.01,hide.ns=T)
###################################################
#recall_MCAR 
stat.test4 <- recall_MCAR %>%
  group_by(prop) %>%
  wilcox_test(recall~method) %>%
  add_significance("p")
bp_rec_MCAR=ggbarplot(recall_MCAR, x = "prop",y ="recall",add ="mean_sd", 
                      color= "black",fill = "method", palette =c("lightgreen","plum1"), 
                      position=position_dodge(0.8),title = "MCAR")+theme(legend.position ="none")+                       
  scale_y_continuous(limits = c(0, 1.5), breaks = c(1), labels = c("1"))
stat.test4=stat.test4 %>% add_xy_position(fun ="mean_sd",x="prop",dodge =0.8) 
bp_rec_MCAR=bp_rec_MCAR+
stat_pvalue_manual(stat.test4,label ="p.signif",tip.length = 0.01,hide.ns=T)
####################################################
#recall_MAR 
stat.test5 <- recall_MAR %>%
  group_by(prop) %>%
  wilcox_test(recall~method) %>%
  add_significance("p")
bp_rec_MAR=ggbarplot(recall_MAR, x = "prop", y ="recall",add ="mean_sd", 
                     color= "black",fill = "method", palette =c("lightgreen","plum1"), 
                     position=position_dodge(0.8),title ="MAR")+theme(legend.position ="none")+        
  scale_y_continuous(limits=c(0, 1.5),breaks=c(1),labels=c("1"))
stat.test5=stat.test5 %>% add_xy_position(fun="mean_sd",x="prop",dodge =0.8) 
bp_rec_MAR=bp_rec_MAR+
  stat_pvalue_manual(stat.test5,label="p.signif",tip.length=0.01,hide.ns=T)
####################################################
#recall_MNAR 
stat.test6 <- recall_MNAR %>%
  group_by(prop) %>%
  wilcox_test(recall~method) %>%
  add_significance("p")
bp_rec_MNAR=ggbarplot(recall_MNAR,x="prop",y="recall",add ="mean_sd", 
                      color= "black",fill = "method", palette =c("lightgreen","plum1"), 
                      position=position_dodge(0.8),title ="MNAR")+theme(legend.position ="none")+                
  scale_y_continuous(limits = c(0, 1.5), breaks = c(1), labels = c("1")) 
stat.test6=stat.test6 %>% add_xy_position(fun ="mean_sd",x="prop",dodge=0.8) 
bp_rec_MNAR=bp_rec_MNAR+, 
stat_pvalue_manual(stat.test6,label ="p.signif"tip.length = 0.01,hide.ns=T)
#####################################################
#f_measure_MCAR 
stat.test7 <- f_measure_MCAR %>%
  group_by(prop) %>%
  wilcox_test(f_measure~method) %>%
  add_significance("p")
bp_f_MCAR=ggbarplot(f_measure_MCAR,x ="prop",y ="f_measure",add="mean_sd", 
                    color= "black",fill = "method", palette =c("lightgreen","plum1"), 
                    position=position_dodge(0.8),title = "MCAR")+ theme(legend.position ="none")+ 
  scale_y_continuous(limits = c(0, 1.5), breaks = c(1), labels = c("1"))
stat.test7=stat.test7 %>% add_xy_position(fun ="mean_sd", x="prop",dodge=0.8) 
bp_f_MCAR=bp_f_MCAR+,
stat_pvalue_manual(stat.test7,label= "p.signif"tip.length = 0.01,hide.ns=T)
####################################################
#f_measure_MAR 
stat.test8 <- f_measure_MAR %>%
  group_by(prop) %>%
  wilcox_test(f_measure~method) %>%
  add_significance("p")
bp_f_MAR=ggbarplot(f_measure_MAR, x="prop",y="f_measure",add="mean_sd", 
                   color= "black",fill = "method", palette =c("lightgreen","plum1"), 
                   position = position_dodge(0.8),title = "MAR")+theme(legend.position ="none")+                        
  scale_y_continuous(limits = c(0, 1.5), breaks = c(1), labels = c("1"))
stat.test8=stat.test8 %>% add_xy_position(fun="mean_sd",x="prop",dodge =0.8) 
bp_f_MAR=bp_f_MAR+
stat_pvalue_manual(stat.test8,label="p.signif",tip.length = 0.01,hide.ns=T)
##################################################
#f_measure_MNAR 
stat.test9 <- f_measure_MNAR %>%
  group_by(prop) %>%
  wilcox_test(f_measure~method) %>%
  add_significance("p")
bp_f_MNAR=ggbarplot(f_measure_MNAR,x="prop",y="f_measure",add="mean_sd", 
                    color= "black",fill = "method", palette =c("lightgreen","plum1"), 
                    position = position_dodge(0.8),title = "MNAR")+  theme(legend.position = "none")+                      
  scale_y_continuous(limits = c(0, 1.5), breaks = c(1), labels = c("1"))
stat.test9=stat.test9 %>% add_xy_position(fun ="mean_sd",x="prop",dodge=0.8) 
bp_f_MNAR=bp_f_MNAR + stat_pvalue_manual(stat.test9,label = "p.signif",
                                         tip.length = 0.01,hide.ns=T)
asia_plot=grid.arrange(bp_prec_MCAR,bp_rec_MCAR,bp_f_MCAR, bp_prec_MAR,
                       bp_rec_MAR,bp_f_MAR,bp_prec_MNAR,bp_rec_MNAR,bp_f_MNAR, ncol = 3) 

################ synthetic data   ###########################
randomBN=function(nodes, maxin=5,numlevels=4,alphaparent=1,alphanone=1, ...){
  # ic-dag: Ide's and Cozman's Generating Multi-connected DAGs algorithm
  if(is.numeric(nodes)){
    numnodes <- nodes
    ranBN <- random.graph(as.character(1:numnodes), method = "ic-dag", 
                          max.in.degree = maxin, ... = ...)
  }else {
    numnodes <- length(nodes)
    ranBN <- random.graph(nodes, method="ic-dag",max.in.degree=maxin, ...=...)
  }
  # now get the CPTs for each node
  CPTs <- list()
  for(i in 1:numnodes){
    indegree <- in.degree(ranBN, nodes(ranBN)[i])
    # first need to get the arrays for the dirichlet
    parents <- rep(alphaparent, numlevels)
    noparents <- rep(alphanone, numlevels)
    if(indegree == 0 ) { # no parents
      CPTs[[ nodes(ranBN)[i] ]] <- rdirichlet(1,noparents) 
    }
    else{ # has parents
      tempdirichlet <- rdirichlet(numlevels^indegree,parents)
      transdirichlet <- t(tempdirichlet)
      dimsdirichlet <- rep(numlevels, indegree+1)
      dim(transdirichlet) <- dimsdirichlet
      CPTs[[ nodes(ranBN)[i] ]] <- transdirichlet
    }
  }
  # and now make the parameterised BN
  ranBNparam <- custom.fit(ranBN, dist = CPTs)
  BNinfo <- list("structure" = ranBN, "complete" = ranBNparam)
  return(BNinfo)
}
##############################################
SEM_MICE=function(nvar,npoints,seed,mech, proportion){
  nodes <- letters[seq( from = 1, to = nvar )]
  set.seed(seed)  #create BN
  original_BN=andomBN(nodes,maxin=3,numlevels=3,alphaparent=0.5,alphanone=5)
  original_BN_structure <- original_BN %>% pluck("structure")
  original_BN_complete <- original_BN %>% pluck("complete")
  set.seed(seed+1)  #generate data from BN
  simulated_data=rbn(original_BN_complete,npoints)%>%
  sapply(as.integer)%>%as.data.frame()
  fac_cols=colnames(simulated_data) # Covert to factor
  simulated_data=simulated_data %>% mutate(across(all_of(fac_cols), as.factor))
  set.seed(seed+3)   #create missing
  simulated_missdata=ampute(simulated_data,prop=proportion,mech=mech)%>%
  pluck("amp") 
  names <- names(simulated_missdata)
  simulated_missdata=simulated_missdata%<>% mutate(across(all_of(names),as.factor))
  set.seed(seed+4)   #use SEM and compare
  check_structure_EM_miss=structural.em(simulated_missdata, maximize="tabu",
                                        max.iter=7,maximize.args = list(score = "bde"),debug = T)
  c1=compare(original_BN $structure,check_structure_EM_miss)%>% unlist()
  imputed_mice_data_miss= mice(simulated_missdata, m=3, maxit = 3, printFlag=T,
                               method = 'polyreg', seed = seed+5) %>% complete()
  set.seed(seed+6)
  BN_structure_mice_miss=bnlearn::tabu(imputed_mice_data_miss, score = "bde")
  c2=compare(original_BN $structure,BN_structure_mice_miss)%>% unlist()
  cc=rbind(c1,c2) %>% as.data.frame()
  rownames(cc)=c("SEM", "MICE")
  return(cc)
}
##################################################
result=function(mechanism="MCAR",miss_prop=0.15,nvariables=10,datapoints=5000){
  z=data.frame()
  for( i in 1:100){
    seeds=i
    a=SEM_MICE(mech=mechanism,nvar=nvariables,npoints=datapoints
               ,proportion=miss_prop,seed=seeds)
    z[i,1]=a[1,1]/(a[1,1]+a[1,2])
    z[i,2]=a[1,1]/(a[1,1]+a[1,3])
    z[i,3]=(2*z[i,1]*z[i,2])/(z[i,1]+z[i,2])
    z[i,4]=a[2,1]/(a[2,1]+a[2,2])
    z[i,5]=a[2,1]/(a[2,1]+a[2,3])
    z[i,6]=(2*z[i,4]*z[i,5])/(z[i,4]+z[i,5])
  }
  colnames(z)=c("precision","recall","f_m","precision","recall","f_m")
  return(z) }
###############################
MCAR_0.15_5_1000=result("MCAR",0.15,5,1000)
MAR_0.15_5_1000=result("MAR",0.15,5,1000)
MNAR_0.15_5_1000=result("MNAR",0.15,5,1000)
MCAR_0.15_10_1000=result("MCAR",0.15,10,1000)
MAR_0.15_10_1000=result("MAR",0.15,10,1000)
MNAR_0.15_10_1000=result("MNAR",0.15,10,1000)
MCAR_0.15_20_1000=result("MCAR",0.15,20,1000)
MAR_0.15_20_1000=result("MAR",0.15,20,1000)
MNAR_0.15_20_1000=result("MNAR",0.15,20,1000)

MCAR_0.15_5_5000=result("MCAR",0.15,5,5000)
MAR_0.15_5_5000=result("MAR",0.15,5,5000)
MNAR_0.15_5_5000=result("MNAR",0.15,5,5000)
MCAR_0.15_10_5000=result("MCAR",0.15,10,5000)
MAR_0.15_10_5000=result("MAR",0.15,10,5000)
MNAR_0.15_10_5000=result("MNAR",0.15,10,5000)
MCAR_0.15_20_5000=result("MCAR",0.15,20,5000)
MAR_0.15_20_5000=result("MAR",0.15,20,5000)
MNAR_0.15_20_5000=result("MNAR",0.15,20,5000)

MCAR_0.3_5_1000=result("MCAR",0.3,5,1000)
MAR_0.3_5_1000=result("MAR",0.3,5,1000)
MNAR_0.3_5_1000=result("MNAR",0.3,5,1000)
MCAR_0.3_10_1000=result("MCAR",0.3,10,1000)
MAR_0.3_10_1000=result("MAR",0.3,10,1000)
MNAR_0.3_10_1000=result("MNAR",0.3,10,1000)
MCAR_0.3_20_1000=result("MCAR",0.3,20,1000)
MAR_0.3_20_1000=result("MAR",0.3,20,1000)
MNAR_0.3_20_1000=result("MNAR",0.3,20,1000)

MCAR_0.3_5_5000=result("MCAR",0.3,5,5000)
MAR_0.3_5_5000=result("MAR",0.3,5,5000)
MNAR_0.3_5_5000=result("MNAR",0.3,5,5000)
MCAR_0.3_10_5000=result("MCAR",0.3,10,5000)
MAR_0.3_10_5000=result("MAR",0.3,10,5000)
MNAR_0.3_10_5000=result("MNAR",0.3,10,5000)
MCAR_0.3_20_5000=result("MCAR",0.3,20,5000)
MAR_0.3_20_5000=result("MAR",0.3,20,5000)
MNAR_0.3_20_5000=result("MNAR",0.3,20,5000)

MCAR_0.5_5_1000=result("MCAR",0.5,5,1000)
MAR_0.5_5_1000=result("MAR",0.5,5,1000)
MNAR_0.5_5_1000=result("MNAR",0.5,5,1000)
MCAR_0.5_10_1000=result("MCAR",0.5,10,1000)
MAR_0.5_10_1000=result("MAR",0.5,10,1000)
MNAR_0.5_10_1000=result("MNAR",0.5,10,1000)
MCAR_0.5_20_1000=result("MCAR",0.5,20,1000)
MAR_0.5_20_1000=result("MAR",0.5,20,1000)
MNAR_0.5_20_1000=result("MNAR",0.5,20,1000)

MCAR_0.5_5_5000=result("MCAR",0.5,5,5000)
MAR_0.5_5_5000=result("MAR",0.5,5,5000)
MNAR_0.5_5_5000=result("MNAR",0.5,5,5000)
MCAR_0.5_10_5000=result("MCAR",0.5,10,5000)
MAR_0.5_10_5000=result("MAR",0.5,10,5000)
MNAR_0.5_10_5000=result("MNAR",0.5,10,5000)
MCAR_0.5_20_5000=result("MCAR",0.5,20,5000)
MAR_0.5_20_5000=result("MAR",0.5,20,5000)
MNAR_0.5_20_5000=result("MNAR",0.5,20,5000)

method=c(rep(c("SEM","MICE"),each=100,times=3))
missing=c(rep(c("0.15","0.3","0.5"),each=200))

MCAR_5_1000=data.frame()
MCAR_5_1000$precision=rbind(MCAR_0.15_5_1000[,1],MCAR_0.15_5_1000[,4],
                            MCAR_0.3_5_1000[,1],MCAR_0.3_5_1000[,4],MCAR_0.5_5_1000[,1],MCAR_0.5_5_1000[,4])
MCAR_5_1000$recall=rbind(MCAR_0.15_5_1000[,2],MCAR_0.15_5_1000[,5],
                         MCAR_0.3_5_1000[,2],MCAR_0.3_5_1000[,5],MCAR_0.5_5_1000[,2],MCAR_0.5_5_1000[,5])
MCAR_5_1000$f_m=rbind(MCAR_0.15_5_1000[,3],MCAR_0.15_5_1000[,6],
                      MCAR_0.3_5_1000[,3],MCAR_0.3_5_1000[,6],MCAR_0.5_5_1000[,3],MCAR_0.5_5_1000[,6])
MCAR_5_1000$method=method
MCAR_5_1000$missing=missing

MAR_5_1000=data.frame()
MAR_5_1000$precision=rbind(MAR_0.15_5_1000[,1],MAR_0.15_5_1000[,4],
                           MAR_0.3_5_1000[,1],MAR_0.3_5_1000[,4],MAR_0.5_5_1000[,1],MAR_0.5_5_1000[,4])
MAR_5_1000$recall=rbind(MAR_0.15_5_1000[,2],MAR_0.15_5_1000[,5],
                        MAR_0.3_5_1000[,2],MAR_0.3_5_1000[,5],MAR_0.5_5_1000[,2],MAR_0.5_5_1000[,5])
MAR_5_1000$f_m=rbind(MAR_0.15_5_1000[,3],MAR_0.15_5_1000[,6],
                     MAR_0.3_5_1000[,3],MAR_0.3_5_1000[,6],MAR_0.5_5_1000[,3],MAR_0.5_5_1000[,6])
MAR_5_1000$method=method
MAR_5_1000$missing=missing

MNAR_5_1000=data.frame()
MNAR_5_1000$precision=rbind(MNAR_0.15_5_1000[,1],MNAR_0.15_5_1000[,4],
                            MNAR_0.3_5_1000[,1],MNAR_0.3_5_1000[,4],MNAR_0.5_5_1000[,1],MNAR_0.5_5_1000[,4])
MNAR_5_1000$recall=rbind(MNAR_0.15_5_1000[,2],MNAR_0.15_5_1000[,5],
                         MNAR_0.3_5_1000[,2],MNAR_0.3_5_1000[,5],MNAR_0.5_5_1000[,2],MNAR_0.5_5_1000[,5])
MNAR_5_1000$f_m=rbind(MNAR_0.15_5_1000[,3],MNAR_0.15_5_1000[,6],
                      MNAR_0.3_5_1000[,3],MNAR_0.3_5_1000[,6],MNAR_0.5_5_1000[,3],MNAR_0.5_5_1000[,6])
MNAR_5_1000$method=method
MNAR_5_1000$missing=missing

#precision_MAR_5_1000
stat.test1 <- MCAR_5_1000 %>%
  group_by(missing) %>%
  wilcox_test(precision~method)  %>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "****",p <= 0.01  ~ "***",
    p <= 0.05  ~ "**",p <= 0.1  ~ "*",TRUE      ~ "ns"))
stat.test1
bp_prec_MCAR_5_1000=ggbarplot(MCAR_5_1000,x="missing",y ="precision", 
                              add ="mean_sd", color= "black",fill ="method", palette=c("#ff9966","#00cc99"),  
                              position = position_dodge(0.8),title = "MCAR")+
  scale_y_continuous(limits = c(0, 1.5), breaks = c(1), labels = c("1"))+
  theme(legend.position="none")
stat.test1=stat.test1 %>% add_xy_position(fun="mean_sd",x="missing",dodge=0.8) 
bp_prec_MCAR_5_1000=bp_prec_MCAR_5_1000+
stat_pvalue_manual(stat.test1,label="p.signif",tip.length =0.01,hide.ns=T)
################################################
#recall_MCAR_5_1000
stat.test2 <- MCAR_5_1000 %>%
  group_by(missing) %>%
  wilcox_test(recall~method)%>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "****",p <= 0.01  ~ "***",
    p <= 0.05  ~ "**",p <= 0.1  ~ "*",TRUE      ~ "ns"))
bp_rec_MCAR_5_1000 <- ggbarplot(MCAR_5_1000, x = "missing", y = "recall",
                                add = "mean_sd",color= "black",fill = "method", palette =c("#ff9966","#00cc99"), 
                                position = position_dodge(0.8),title = "MCAR")+theme(legend.position = "none")+                                
  scale_y_continuous(limits = c(0, 1.5), breaks = c(1), labels = c("1"))
stat.test2=stat.test2 %>% add_xy_position(fun = "mean_sd",x="missing",dodge=0.8) 
bp_rec_MCAR_5_1000=bp_rec_MCAR_5_1000 +
  stat_pvalue_manual(stat.test2,label = "p.signif",tip.length = 0.01,hide.ns=T)
##################################################
#f_measure_MCAR_5_1000
stat.test3 <- MCAR_5_1000 %>%
  group_by(missing) %>%
  wilcox_test(f_m~method)%>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "****",p <= 0.01  ~ "***",
    p <= 0.05  ~ "**",p <= 0.1  ~ "*",TRUE      ~ "ns"))
bp_f_MCAR_5_1000 <- ggbarplot(MCAR_5_1000, x = "missing", y = "f_m",  
                              add = "mean_sd",color= "black",fill = "method", palette =c("#ff9966","#00cc99"), 
                              position = position_dodge(0.8),title = "MCAR")+theme(legend.position = "none")+  
  scale_y_continuous(limits = c(0, 1.5), breaks = c(1), labels = c("1"))+
  theme(legend.position = "none")
stat.test3=stat.test3 %>% add_xy_position(fun="mean_sd",x="missing",dodge=0.8) 
bp_f_MCAR_5_1000=bp_f_MCAR_5_1000 + 
  stat_pvalue_manual(stat.test3,label = "p.signif",tip.length = 0.01,hide.ns=T)
####################################################
#precision_MAR_5_1000 
stat.test4 <- MAR_5_1000 %>%
  group_by(missing) %>%
  wilcox_test(precision~method)%>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "****",p <= 0.01  ~ "***",
    p <= 0.05  ~ "**",p <= 0.1  ~ "*",TRUE      ~ "ns"))
bp_prec_MAR_5_1000 <- ggbarplot(MAR_5_1000, x = "missing", y = "precision",  
                                add = "mean_sd",color= "black",fill = "method", palette =c("#ff9966","#00cc99"),  
                                position = position_dodge(0.8),title = "MAR")+theme(legend.position = "none")
scale_y_continuous(limits = c(0, 1.5), breaks = c(1), labels = c("1"))+
  stat.test4=tat.test4 %>% add_xy_position(fun="mean_sd",x="missing",dodge=0.8) 
bp_prec_MAR_5_1000=bp_prec_MAR_5_1000 +  
stat_pvalue_manual(stat.test4,label = "p.signif", tip.length =0.01,hide.ns=T)
###################################################
#recall_MAR_5_1000
stat.test5 <- MAR_5_1000 %>%
  group_by(missing) %>%
  wilcox_test(recall~method)%>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "****",p <= 0.01  ~ "***",
    p <= 0.05  ~ "**",p <= 0.1  ~ "*",TRUE      ~ "ns"))
bp_rec_MAR_5_1000 <- ggbarplot(MAR_5_1000, x = "missing", y = "recall", 
                               add = "mean_sd", color= "black",fill = "method", palette =c("#ff9966","#00cc99"), 
                               position = position_dodge(0.8),title = "MAR")+theme(legend.position = "none")+   
  scale_y_continuous(limits = c(0, 1.5), breaks = c(1), labels = c("1"))+
  theme(legend.position = "none")
stat.test5=stat.test5 %>% add_xy_position(fun="mean_sd",x="missing",dodge=0.8) 
bp_rec_MAR_5_1000=bp_rec_MAR_5_1000 + 
  stat_pvalue_manual(stat.test5,label = "p.signif",tip.length = 0.01,hide.ns=T)
####################################################
#f_measure_MAR_5_1000
stat.test6 <- MAR_5_1000 %>%
  group_by(missing) %>%
  wilcox_test(f_m~method)%>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "****",p <= 0.01  ~ "***",
    p <= 0.05  ~ "**",p <= 0.1  ~ "*",TRUE      ~ "ns"))
bp_f_MAR_5_1000 <- ggbarplot(MAR_5_1000, x = "missing", y = "f_m",
                             add = "mean_sd",color= "black",fill = "method", palette =c("#ff9966","#00cc99"), 
                             position = position_dodge(0.8),title = "MAR")+theme(legend.position = "none")+  
  scale_y_continuous(limits = c(0, 1.5), breaks = c(1), labels = c("1"))
stat.test6=stat.test6 %>% add_xy_position(fun="mean_sd",x="missing",dodge=0.8) 
bp_f_MAR_5_1000=bp_f_MAR_5_1000 + 
  stat_pvalue_manual(stat.test6,label = "p.signif",tip.length = 0.01,hide.ns=T)
####################################################
#precision_MNAR_5_1000 
stat.test7 <- MNAR_5_1000 %>%
  group_by(missing) %>%
  wilcox_test(precision~method)%>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "****",p <= 0.01  ~ "***",
    p <= 0.05  ~ "**",p <= 0.1  ~ "*",TRUE      ~ "ns"))
bp_prec_MNAR_5_1000 <- ggbarplot(MNAR_5_1000, x = "missing", y = "precision",  
                                 add = "mean_sd",color= "black",fill = "method", palette =c("#ff9966","#00cc99"),  
                                 position = position_dodge(0.8),title = "MNAR")+ theme(legend.position = "none")
scale_y_continuous(limits = c(0, 1.5), breaks = c(1), labels = c("1"))
stat.test7=stat.test7 %>% add_xy_position(fun="mean_sd",x="missing",dodge=0.8) 
bp_prec_MNAR_5_1000=bp_prec_MNAR_5_1000 +  
stat_pvalue_manual(stat.test7,label = "p.signif", tip.length =0.01,hide.ns=T)
####################################################
#recall_MNAR_5_1000
stat.test8 <- MNAR_5_1000 %>%
  group_by(missing) %>%
  wilcox_test(recall~method)%>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "****",p <= 0.01  ~ "***",
    p <= 0.05  ~ "**",p <= 0.1  ~ "*",TRUE      ~ "ns"))
bp_rec_MNAR_5_1000 <- ggbarplot(MNAR_5_1000, x = "missing", y = "recall",  
                                add = "mean_sd",color= "black",fill = "method", palette =c("#ff9966","#00cc99"),
                                position = position_dodge(0.8),title = "MNAR")+theme(legend.position = "none")+  
  scale_y_continuous(limits = c(0, 1.5), breaks = c(1), labels = c("1"))
stat.test8=stat.test8 %>% add_xy_position(fun="mean_sd",x="missing",dodge=0.8) 
bp_rec_MNAR_5_1000=bp_rec_MNAR_5_1000 +
  stat_pvalue_manual(stat.test8,label = "p.signif",tip.length = 0.01,hide.ns=T)
####################################################
#f_measure_MNAR_5_1000
stat.test9 <- MNAR_5_1000 %>%
  group_by(missing) %>%
  wilcox_test(f_m~method)%>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "****",p <= 0.01  ~ "***",
    p <= 0.05  ~ "**",p <= 0.1  ~ "*",TRUE      ~ "ns"))
bp_f_MNAR_5_1000 <- ggbarplot(MNAR_5_1000, x = "missing", y = "f_m", 
                              add = "mean_sd",color= "black",fill ="method",palette =c("#ff9966","#00cc99"), 
                              position = position_dodge(0.8),title ="MNAR")+theme(legend.position = "none")+  
  scale_y_continuous(limits = c(0, 1.5), breaks = c(1), labels = c("1"))
stat.test9=stat.test9 %>% add_xy_position(fun="mean_sd", x="missing",dodge=0.8) 
bp_f_MNAR_5_1000=bp_f_MNAR_5_1000 +
  stat_pvalue_manual(stat.test9,label = "p.signif",tip.length = 0.01,hide.ns=T)
plot=ggarrange(bp_prec_MCAR_5_1000,bp_rec_MCAR_5_1000,bp_f_MCAR_5_1000, 
               bp_prec_MAR_5_1000, bp_rec_MAR_5_1000,
               bp_f_MAR_5_1000,bp_prec_MNAR_5_1000,
               bp_rec_MNAR_5_1000,bp_f_MNAR_5_1000, common.legend=TRUE, legend="right")
## the same code for the rest of the plots.just the number of
##variables and observations need to be changed. 
###         hamming synthetic data       #############
ham_0.3_mcar=data.frame()
hs_0.3_mcar_20_5000=c()
hm_0.3_mcar_20_5000=c()
for( i in 1:100){
  seed=i
  nodes <- letters[seq( from = 1, to = 20 )]
  set.seed(seed)  #create BN
  original_BN=randomBN(nodes,maxin=3,numlevels=3,alphaparent=0.5,alphanone=5)
  original_BN_structure <- original_BN %>% pluck("structure")
  original_BN_complete <- original_BN %>% pluck("complete")
  set.seed(seed+1)  #generate data from BN
  simulated_data=rbn(original_BN_complete,5000)%>%sapply(as.integer)%>%
  as.data.frame()
  fac_cols <- colnames(simulated_data) # Covert to factor
  simulated_data <- simulated_data %>% mutate(across(all_of(fac_cols), as.factor))
  set.seed(seed+3)   #create missing
  simulated_missdata=ampute(simulated_data,prop=0.5,mech="MNAR")%>%
  pluck("amp")
  names <- names(simulated_missdata)
  simulated_missdata=simulated_missdata%<>%mutate(across(all_of(names),as.factor))
  set.seed(seed+4)   #use SEM and compare
  check_structure_EM_miss=structural.em(simulated_missdata,maximize="tabu"
                                        ,max.iter=5,return.all =T,maximize.args = list(score = "bde"),debug = T)
  hs_0.3_mcar_20_5000[i]=sum(simulated_data!=check_structure_EM_miss$imputed)
  imputed_mice_data_miss <- mice(simulated_missdata, m=3, maxit = 3,
                                 printFlag=T,method = 'polyreg', seed = seed+5) %>% complete()
  hm_0.3_mcar_20_5000[i]=sum(simulated_data!=imputed_mice_data_miss)
  set.seed(seed+6)
  BN_structure_mice_miss=bnlearn::tabu(imputed_mice_data_miss,score="bde")
}
ham_0.3_mcar$hamming=c(hs_0.3_mcar_20_5000,hm_0.3_mcar_20_5000)
ham_0.3_mcar$method=c(rep(c("SEM","MICE")),each=100)
#the same code for the rest.just need to change the mechanism and missing prop
################  hamming plots ##############
stat.test1 <- ham_0.15_mcar %>%
  wilcox_test(hamming~method) %>%
  add_significance("p")
bp_ham_0.3_MCAR <- ggbarplot(ham_0.15_mcar, x = "method", y ="hamming", ,
                             add="mean_sd",xlab="",color="black",fill="method",palette=c("pink","#cbd5e8"),
                             position=position_dodge(0.8),title="MCAR 0.15")+theme(legend.position="none")  
stat.test1=stat.test1%>%add_xy_position(fun="mean_sd",x="method",dodge=0.8) 
bp_ham_0.3_MCAR=bp_ham_0.3_MCAR + 
  stat_pvalue_manual(stat.test1,label = "p.signif", tip.length = 0.01,hide.ns=T)
#####################################
stat.test2 <- ham_0.3_mar %>%
  wilcox_test(hamming~method) %>%
  add_significance("p")
bp_ham_0.3_MAR <- ggbarplot(ham_0.3_mar, x = "method", y = "hamming", 
                            add = "mean_sd",xlab="",color="black",fill ="method",palette=c("pink","#cbd5e8"),
                            position=position_dodge(0.8),title="MAR 0.15")+theme(legend.position="none")  
stat.test2=stat.test2 %>% add_xy_position(fun="mean_sd",x="method",dodge=0.8) 
bp_ham_0.3_MAR=bp_ham_0.3_MAR + 
  stat_pvalue_manual(stat.test2,label = "p.signif", tip.length = 0.01,hide.ns=T)
#####################################
stat.test3 <- ham_0.3_mnar %>%
  wilcox_test(hamming~method) %>%
  add_significance("p")
bp_ham_0.3_MNAR <- ggbarplot(ham_0.3_mnar, x = "method", y = "hamming",  
                             add ="mean_sd",color="black",fill="method",palette=c("pink","#cbd5e8"),xlab ="",
                             position = position_dodge(0.8),title="MNAR 0.15")+theme(legend.position="none")  
stat.test3=stat.test3%>%add_xy_position(fun="mean_sd",x="method",dodge=0.8) 
bp_ham_0.3_MNAR=bp_ham_0.3_MNAR + 
  stat_pvalue_manual(stat.test3,label = "p.signif" , tip.length = 0.01,hide.ns=T)
###############################################
stat.test4 <- ham_0.5_mcar %>%
  wilcox_test(hamming~method) %>%
  add_significance("p")
bp_ham_0.5_MCAR <- ggbarplot(ham_0.5_mcar, x = "method", y = "hamming",  
                             add="mean_sd",color="black",fill="method",palette=c("pink","#cbd5e8"),xlab ="",
                             position=position_dodge(0.8),title="MCAR 0.3")+theme(legend.position="none")  
stat.test4=stat.test4%>%add_xy_position(fun="mean_sd",x="method",dodge=0.8) 
bp_ham_0.5_MCAR=bp_ham_0.5_MCAR + 
  stat_pvalue_manual(stat.test4,label = "p.signif", tip.length = 0.01,hide.ns=T)
#################################################
stat.test5 <- ham_0.5_mar %>%
  wilcox_test(hamming~method) %>%
  add_significance("p")
bp_ham_0.5_MAR <- ggbarplot(ham_0.5_mar, x = "method", y = "hamming",  
                            add="mean_sd",color="black",fill="method",palette=c("pink","#cbd5e8"),xlab ="",
                            position=position_dodge(0.8),title = "MAR 0.3")+theme(legend.position = "none") 
stat.test5=stat.test5%>%add_xy_position(fun="mean_sd",x="method",dodge=0.8) 
bp_ham_0.5_MAR=bp_ham_0.5_MAR + 
  stat_pvalue_manual(stat.test5,label = "p.signif", tip.length = 0.01,hide.ns=T)
#####################################
stat.test6 <- ham_0.5_mnar %>%
  wilcox_test(hamming~method) %>%
  add_significance("p")
bp_ham_0.5_MNAR <- ggbarplot(ham_0.5_mnar, x = "method", y = "hamming",  
                             add="mean_sd",color="black",fill="method",palette=c("pink","#cbd5e8"),xlab ="",
                             position=position_dodge(0.8),title="MNAR 0.3")+theme(legend.position="none")  
stat.test6=stat.test6%>%add_xy_position(fun="mean_sd",x="method",dodge =0.8) 
bp_ham_0.5_MNAR=bp_ham_0.5_MNAR + 
  stat_pvalue_manual(stat.test6,label = "p.signif",tip.length = 0.01,hide.ns=T)
hamming_plot=grid.arrange(bp_ham_0.3_MCAR,bp_ham_0.5_MCAR,
                          bp_ham_0.3_MAR,bp_ham_0.5_MAR,bp_ham_0.3_MNAR,bp_ham_0.5_MNAR,ncol=2)
###      HRS     dataset       ##############
#HRSdata.csv
data4=read.csv(choose.files())
dim(data4)
for( i in 2:ncol(data4)){
  data4[,i]=factor(data4[,i])
}
#missing prop in each column
apply(data4,2,function(x)round(sum(is.na(x))/nrow(data4),2))
#missing percentage in data
sum(is.na(data4[,-1]))/(nrow(data4)*ncol(data4[,-1]))*100
sum(complete.cases(data4))
data4[complete.cases(data4),]
library(VIM)
missing_plot=aggr(data4[,-1], col=c('blue','orange'),gap=0.2,cex.axis=0.53,bars=T,
                  numbers=T, sortVars=F,sortCombs = TRUE,digits=15,
                  labels=names(data4[,-1]),  ylab=c("Missing data","Pattern"))
miss_prop=data.frame(variables=names(data4[-1])
                     ,missing=round((missing_plot$missings[,2]/nrow(data4)),4))
barplot(miss_prop$missing,names.arg = names(data4[-1])
        ,col="orange",las=2,cex.names= 0.7,ylab = "missing proportion",cex.lab =1.4)
# Get 50 SEM outputs
set.seed(10001)
seed_numbers <- sample.int(10000, 50)
SEM_output_networks = lapply(seed_numbers, function(number) {
  set.seed(number)
  structural.em(data4[,-1],maximize="tabu",maximize.args=list(score="bde"),debug=T)
})
a=list()
b=list()
for(i in 1:100){
  set.seed(seed_numbers[i])
  a[[i]]=structural.em(data4[,-1],maximize="tabu",
                       maximize.args=list(score="bde"),debug=T)
  mice_data <- mice(data4[,-1], m=3, maxit = 3, printFlag=T,
                    method = 'polyreg', seed = seed_numbers[i]+5) %>% complete()
  set.seed(seed_numbers[i]+6)
  b[[i]]= bnlearn::tabu(mice_data, score = "bde")
}
nodes = colnames(data4[,-1])
# Get strength and relative probability for each arc & direction
arcs_strength = custom.strength(a, nodes = nodes, cpdag = F)
arcs_strength2 = custom.strength(b, nodes = nodes, cpdag = F)
# Get averaged network
average_network <- averaged.network(arcs_strength)
average_network2 <- averaged.network(arcs_strength2)
attr(arcs_strength, "threshold")
attr(arcs_strength2, "threshold")
shd(average_network,average_network2,debug = T)
graphviz.compare(average_network,average_network2, layout = "fdp",
                 shape = "rectangle",fontsize = 22 )
av1=data.frame(average_network$arcs)
av2=data.frame(average_network2$arcs)
#same arcs
intersect(av1,av2)
#different arcs
anti_join(av1,av2,by=c("from","to"))
anti_join(av2,av1,by=c("from","to"))
gR <- graphviz.plot(average_network,fontsize = 22,layout = "fdp")
gR2 <- graphviz.plot(average_network2,fontsize = 22,layout = "fdp")                  
node.attrs = nodeRenderInfo(gR)
node.attrs$col["gender"] = "#f4cae4"
node.attrs$fill["gender"] = "#f4cae4"
node.attrs$col["drinking"] = "#f4cae4"
node.attrs$fill["drinking"] = "#f4cae4"
node.attrs$col["partnership"] = "#f4cae4"
node.attrs$fill["partnership"] = "#f4cae4"
node.attrs$col["hhincome"] = "#f4cae4"
node.attrs$fill["hhincome"] = "#f4cae4"
node.attrs$col["education"] = "#f4cae4"
node.attrs$fill["education"] = "#f4cae4"
node.attrs$col["exercise"] = "#f4cae4"
node.attrs$fill["exercise"] = "#f4cae4"
node.attrs$col["BMI"] = "#f4cae4"
node.attrs$fill["BMI"] = "#f4cae4"
node.attrs$col["arthritis"] = "#f4cae4"
node.attrs$fill["arthritis"] = "#f4cae4"
node.attrs$col["memory"] = "#f4cae4"
node.attrs$fill["memory"] = "#f4cae4"
node.attrs$col["race"] = "#f4cae4"
node.attrs$fill["race"] = "#f4cae4"
node.attrs$col["TICS_M"] = "#f4cae4"
node.attrs$fill["TICS_M"] = "#f4cae4"
node.attrs$col["age"] = "#f4cae4"
node.attrs$fill["age"] = "#f4cae4"
node.attrs$col["idincome"] = "#f4cae4"
node.attrs$fill["idincome"] = "#f4cae4"
node.attrs$col["smoking"] = "#f4cae4"
node.attrs$fill["smoking"] = "#f4cae4"
node.attrs$col["lung"] = "#f4cae4"
node.attrs$fill["lung"] = "#f4cae4"
node.attrs$col["cancer"] = "#f4cae4"
node.attrs$fill["cancer"] = "#f4cae4"
node.attrs$col["cancer_treatment"] = "#f4cae4"
node.attrs$fill["cancer_treatment"] = "#f4cae4"
node.attrs$col["diabetes_treatment"] = "#f4cae4"
node.attrs$fill["diabetes_treatment"] = "#f4cae4"
node.attrs$col["diabetes"] = "#f4cae4"
node.attrs$fill["diabetes"] = "#f4cae4"
node.attrs$col["hbp_treatment"] = "#f4cae4"
node.attrs$fill["hbp_treatment"] = "#f4cae4"
node.attrs$col["hbp"] = "#f4cae4"
node.attrs$fill["hbp"] = "#f4cae4"
node.attrs$col["cholestrol_treatment"] = "#f4cae4"
node.attrs$fill["cholestrol_treatment"] = "#f4cae4"
node.attrs$col["heart"] = "#f4cae4"
node.attrs$fill["heart"] = "#f4cae4"
node.attrs$col["heart_treatment"] = "#f4cae4"
node.attrs$fill["heart_treatment"] = "#f4cae4"
node.attrs$col["cholestrol"] = "#f4cae4"
node.attrs$fill["cholestrol"] = "#f4cae4"
node.attrs$col["stroke"] = "#f4cae4"
node.attrs$fill["stroke"] = "#f4cae4"
node.attrs$col["stroke_treatment"] = "#f4cae4"
node.attrs$fill["stroke_treatment"] = "#f4cae4"
nodeRenderInfo(gR) = node.attrs
renderGraph(gR)
# output to pdf file
pdf(file = "SEM_HRS.pdf",
    width = 10,
    height = 7)
renderGraph(gR)
dev.off()
#########################################
node.attrs2 = nodeRenderInfo(gR2)
node.attrs2$col["gender"] = "#cbd5e8"
node.attrs2$fill["gender"] = "#cbd5e8"
node.attrs2$col["drinking"] = "#cbd5e8"
node.attrs2$fill["drinking"] = "#cbd5e8"
node.attrs2$col["partnership"] = "#cbd5e8"
node.attrs2$fill["partnership"] = "#cbd5e8"
node.attrs2$col["hhincome"] = "#cbd5e8"
node.attrs2$fill["hhincome"] = "#cbd5e8"
node.attrs2$col["education"] = "#cbd5e8"
node.attrs2$fill["education"] = "#cbd5e8"
node.attrs2$col["exercise"] = "#cbd5e8"
node.attrs2$fill["exercise"] = "#cbd5e8"
node.attrs2$col["BMI"] = "#cbd5e8"
node.attrs2$fill["BMI"] = "#cbd5e8"
node.attrs2$col["arthritis"] = "#cbd5e8"
node.attrs2$fill["arthritis"] = "#cbd5e8"
node.attrs2$col["memory"] = "#cbd5e8"
node.attrs2$fill["memory"] = "#cbd5e8"
node.attrs2$col["race"] = "#cbd5e8"
node.attrs2$fill["race"] = "#cbd5e8"
node.attrs2$col["TICS_M"] = "#cbd5e8"
node.attrs2$fill["TICS_M"] = "#cbd5e8"
node.attrs2$col["age"] = "#cbd5e8"
node.attrs2$fill["age"] = "#cbd5e8"
node.attrs2$col["idincome"] = "#cbd5e8"
node.attrs2$fill["idincome"] = "#cbd5e8"
node.attrs2$col["smoking"] = "#cbd5e8"
node.attrs2$fill["smoking"] = "#cbd5e8"
node.attrs2$col["lung"] = "#cbd5e8"
node.attrs2$fill["lung"] = "#cbd5e8"
node.attrs2$col["cancer"] = "#cbd5e8"
node.attrs2$fill["cancer"] = "#cbd5e8"
node.attrs2$col["cancer_treatment"] = "#cbd5e8"
node.attrs2$fill["cancer_treatment"] = "#cbd5e8"
node.attrs2$col["diabetes_treatment"] = "#cbd5e8"
node.attrs2$fill["diabetes_treatment"] = "#cbd5e8"
node.attrs2$col["diabetes"] = "#cbd5e8"
node.attrs2$fill["diabetes"] = "#cbd5e8"
node.attrs2$col["hbp_treatment"] = "#cbd5e8"
node.attrs2$fill["hbp_treatment"] = "#cbd5e8"
node.attrs2$col["hbp"] = "#cbd5e8"
node.attrs2$fill["hbp"] = "#cbd5e8"
node.attrs2$col["cholestrol_treatment"] = "#cbd5e8"
node.attrs2$fill["cholestrol_treatment"] = "#cbd5e8"
node.attrs2$col["heart"] = "#cbd5e8"
node.attrs2$fill["heart"] = "#cbd5e8"
node.attrs2$col["heart_treatment"] = "#cbd5e8"
node.attrs2$fill["heart_treatment"] = "#cbd5e8"
node.attrs2$col["cholestrol"] = "#cbd5e8"
node.attrs2$fill["cholestrol"] = "#cbd5e8"
node.attrs2$col["stroke"] = "#cbd5e8"
node.attrs2$fill["stroke"] = "#cbd5e8"
node.attrs2$col["stroke_treatment"] = "#cbd5e8"
node.attrs2$fill["stroke_treatment"] = "#cbd5e8"
nodeRenderInfo(gR2) = node.attrs2
renderGraph(gR2)