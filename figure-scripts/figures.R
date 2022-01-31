rm(list=ls())
source(here::here("0-config.R"))
library(cowplot)
library(patchwork)
theme_set(theme_ki())

#load spline data
H1_spline <- readRDS(here("figure-data/H1_adj_spline_data.RDS"))
H2_spline <- readRDS(here("figure-data/H2_adj_spline_data.RDS"))
H3_spline <- readRDS(here("figure-data/H3_adj_spline_data.RDS"))
H4_spline <- readRDS(here("figure-data/H4_adj_spline_data.RDS"))

#load results for quartiles
H1_quartiles <- readRDS(here("results/adjusted/H1_adj_res.RDS"))
H2_quartiles <- readRDS(here("results/adjusted/H2_adj_res.RDS"))
H3_quartiles <- readRDS(here("results/adjusted/H3_adj_res.RDS"))
H4_quartiles <- readRDS(here("results/adjusted/H4_adj_res.RDS"))

d_for_plot <- function(x_name, y_name, x_var, y_var, spline, quart){
  d <- NULL
  for (i in 1: length(x_var)) {
    for (j in 1:length(y_var)){
      exists <- (quart%>%filter(X==x_var[i], Y==y_var[j]) %>% nrow()) != 0
      if (exists){
        new <- data.frame(x=x_name[i], y=y_name[j], spline%>%filter(Xvar==x_var[i], Yvar==y_var[j]), quart%>%filter(X==x_var[i], Y==y_var[j])%>%select(q1, q3))
        d <- rbind(d, new)
      }
    }
  }
  d
}

color_levels <- c("Sum of WHO motor milestones", "Comprehension", "Expressive language",
                  "Communication", "Gross motor", "Personal social", "Combined")

d1 <- d_for_plot(c("Ln CRP Year 1", "Ln CRP Year 2"),
                 c("Sum of WHO motor milestones", 
                   "Comprehension", "Expressive language",
                   "Comprehension", "Expressive language",
                   "Communication", "Gross motor", "Personal social", "Combined"), 
                 c("t2_ln_crp", "t3_ln_crp"), 
                 c("sum_who", 
                   "z_cdi_und_t2", "z_cdi_say_t2",
                   "z_cdi_und_t3", "z_cdi_say_t3",
                   "z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq"),
                 H1_spline %>% rbind(H2_spline), H1_quartiles %>% rbind(H2_quartiles))

d1$x <- d1$x %>% as.factor()
d1$y <- factor(d1$y, levels=c("Sum of WHO motor milestones", 'Comprehension', "Expressive language",
                              "Communication", "Gross motor", "Personal social", "Combined"))

plot_facet <- function(d, theme, title, x_label, type) {
  if (type=="solid") {
    p <- d %>% ggplot(aes(x=X))+
      geom_smooth(aes(y = fit, color=y, linetype="Year 1"), se = F) 
  } else {
    p <- d %>% ggplot(aes(x=X))+
      geom_smooth(aes(y = fit, color=y, linetype="Year 2"), se = F)
  }
  p <- p +
    geom_vline(aes(xintercept=q1), size=.5, color="grey30", linetype="dashed") +
    geom_vline(aes(xintercept=q3), size=.5, color="grey30", linetype="dashed") +
    #geom_point(aes(y=Y), alpha=0.5) +
    scale_colour_manual(values=tableau10[c(1:7)], limits=color_levels) + 
    scale_fill_manual(values=tableau10[c(1:7)], limits=color_levels) + 
    scale_linetype_manual("Outcome measurement",values=c("Year 1"=1,"Year 2"=6)) +
    geom_ribbon(aes(ymin=lwrS, ymax=uprS, fill=y, color=y), alpha=0.5) + 
    xlab(" ") + ylab(element_blank()) + 
    facet_grid(rows = vars(y), scales = "free_y") 

  
  if(!theme) {p <- p + theme(strip.text.y = element_blank())}
  if(!is.null(title)) {p <- p + ylab(title)}
  if(x_label) {p <- p + xlab(d$x[1])}
  
  p + theme(legend.position = "right") + guides(color="none", fill="none", linetype=guide_legend(override.aes = aes(color="black")))
}

fig <- function(width, heigth){
  options(repr.plot.width = width, repr.plot.height = heigth)
}

p <- plot_facet(d1 %>% filter(grepl("sum_who", Yvar)), F, "Sum of WHO motor\nmilestones", F, type = "solid")
p1 <- plot_facet(d1 %>% filter(grepl("z_cdi_.{3}_t2", Yvar)), F, "CDI Z-scores", T, type = "solid")
p2 <- plot_facet(d1 %>% filter(grepl("z_cdi_.{3}_t3", Yvar) & Xvar=="t2_ln_crp"), F, NULL, F, "longdash")
p3 <- plot_facet(d1 %>% filter(grepl("z_cdi_.{3}_t3", Yvar) & Xvar=="t3_ln_crp"), T, NULL, F, "longdash") 
p4 <- ggplot() + theme_nothing()
p5 <- plot_facet(d1 %>% filter(grepl("easq", Yvar) & Xvar=="t2_ln_crp"), F, "EASQ Z-scores", T, "longdash")
p6 <- plot_facet(d1 %>% filter(grepl("easq", Yvar) & Xvar=="t3_ln_crp"), T, NULL, T, "longdash")

ggpubr::ggarrange(p,ggplot()+theme_nothing(),ggplot()+theme_nothing(),
          p1,p2,p3, p4,p5,p6, common.legend = T,
          ncol = 3, nrow = 3,
          heights = c(1, 1.7, 3.5), widths=c(1, 1, 1)) %>% ggsave(filename = here("figures/splines_crp.jpg"), width = 10, height=15)


d2 <- d_for_plot(c("Ln AGP Year 1", "Ln AGP Year 2"),
                 c("Sum of WHO motor milestones", 
                   "Comprehension", "Expressive language",
                   "Comprehension", "Expressive language",
                   "Communication", "Gross motor", "Personal social", "Combined"), 
                 c("t2_ln_agp", "t3_ln_agp"), 
                 c("sum_who", 
                   "z_cdi_und_t2", "z_cdi_say_t2",
                   "z_cdi_und_t3", "z_cdi_say_t3",
                   "z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq"),
                 H1_spline %>% rbind(H2_spline), H1_quartiles %>% rbind(H2_quartiles))

d2$x <- d2$x %>% as.factor()
d2$y <- factor(d2$y, levels=c("Sum of WHO motor milestones", 'Comprehension', "Expressive language",
                              "Communication", "Gross motor", "Personal social", "Combined"))
p <- plot_facet(d2 %>% filter(grepl("sum_who", Yvar)), F, "Sum of WHO motor\nmilestones", F, type = "solid")
p1 <- plot_facet(d2 %>% filter(grepl("z_cdi_.{3}_t2", Yvar)), F, "CDI Z-scores", T, type = "solid")
p2 <- plot_facet(d2 %>% filter(grepl("z_cdi_.{3}_t3", Yvar) & Xvar=="t2_ln_agp"), F, NULL, F, "longdash")
p3 <- plot_facet(d2 %>% filter(grepl("z_cdi_.{3}_t3", Yvar) & Xvar=="t3_ln_agp"), T, NULL, F, "longdash") 
p4 <- ggplot() + theme_nothing()
p5 <- plot_facet(d2 %>% filter(grepl("easq", Yvar) & Xvar=="t2_ln_agp"), F, "EASQ Z-scores", T, "longdash")
p6 <- plot_facet(d2 %>% filter(grepl("easq", Yvar) & Xvar=="t3_ln_agp"), T, NULL, T, "longdash")

ggpubr::ggarrange(p,ggplot()+theme_nothing(),ggplot()+theme_nothing(),
                  p1,p2,p3, p4,p5,p6, common.legend = T,
                  ncol = 3, nrow = 3,
                  heights = c(1, 1.7, 3.5), widths=c(1, 1, 1)) %>% ggsave(filename = here("figures/splines_agp.jpg"), width = 10, height=15)


d3 <- d_for_plot(c("Th1/Th2 Year 1", "Th1/Th2 Year 2"),
                 c("Sum of WHO motor milestones", 
                   "Comprehension", "Expressive language",
                   "Comprehension", "Expressive language",
                   "Communication", "Gross motor", "Personal social", "Combined"), 
                 c("t2_ratio_th1_th2", "t3_ratio_th1_th2"), 
                 c("sum_who", 
                   "z_cdi_und_t2", "z_cdi_say_t2",
                   "z_cdi_und_t3", "z_cdi_say_t3",
                   "z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq"),
                 H1_spline %>% rbind(H2_spline), H1_quartiles %>% rbind(H2_quartiles))


d3$x <- d3$x %>% as.factor()
d3$y <- factor(d3$y, levels=c("Sum of WHO motor milestones", 'Comprehension', "Expressive language",
                              "Communication", "Gross motor", "Personal social", "Combined"))
p <- plot_facet(d3 %>% filter(grepl("sum_who", Yvar)), F, "Sum of WHO motor\nmilestones", F, type = "solid")
p1 <- plot_facet(d3 %>% filter(grepl("z_cdi_.{3}_t2", Yvar)), F, "CDI Z-scores", T, type = "solid")
p2 <- plot_facet(d3 %>% filter(grepl("z_cdi_.{3}_t3", Yvar) & Xvar=="t2_ratio_th1_th2"), F, NULL, F, "longdash")
p3 <- plot_facet(d3 %>% filter(grepl("z_cdi_.{3}_t3", Yvar) & Xvar=="t3_ratio_th1_th2"), T, NULL, F, "longdash") 
p4 <- ggplot() + theme_nothing()
p5 <- plot_facet(d3 %>% filter(grepl("easq", Yvar) & Xvar=="t2_ratio_th1_th2"), F, "EASQ Z-scores", T, "longdash")
p6 <- plot_facet(d3 %>% filter(grepl("easq", Yvar) & Xvar=="t3_ratio_th1_th2"), T, NULL, T, "longdash")

ggpubr::ggarrange(p,ggplot()+theme_nothing(),ggplot()+theme_nothing(),
                  p1,p2,p3, p4,p5,p6, common.legend = T,
                  ncol = 3, nrow = 3, 
                  heights = c(1, 1.7, 3.5), widths=c(1, 1, 1)) %>% ggsave(filename = here("figures/splines_th1_th2.jpg"), width = 10, height=15)




d3 <- H1_quartiles %>% rbind(H2_quartiles) %>% filter(grepl("th1_th2", X))
d3 <- d3 %>% mutate(group=factor(ifelse(Y=="sum_who", "WHO motor milestones", ifelse(grepl("cdi", Y), "CDI Z-scores", "EASQ Z-scores")),
                                 levels=c("WHO motor milestones", "CDI Z-scores", "EASQ Z-scores")), 
                    Xvar = factor(ifelse(grepl("t2", X), "Th1/Th2 Year 1", "Th1/Th2 Year 2"), levels=c("Th1/Th2 Year 1", "Th1/Th2 Year 2")),
                    Yvar = case_when(grepl("who", Y) ~ "Sum of WHO motor milestones",
                                     grepl("und", Y) ~ "CDI comprehension Z-score",
                                     grepl("say", Y) ~ "CDI expressive language Z-score",
                                     grepl("comm", Y) ~ "EASQ communication Z-score",
                                     grepl("motor", Y) ~ "EASQ gross motor Z-score",
                                     grepl("personal", Y) ~ "EASQ personal social Z-score",
                                     grepl("combined", Y) ~ "EASQ combined Z-score"))
d3$Yvar <- factor(d3$Yvar, levels=c("Sum of WHO motor milestones", paste(c("CDI comprehension", "CDI expressive language",
                                   "EASQ communication", "EASQ gross motor", "EASQ personal social", "EASQ combined"), "Z-score")))
d3$Ytime <- factor(ifelse(grepl("t3|easq", d3$Y), "Year 2", "Year 1"))

p <- ggplot(d3, aes(x=Yvar, y=point.diff)) + 
  geom_pointrange(aes(ymin=lb.diff , ymax=ub.diff, color=group, group=Yvar, shape=Ytime),
                  position = position_dodge2(reverse=T,width = 0.5),
                  size = 1) +
  #geom_text(aes(label=ref), position = position_nudge(y = (abs(yrange[1])+abs(yrange[2]))/10)) +
  scale_x_discrete(limits = rev(levels(d3$Yvar))) +
  coord_flip() + 
  facet_grid(~Xvar) + 
  labs(y = "Adjusted difference in mean development outcome\nbetween 25th and 75th percentile of Th1/Th2 ratio", 
       x=element_blank()) +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_shape_manual(breaks =c("Year 1", "Year 2"), values=c(21,16)) +
  scale_colour_manual(values=tableau10) + 
  theme_ki() +
  theme(plot.title = element_text(hjust = 0),
        panel.spacing = unit(0, "lines"),
        legend.position = "right") +
  guides(color = guide_legend(title="Development measurement"), 
         shape=guide_legend(title="Development\nmeasurement time"))
  
p


HR <- readRDS(here("results/adjusted/HR_adj_res.RDS"))
HR_th1_th2 <- HR %>% filter(grepl("th1_th2", X) & !grepl("sit", Y))
HR_th1_th2$Yvar <- factor(c("Hands-and-knees crawling","Standing with assistance",
                     "Walking with assistance","Standing alone","Walking alone"),
                     levels=c("Hands-and-knees crawling","Standing with assistance","Walking with assistance","Standing alone","Walking alone"))
HR_th1_th2$Xvar <- factor(ifelse(grepl("t2", HR_th1_th2$X), "Th1/Th2 Year 1", "Th1/Th2 Year 2"), levels=c("Th1/Th2 Year 1", "Th1/Th2 Year 2"))

p1 <- ggplot(HR_th1_th2, aes(x=Yvar, y=point.HR, color="WHO motor milestones", shape="Year 1")) + 
  geom_pointrange(aes(ymin=lb.HR, ymax=ub.HR),
                  size = 1) +
  scale_x_discrete(limits = rev(levels(HR_th1_th2$Yvar))) +
  facet_grid(~Xvar, drop = T) +
  coord_flip() + 
  labs(y = "Adjusted hazard ratio of motor milestone achievement\nbetween 25th and 75th percentile of Th1/Th2 ratio", 
       x=element_blank()) +
  geom_hline(yintercept = 1, linetype="dashed") +
  scale_shape_manual(values=c(21,16)) +
  scale_colour_manual(values=tableau10) + 
  theme_ki() +
  theme(plot.title = element_text(hjust = 0),
        panel.spacing = unit(0, "lines"))

ggpubr::ggarrange(p1, p, widths = c(1.5, 2.5)) %>%  ggsave(filename = here("figures/gam_th1_th2.jpg"), width = 16, height=7)




d1 <- H1_quartiles %>% rbind(H2_quartiles) %>% filter(grepl("crp", X))
d1 <- d1 %>% mutate(group=factor(ifelse(Y=="sum_who", "WHO motor milestones", ifelse(grepl("cdi", Y), "CDI Z-scores", "EASQ Z-scores")),
                                 levels=c("WHO motor milestones", "CDI Z-scores", "EASQ Z-scores")), 
                    Xvar = factor(ifelse(grepl("t2", X), "Ln CRP Year 1", "Ln CRP Year 2"), levels=c("Ln CRP Year 1", "Ln CRP Year 2")),
                    Yvar = case_when(grepl("who", Y) ~ "Sum of WHO motor milestones",
                                     grepl("und", Y) ~ "CDI comprehension Z-score",
                                     grepl("say", Y) ~ "CDI expressive language Z-score",
                                     grepl("comm", Y) ~ "EASQ communication Z-score",
                                     grepl("motor", Y) ~ "EASQ gross motor Z-score",
                                     grepl("personal", Y) ~ "EASQ personal social Z-score",
                                     grepl("combined", Y) ~ "EASQ combined Z-score"))
d1$Yvar <- factor(d1$Yvar, levels=c("Sum of WHO motor milestones", paste(c("CDI comprehension", "CDI expressive language",
                                                                           "EASQ communication", "EASQ gross motor", "EASQ personal social", "EASQ combined"), "Z-score")))
d1$Ytime <- factor(ifelse(grepl("t3|easq", d1$Y), "Year 2", "Year 1"))

p <- ggplot(d1, aes(x=Yvar, y=point.diff)) + 
  geom_pointrange(aes(ymin=lb.diff , ymax=ub.diff, color=group, group=Yvar, shape=Ytime),
                  position = position_dodge2(reverse = TRUE, width = .5),
                  size = 1) +
  scale_x_discrete(limits = rev(levels(d1$Yvar)))+
  #geom_text(aes(label=ref), position = position_nudge(y = (abs(yrange[1])+abs(yrange[2]))/10)) +
  facet_grid(~Xvar, drop = T) +
  coord_flip() + 
  labs(y = "Adjusted difference in mean development outcome\nbetween 25th and 75th percentile of Th1/Th2 ratio", 
       x=element_blank()) +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_shape_manual(values=c(21,16)) +
  scale_colour_manual(values=tableau10) + 
  theme_ki() +
  theme(plot.title = element_text(hjust = 0),
        panel.spacing = unit(0, "lines"),
        legend.position = "right")+
  guides(color = guide_legend(title="Development measurement"), shape=guide_legend(title="Development\nmeasurement time", reverse = T))
p


HR_crp <- HR %>% filter(grepl("crp", X) & !grepl("sit", Y))
HR_crp$Yvar <- factor(c("Hands-and-knees crawling","Standing with assistance",
                            "Walking with assistance","Standing alone","Walking alone"),
                          levels=c("Hands-and-knees crawling","Standing with assistance","Walking with assistance","Standing alone","Walking alone"))
HR_crp$Xvar <- factor(ifelse(grepl("t2", HR_crp$X), "Ln CRP Year 1", "Ln CRP Year 2"), levels=c("Ln CRP Year 1", "Ln CRP Year 2"))

p1 <- ggplot(HR_crp, aes(x=Yvar, y=point.HR, color="WHO motor milestones", shape="Year 1")) + 
  geom_pointrange(aes(ymin=lb.HR, ymax=ub.HR),
                  size = 1) +
  scale_x_discrete(limits = rev(levels(HR_crp$Yvar)))+
  facet_grid(~Xvar, drop = T, scales = "free") +
  coord_flip() + 
  labs(y = "Adjusted hazard ratio of motor milestone achievement\nbetween 25th and 75th percentile of Th1/Th2 ratio", 
       x=element_blank()) +
  geom_hline(yintercept = 1, linetype="dashed") +
  scale_shape_manual(values=c(21,16)) +
  scale_colour_manual(values=tableau10) + 
  theme_ki() +
  theme(plot.title = element_text(hjust = 0),
        panel.spacing = unit(0, "lines"))
p1 

ggpubr::ggarrange(p1, p, widths = c(1.5, 2.5)) %>%  ggsave(filename = here("figures/gam_crp.jpg"), width = 16, height=7)



d2 <- H1_quartiles %>% rbind(H2_quartiles) %>% filter(grepl("agp", X))
d2 <- d2 %>% mutate(group=factor(ifelse(Y=="sum_who", "WHO motor milestones", ifelse(grepl("cdi", Y), "CDI Z-scores", "EASQ Z-scores")),
                                 levels=c("WHO motor milestones", "CDI Z-scores", "EASQ Z-scores")), 
                    Xvar = factor(ifelse(grepl("t2", X), "Ln AGP Year 1", "Ln AGP Year 2"), levels=c("Ln AGP Year 1", "Ln AGP Year 2")),
                    Yvar = case_when(grepl("who", Y) ~ "Sum of WHO motor milestones",
                                     grepl("und", Y) ~ "CDI comprehension Z-score",
                                     grepl("say", Y) ~ "CDI expressive language Z-score",
                                     grepl("comm", Y) ~ "EASQ communication Z-score",
                                     grepl("motor", Y) ~ "EASQ gross motor Z-score",
                                     grepl("personal", Y) ~ "EASQ personal social Z-score",
                                     grepl("combined", Y) ~ "EASQ combined Z-score"))
d2$Yvar <- factor(d2$Yvar, levels=c("Sum of WHO motor milestones", paste(c("CDI comprehension", "CDI expressive language",
                                                                           "EASQ communication", "EASQ gross motor", "EASQ personal social", "EASQ combined"), "Z-score")))
d2$Ytime <- factor(ifelse(grepl("t3|easq", d2$Y), "Year 2", "Year 1"))

p <- ggplot(d2, aes(x=Yvar, y=point.diff)) + 
  geom_pointrange(aes(ymin=lb.diff , ymax=ub.diff, color=group, group=Yvar, shape=Ytime),
                  position = position_dodge2(reverse = TRUE, width = .5),
                  size = 1) +
  scale_x_discrete(limits = rev(levels(d1$Yvar)))+
  #geom_text(aes(label=ref), position = position_nudge(y = (abs(yrange[1])+abs(yrange[2]))/10)) +
  facet_grid(~Xvar, drop = T) +
  coord_flip() + 
  labs(y = "Adjusted difference in mean development outcome\nbetween 25th and 75th percentile of Th1/Th2 ratio", 
       x=element_blank()) +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_shape_manual(values=c(21,16)) +
  scale_colour_manual(values=tableau10) + 
  theme_ki() +
  theme(plot.title = element_text(hjust = 0),
        panel.spacing = unit(0, "lines"),
        legend.position = "right")+
  guides(color = guide_legend(title="Development measurement"), shape=guide_legend(title="Development\nmeasurement time", reverse = T))
p


HR_agp <- HR %>% filter(grepl("agp", X) & !grepl("sit", Y))
HR_agp$Yvar <- factor(c("Hands-and-knees crawling","Standing with assistance",
                        "Walking with assistance","Standing alone","Walking alone"),
                      levels=c("Hands-and-knees crawling","Standing with assistance","Walking with assistance","Standing alone","Walking alone"))
HR_agp$Xvar <- factor(ifelse(grepl("t2", HR_agp$X), "Ln AGP Year 1", "Ln AGP Year 2"), levels=c("Ln AGP Year 1", "Ln AGP Year 2"))

p1 <- ggplot(HR_agp, aes(x=Yvar, y=point.HR, color="WHO motor milestones", shape="Year 1")) + 
  geom_pointrange(aes(ymin=lb.HR, ymax=ub.HR),
                  size = 1) +
  scale_x_discrete(limits = rev(levels(HR_agp$Yvar)))+
  facet_grid(~Xvar, drop = T, scales = "free") +
  coord_flip() + 
  labs(y = "Adjusted hazard ratio of motor milestone achievement\nbetween 25th and 75th percentile of Th1/Th2 ratio", 
       x=element_blank()) +
  geom_hline(yintercept = 1, linetype="dashed") +
  scale_shape_manual(values=c(21,16)) +
  scale_colour_manual(values=tableau10) + 
  theme_ki() +
  theme(plot.title = element_text(hjust = 0),
        panel.spacing = unit(0, "lines"))
p1 

ggpubr::ggarrange(p1, p, widths = c(1.5, 2.5)) %>%  ggsave(filename = here("figures/gam_agp.jpg"), width = 16, height=7)

