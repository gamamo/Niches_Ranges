# R version 4.5.1 (2025-06-13 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 26200)
# Rstudio 2025.9.1.401 (Cucumberleaf Sunflower)
# 

# Load packages -----------------------------------------------------------
library(here)          # CRAN v1.0.2
library(tidyverse)     # CRAN v2.0.0
library(fs)            # CRAN v1.6.6
library(terra)         # CRAN v1.8-70
library(tidyterra)     # CRAN v0.7.2
library(rnaturalearth) # CRAN v1.1.0
library(rphylopic)     # [::palaeoverse/rphylopic] v1.5.0.9000
library(patchwork)     # CRAN v1.3.2
library(scales)        # CRAN v1.4.0
library(lme4)          # CRAN v1.1-37
library(MuMIn)         # CRAN v1.48.11
library(broom.mixed)   # CRAN v0.2.9.6
library(janitor)       # CRAN v2.2.1
library(emmeans)       # CRAN v1.11.2-8
library(segmented)     # CRAN v2.1-4
library(merTools)      # CRAN v0.6.3
library(smatr)         # CRAN v3.4-8
library(ggridges)      # CRAN v0.5.7
library(rstatix)       # CRAN v0.7.2
library(forcats)       # CRAN v1.0.1

# Load dataset --------------------------------------
d<- read_csv("data.csv")

# Analysis -----------------------------------------------------------
### Graphs by longitudinal band ------------------------------------------

pal <- c("#56b4e9","#d55e00", "#009e73","#cc79a7")


# latitude vs range sizes
d |>  dplyr::filter(convex !=0 ) |> 
  dplyr::group_by(cLat,higher_plant_group,band) |> 
  dplyr::summarise(maxConvex = max(area),mConvex=mean(area)) |> 
  dplyr::filter(maxConvex !=0) |> 
  
  ggplot(aes(x=cLat,y=maxConvex))+
  geom_col(aes(fill = higher_plant_group))+
  geom_line(aes(x=cLat,y=mConvex))+
  facet_grid(band~higher_plant_group)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        panel.grid = element_blank())+
  scale_x_continuous(limits = c(-60,90),breaks = c(-60,-40,-20,0,20,40,60,80))+
  #scale_y_log10(labels=trans_format("log10", math_format(10^.x)))+
  scale_y_continuous(labels=comma,
                     )+
  scale_fill_manual(values=pal)+
  #labs(x="Latitude",y=expression(log["10"]~"[(Range size)"~km^2~"]"))+
  labs(x="Latitude",y=expression("Range size"~"(km2)"))+
  coord_flip() -> lrs
lrs

ggsave(here("products","RESU_general_latitude_RANGES_nolog.jpeg"),
       units="cm",dpi = 300,height = 25,width = 40)


# make a composition using range size map and latitudinal/longitudinal data
g5/lrs +
  plot_annotation(tag_levels = "a",tag_suffix = ")")
ggsave(here("products","RESU_general_latitude_RANGES_wMAP.jpeg"),
       units="cm",dpi = 300,height = 40,width = 35)

# latitude vs niche breadth
d |>  filter(convex !=0 ) |> 
  group_by(cLat,higher_plant_group,band) |> 
  summarise(maxConvex = max(convex),mConvex=mean(convex)) |> 
  filter(maxConvex !=0) |> 
  ggplot(aes(x=cLat,y=maxConvex))+
  geom_col(aes(fill = higher_plant_group))+
  geom_line(aes(x=cLat,y=mConvex))+
  facet_grid(band~higher_plant_group)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        panel.grid = element_blank())+
  scale_y_continuous(limits = c(0,80),breaks = c(5,25,50,75),expand=c(0,0))+
  scale_x_continuous(limits = c(-60,90),breaks = c(-60,-40,-20,0,20,40,60,80))+
  scale_fill_manual(values=pal)+
  labs(x="Latitude",y="Niche breadth")+
  coord_flip()
ggsave(here("products","RESU_general_latitude.jpeg"),
       units="cm",dpi = 300,height = 25,width = 30)

## 0) model 0: general patterns ---------------------------------------------------------------------
pal <- c('grey25','gold','red','steelblue1')

b <- d |> dplyr::filter(higher_plant_group == "Bryophytes") |> 
  filter(area>100) |>
  dplyr::select(convex2) |> 
  pull()
g <- d |> dplyr::filter(higher_plant_group == "Gymnosperms") |> 
  filter(area>100) |>
  dplyr::select(convex2)|> 
  pull()
f <- d |> dplyr::filter(higher_plant_group == "Ferns and lycophytes") |> 
  filter(area>100) |>
  dplyr::select(convex2)|> 
  pull()
a <- d |> dplyr::filter(higher_plant_group == "Flowering plants") |> 
  filter(area>100) |>
  dplyr::select(convex2)|> 
  pull()
ks.test(b,g) |> tidy() |> mutate(group = "Bryophytes - Gymnosperms") -> a1
ks.test(b,f) |> tidy() |> mutate(group = "Bryophytes - Ferns and lycophytes") -> a2
ks.test(b,a) |> tidy() |> mutate(group = "Bryophytes - Flowering plants") -> a3
ks.test(g,f) |> tidy() |> mutate(group = "Gymnosperms - Ferns and lycophytes") -> a4
ks.test(g,a) |> tidy() |> mutate(group = "Gymnosperms - Flowering plants") -> a5
ks.test(f,a) |> tidy() |> mutate(group = "Ferns and lycophytes - Flowering plants") -> a6

kshpg <- rbind(a1,a2,a3,a4,a5,a6) |> select(statistic, p.value, group) |> 
  dplyr::mutate(p.value = round(p.value,4)) |>
  dplyr::mutate(test = "KS") |> 
  dplyr::relocate(group,.before = statistic) |> 
  dplyr::relocate(test,.before = statistic) 

dunnhpg <- d  |>  filter(area>100) |>
  dunn_test(convex2~higher_plant_group) |> 
  unite(group, c("group1", "group2")) |> 
  dplyr::mutate(p.value = round(p.adj,4)) |> 
  dplyr::mutate(test = "KW") |>
  dplyr::select(test, statistic,p.value, group) |> 
  dplyr::relocate(group,.before = statistic) |> 
  dplyr::relocate(test,.before = statistic) 

x1 <- rbind(kshpg,dunnhpg) |> 
rename(Niche_statistic = statistic, Niche_p.value = p.value)

dars <- d |> 
  filter(area>100) |> 
  ggplot(aes(x=area10,y=hpg))+
  geom_density_ridges(aes(scale=5,color=paste(br,fl,fp,gy)),fill="transparent",linewidth=1,alpha=0.2)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        rect = element_rect(fill = "transparent"),
        plot.background = element_rect(fill="transparent"),
        panel.background = element_rect(fill="transparent"),
        panel.grid = element_blank(),
        legend.text = element_text(size=12))+
  scale_color_manual(values=pal,expand=c(0,0),
                     label=c("Bryophyes","Ferns and lycophytes","Flowering plants","Gymnosperms"),
                     guide="legend")+
  labs(x=expression(log["10"]~"[Range size"~(km^2)~"]"))+
  scale_x_continuous(limits=c(0,10),breaks = c(0,2,4,6,8,10))+
  scale_y_discrete(expand=c(0,0))
dars

b <- d |> dplyr::filter(higher_plant_group == "Bryophytes") |> 
  filter(area>100) |> 
  dplyr::select(area10) |> 
  pull()
g <- d |> dplyr::filter(higher_plant_group == "Gymnosperms") |> 
  filter(area>100) |> 
  dplyr::select(area10)|> 
  pull()
f <- d |> dplyr::filter(higher_plant_group == "Ferns and lycophytes") |>
  filter(area>100) |> 
  dplyr::select(area10)|> 
  pull()
a <- d |> dplyr::filter(higher_plant_group == "Flowering plants") |> 
  filter(area>100) |> 
  dplyr::select(area10)|> 
  pull() 
ks.test(b,g) |> tidy() |> mutate(group = "Bryophytes - Gymnosperms") -> a1
ks.test(b,f) |> tidy() |> mutate(group = "Bryophytes - Ferns and lycophytes") -> a2
ks.test(b,a) |> tidy() |> mutate(group = "Bryophytes - Flowering plants") -> a3
ks.test(g,f) |> tidy() |> mutate(group = "Gymnosperms - Ferns and lycophytes") -> a4
ks.test(g,a) |> tidy() |> mutate(group = "Gymnosperms - Flowering plants") -> a5
ks.test(f,a) |> tidy() |> mutate(group = "Ferns and lycophytes - Flowering plants") -> a6

kshpg <- rbind(a1,a2,a3,a4,a5,a6) |> select(statistic, p.value, group) |> 
  dplyr::mutate(p.value = round(p.value,4)) |>
  dplyr::mutate(test = "KS") |> 
  dplyr::relocate(group,.before = statistic) |> 
  dplyr::relocate(test,.before = statistic) 

dunnhpg <- d  |>  filter(area>100) |>
  dunn_test(area10~higher_plant_group) |> 
  unite(group, c("group1", "group2")) |> 
  dplyr::mutate(p.value = round(p.adj,4)) |> 
  dplyr::mutate(test = "KW") |>
  dplyr::select(test, statistic,p.value, group) |> 
  dplyr::relocate(group,.before = statistic) |> 
  dplyr::relocate(test,.before = statistic) 

x2<- rbind(kshpg,dunnhpg) |> 
  rename(Range_statistic = statistic, Range_p.value = p.value)

left_join(x1,x2) |> 
  flextable::flextable() |>                     
  flextable::autofit() |>     
  flextable::save_as_docx(path = here("products","Table_diff_plantgroups.docx"))

danb <- d |> 
  filter(area>100) |> 
  ggplot(aes(x=convex2,y=hpg))+
  geom_density_ridges(aes(scale=5,color=paste(br,fl,fp,gy)),fill="transparent",linewidth=1,alpha=0.2)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        rect = element_rect(fill = "transparent"),
        plot.background = element_rect(fill="transparent"),
        panel.background = element_rect(fill="transparent"),
        panel.grid = element_blank(),
        legend.text = element_text(size=12))+
  scale_color_manual(values=pal,expand=c(0,0),
                     label=c("Bryophyes","Ferns and lycophytes","Flowering plants","Gymnosperms"),
                     guide="legend")+
  labs(x="sqrt (Niche breadth)")+
  scale_x_continuous(limits=c(0,10),expand=c(0,0.1),breaks = c(0,2,4,6,8,10))+
  scale_y_discrete(expand=c(0,0))
danb


### make lat plots
aa0<- d |>  ggplot(aes(y=convex,x=midY84))+
  geom_point(color="gray70",alpha = 0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = 'none',
        rect = element_rect(fill = "transparent"),
        plot.background = element_rect(fill="transparent"),
        panel.background = element_rect(fill="transparent"),
        strip.background = element_rect(colour="white", fill="white"),
        panel.grid = element_blank(),
        legend.text = element_text(size=12))+
  scale_color_manual(values=pal)+
  labs(y="Niche breadth",x="Latitude")
aa0

pal <- c("#56b4e9","#d55e00", "#009e73","#cc79a7")
aa <- d |>  ggplot(aes(y=convex,x=midY84,color=higher_plant_group))+
  geom_point(alpha = 0.6)+
  theme_bw()+
  facet_grid(band~higher_plant_group)+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = 'none',
        rect = element_rect(fill = "transparent"),
        plot.background = element_rect(fill="transparent"),
        panel.background = element_rect(fill="transparent"),
        strip.background = element_rect(colour="white", fill="white"),
        panel.grid = element_blank(),
        legend.text = element_text(size=12))+
  scale_color_manual(values=pal)+
  labs(y="Niche breadth",x="Latitude")
aa

bb0<- d |>  ggplot(aes(y=area,x=midY84))+
  geom_point(color="gray70",alpha = 0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = 'none',
        rect = element_rect(fill = "transparent"),
        plot.background = element_rect(fill="transparent"),
        panel.background = element_rect(fill="transparent"),
        strip.background = element_rect(colour="white", fill="white"),
        panel.grid = element_blank(),
        legend.text = element_text(size=12))+
  scale_color_manual(values=pal)+
  scale_y_continuous(labels=comma)+ 
  labs(y=expression("Range size"~km^2),x="Latitude")
bb0

bb <- d |>  ggplot(aes(y=area,x=midY84,color=higher_plant_group))+
  geom_point(alpha = 0.6)+
  facet_grid(band~higher_plant_group,)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = 'none',
        rect = element_rect(fill = "transparent"),
        plot.background = element_rect(fill="transparent"),
        panel.background = element_rect(fill="transparent"),
        strip.background = element_rect(colour="white", fill="white"),
        panel.grid = element_blank(),
        legend.text = element_text(size=12))+
  scale_y_continuous(labels=comma)+  
  scale_color_manual(values=pal)+
  labs(y=expression("Range size"~km^2),x="Latitude")
bb

##### make composite figure -----------------------------------------------------------
dars + danb+
  plot_annotation(tag_levels = "a",tag_suffix = ")")+
  plot_layout(guides = 'collect',ncol=2)&
  theme(legend.position="bottom")
ggsave(here("products","RESU_dist_groups.jpeg"),
       units="cm",dpi = 300,height = 15,width = 30)


cp=aa0 / aa+
  plot_annotation(tag_levels = "a",tag_suffix = ")")+
  plot_layout(guides = 'collect',heights = c(1,2))
ggsave(here("products","RESU_latGrad_Niche_groups.jpeg"),
       units="cm",dpi = 300,height = 30,width = 30,plot=cp)


cp=bb0 / bb+
  plot_annotation(tag_levels = "a",tag_suffix = ")")+
  plot_layout(guides = 'collect',heights = c(1,2))
ggsave(here("products","RESU_latGrad_Range_groups.jpeg"),
       units="cm",dpi = 300,height = 30,width = 30,plot=cp)



## 1) model 1: range size per region -----------------------------------------------
pal2 <- c("#56b4e9","#d55e00","#009e73")

### All species ---------------------------------------------------------
#### NICHE  -------------------------------------------------------------
t <- d |> dplyr::filter(region == "Tropical") |> 
  dplyr::select(convex) |> 
  pull()
e <- d |> dplyr::filter(region == "Extratropical") |> 
  dplyr::select(convex)|> 
  pull()
b <- d |> dplyr::filter(region == "Latitudinal generalists") |> 
  dplyr::select(convex)|> 
  pull()
ks.test(t,e)
ks.test(t,b)
ks.test(e,b)
d |> dunn_test(convex~region)

nd1 <- d |> 
  dplyr::filter(convex>1)|> 
  ggplot(aes(x=convex2,y=region,fill=region))+
  #ggplot(aes(x=area10,y=r))+
  #geom_density_ridges(aes(scale=5,color=paste(tr,et,bo)),fill="transparent",linewidth=1,alpha=0.2)+
 geom_density_ridges2(aes(fill=region),scale=2,color="white")+
 stat_density_ridges(quantile_lines = TRUE,
                     quantiles=2,color="white",scale=2)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        panel.grid = element_blank())+
  scale_fill_manual(values=pal2)+
  scale_y_discrete(expand=c(0,0.1),labels=c("Latitudinal \n generalists",
                                            "Extratropical", "Tropical"))+
  scale_x_continuous(limits=c(0,10), breaks=c(0,2,4,6,8,10))+
  labs(x="sqrt (Niche breadth)")+
  annotate("segment",x=8.6,xend=8.6, y=2,yend = 3)+
  annotate("segment",x=9.3,xend=9.3, y=2,yend = 1)+
  annotate("segment",x=10,xend=10, y=3,yend = 1)+
  
  annotate("segment",x=8.4,xend=8.6, y=3,yend = 3)+
  annotate("segment",x=8.4,xend=8.6, y=2,yend = 2)+
  
  annotate("segment",x=9.1,xend=9.3, y=2,yend = 2)+
  annotate("segment",x=9.1,xend=9.3, y=1,yend = 1)+
  
  annotate("segment",x=9.8,xend=10, y=3,yend = 3)+
  annotate("segment",x=9.8,xend=10, y=1,yend = 1)+
  
  annotate("text",label="ks ***", x=7.8, y=2.8,size=5)+
  annotate("text",label="ks ***", x=8.6, y=1.8,size=5)+
  annotate("text",label="ks ***", x=9.3, y=2.8  ,size=5)+
  
  annotate("text",label="kw ***", x=7.8, y=2.5,size=5)+
  annotate("text",label="kw ***", x=8.6, y=1.5,size=5)+
  annotate("text",label="kw ***", x=9.3, y=2.5  ,size=5)
nd1

# inset graph
pal <- c("#56b4e9","#d55e00", "#009e73")

ind1 <- d |>  filter(convex !=0 ) |> 
  group_by(cLat,region) |> 
  summarise(maxConvex = max(convex2),mConvex=median(convex2)) |> 
  filter(maxConvex !=0) |> 
  ggplot(aes(x=cLat,y=mConvex,fill=region))+
  geom_col(alpha=0.6)+
  #geom_line(aes(x=cLat,y=mConvex,color=region),linewidth=1)+
  facet_wrap(~fct_rev(region))+
  theme_bw()+
  theme(text=element_text(size=12,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        rect = element_rect(fill = "transparent"),
        axis.title.x = element_blank())+
  scale_x_continuous(limits = c(-60,90),breaks = c(-60,-40,-20,0,20,40,60,80),expand = c(0,0))+
  scale_y_continuous(limits=c(0,8),breaks = c(0,4,8))+
  #scale_y_log10(labels=trans_format("log10", math_format(10^.x)))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(x="Latitude",y="Niche breadth")+
  coord_flip()
ind1


#### RANGE --------------------------------------------------------------------
t <- d |> dplyr::filter(region == "Tropical") |> 
  dplyr::filter(area>100) |> 
  dplyr::select(area10) |> 
  pull()
e <- d |> dplyr::filter(region == "Extratropical") |> 
  dplyr::filter(area>100) |> 
  dplyr::select(area10) |> 
  pull()
b <- d |> dplyr::filter(region == "Latitudinal generalists") |> 
  dplyr::filter(area>100) |> 
  dplyr::select(area10) |> 
  pull()
ks.test(t,e)
ks.test(t,b)
ks.test(e,b)
d |> 
  dplyr::filter(area>100) |> 
  dunn_test(area10~region)

ad1 <- d |> 
  dplyr::filter(area>100) |> 
  ggplot(aes(x=area10,y=region,fill=region))+
  geom_density_ridges2(aes(fill=region),scale=2,color="white")+
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles=2,color="white",scale=2)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())+
  scale_fill_manual(values=pal2,expand=c(0,0))+
  scale_y_discrete(expand=c(0,0.1),labels=c("Latitudinal \n generalists",
                                            "Extratropical", "Tropical"))+
  #scale_y_continuous(limits=c(1,8))+
  labs(x=expression(log["10"]~ "[Range size"~(km^2)~"]"))  +
  scale_x_continuous(limits=c(2,10))+
  annotate("segment",x=8.6,xend=8.6, y=2,yend = 3)+
  annotate("segment",x=9.3,xend=9.3, y=2,yend = 1)+
  annotate("segment",x=10,xend=10, y=3,yend = 1)+
  
  annotate("segment",x=8.4,xend=8.6, y=3,yend = 3)+
  annotate("segment",x=8.4,xend=8.6, y=2,yend = 2)+
  
  annotate("segment",x=9.1,xend=9.3, y=2,yend = 2)+
  annotate("segment",x=9.1,xend=9.3, y=1,yend = 1)+
  
  annotate("segment",x=9.8,xend=10, y=3,yend = 3)+
  annotate("segment",x=9.8,xend=10, y=1,yend = 1)+
  
  annotate("text",label="ks ***", x=7.8, y=2.8,size=5)+
  annotate("text",label="ks ***", x=8.6, y=1.8,size=5)+
  annotate("text",label="ks ***", x=9.3, y=2.8  ,size=5)+
  
  annotate("text",label="kw ***", x=7.8, y=2.5,size=5)+
  annotate("text",label="kw ***", x=8.6, y=1.5,size=5)+
  annotate("text",label="kw ***", x=9.3, y=2.5  ,size=5)
ad1

# inset graph
pal <- c("#56b4e9","#d55e00", "#009e73")

iad1 <- d |>   
  dplyr::filter(area>100) |>  
  group_by(cLat,region) |> 
  summarise(maxConvex = max(area10),mConvex=median(area10)) |> 
  ggplot(aes(x=cLat,y=mConvex,fill=region))+
  geom_col(alpha=0.6)+
  #geom_line(aes(x=cLat,y=mConvex,color=region),linewidth=1)+
  facet_wrap(~fct_rev(region))+
  theme_bw()+
  theme(text=element_text(size=12,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        rect = element_rect(fill = "transparent"),
        axis.title.x = element_blank())+
  scale_x_continuous(limits = c(-60,90),breaks = c(-60,-40,-20,0,20,40,60,80),
                     expand = c(0,0))+
  scale_y_continuous(limits = c(0,8),breaks = c(0,4,8))+
  #scale_y_log10(labels=trans_format("log10", math_format(10^.x)))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(x="Latitude",y=expression(log["10"]~ "[Range size"~(km^2)~"]"))+
  coord_flip()
iad1

### modelling -------------------------------------------------------------------------
#### SMA -----------------------------------------------------------------------

dsma <- d |> select(convex2,region,area10) |> drop_na() |> as.data.frame()
s1 <- sma(area10 ~ convex2+region,data=dsma, multcomp=TRUE,type="elevation")
cc <- as.data.frame.table(multcompmatrix(s1)) |> mutate(test="elevation")
s2 <- sma(area10 ~ convex2+region,data=dsma, multcomp=TRUE,type="shift")
multcompmatrix(s2)
cc2 <- as.data.frame.table(multcompmatrix(s2)) |> mutate(test="shift")
s3 <- sma(area10 ~ convex2*region,data=dsma, multcomp=TRUE)
multcompmatrix(s3)
cc3 <- as.data.frame.table(multcompmatrix(s2)) |> mutate(test="slope")

resu <- rbind(cc,cc2,cc3)
resu <- resu[resu$Freq!= "-",]
resu <- resu |> unite(comparison, c("region_1","region_2"),sep = "-") |> 
  flextable::flextable() %>%                     
  flextable::autofit() %>%    
  flextable::save_as_docx(path = here("products","Table_SMAmodel1.docx"))

#### GLMM --------------------------------------------------------------------
d1 <- d 

m1 <- lmer(data=d1, area10 ~ convex2:region +(1|region)) #random slope and intercept
r.squaredGLMM(m1)
m1coef <- tidy(m1) |> print()

d1 <- d1 |> mutate(fit.m = predict(m1, re.form=NA),
                   fit.c = predict(m1, re.form=NULL))

pairs(emtrends(m1,"region",var="convex2"),pbkrtest.limit = 156135) |>
  tidy()

m1g <-d1 |>
  ggplot(aes(y=area10,x=convex2,color=region))+
  geom_point(color="gray70",alpha = 0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())+
  geom_line(aes(y=fit.c,col=region),linewidth=1)+
  scale_color_manual(values =c("#56b4e9","#d55e00","#009e73"))+
  labs(y=expression("log10 [Range size"~(km^2)~"]"),x="sqrt (Niche breadth)")+
  scale_y_continuous(limits=c(1,8))+
  scale_x_continuous(limits=c(1,8), breaks=c(0,2,4,6,8))+
  annotate("text",label = "conditional~R^2==0.49",
           parse = TRUE,
           #expression(conditional~R^2 ~ "= 0.49"),
           x=1,y=8,size=5,hjust = 0)+

  #legend
  annotate("segment", y = 2.5, yend = 2.5, x = 5.5, xend = 6.1, 
           color = "#56b4e9",linewidth = 1)+
  annotate("segment", y = 2, yend = 2, x = 5.5, xend = 6.1, 
           color = "#d55e00",linewidth = 1)+
  annotate("segment", y = 1.5, yend = 1.5, x = 5.5, xend = 6.1, 
           color = "#009e73",linewidth = 1)+
  
  annotate("text", label="slope = 0.39", x = 7.2, y =2.5,size= 5)+
  annotate("text", label="slope = 0.65", x = 7.2, y =2,size=5)+
  annotate("text", label="slope = 0.69", x = 7.2, y =1.5,size=5)+
  annotate("text", label= "***", x = 7.2, y =2.8,size=5)
m1g


##### make composite figures -------------------------------------------------------------

cp=  (nd1 + 
    inset_element(ind1, left = 0.6, bottom = 0.6, right = 0.98, top = 0.93, align_to = 'full'))+
(ad1+inset_element(iad1, left = 0.6, bottom = 0.6, right = 0.98, top = 0.93, align_to = 'full'))+
  m1g+
  plot_layout(guides = 'collect',ncol = 3,axis_titles = "collect",tag_level = "new")+
  plot_annotation(tag_levels = list(c("a)","","b)","","c)")))
ggsave(here("products","RESU_dist_region_allPlants_v3.jpeg"),
       units="cm",dpi = 300,height = 15,width = 40,plot=cp)
ggsave(here("products","RESU_dist_region_allPlants_v3.tiff"),
       units="cm",dpi = 300,height = 15,width = 40,plot=cp)



(g5 + g0) / 
  (m1g | (l1+n1))+
  plot_annotation(tag_levels = "a", tag_suffix = ")")
  ggsave(here("products","RESU_compositionFigure1_v1.jpeg"),
         units="cm",dpi = 300,height = 30,width = 30)



### Bryophytes ------------------------------------------------------------------
db <- d |> 
  dplyr::filter(higher_plant_group=="Bryophytes")

#### distribution graph and analysis --------------------------------------------
##### NICHE ---------------------------------------------------------------------
t <- db |> dplyr::filter(region == "Tropical") |> 
  dplyr::select(convex) |> 
  pull()
e <- db |> dplyr::filter(region == "Extratropical") |> 
  dplyr::select(convex)|> 
  pull()
b <- db |> dplyr::filter(region == "Latitudinal generalists") |> 
  dplyr::select(convex)|> 
  pull()
ks.test(t,e)
ks.test(t,b)
ks.test(e,b)
db |> dunn_test(convex~region)

ndb <- db |> 
ggplot(aes(x=convex2,y=region,fill=region))+
  geom_density_ridges2(aes(fill=region),scale=2,color="white")+
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles=2,color="white",scale=2)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        panel.grid = element_blank())+
  scale_fill_manual(values=pal2,expand=c(0,0))+
  scale_y_discrete(expand=c(0,0.1),labels=c("Latitudinal \n generalists",
                                            "Extratropical", "Tropical"))+
  scale_x_continuous(limits=c(0,10), breaks=c(0,2,4,6,8,10))+
  labs(x="sqrt (Niche breadth)",title="Bryophytes")+

  annotate("segment",x=8.6,xend=8.6, y=2,yend = 3)+
  annotate("segment",x=9.3,xend=9.3, y=2,yend = 1)+
  annotate("segment",x=10,xend=10, y=3,yend = 1)+
  
  annotate("segment",x=8.4,xend=8.6, y=3,yend = 3)+
  annotate("segment",x=8.4,xend=8.6, y=2,yend = 2)+
  
  annotate("segment",x=9.1,xend=9.3, y=2,yend = 2)+
  annotate("segment",x=9.1,xend=9.3, y=1,yend = 1)+
  
  annotate("segment",x=9.8,xend=10, y=3,yend = 3)+
  annotate("segment",x=9.8,xend=10, y=1,yend = 1)+
  
  annotate("text",label="ks ***", x=7.8, y=2.8,size=5)+
  annotate("text",label="ks ***", x=8.6, y=1.8,size=5)+
  annotate("text",label="ks ***", x=9.3, y=2.8  ,size=5)+
  
  annotate("text",label="kw ***", x=7.8, y=2.5,size=5)+
  annotate("text",label="kw ***", x=8.6, y=1.5,size=5)+
  annotate("text",label="kw ***", x=9.3, y=2.5  ,size=5)
ndb


##### RANGE ---------------------------------------------------------------------
t <- db |> dplyr::filter(region == "Tropical") |> 
  dplyr::filter(area>100) |> 
  dplyr::select(area10) |> 
  pull()
e <- db |> dplyr::filter(region == "Extratropical") |> 
  dplyr::filter(area>100) |> 
  dplyr::select(area10) |> 
  pull()
b <- db |> dplyr::filter(region == "Latitudinal generalists") |> 
  dplyr::filter(area>100) |> 
  dplyr::select(area10) |> 
  pull()
ks.test(t,e)
ks.test(t,b)
ks.test(e,b)
db |> 
  dplyr::filter(area>100) |> 
  dunn_test(area10~region)

adb <- db |> 
  dplyr::filter(area>100) |> 
  ggplot(aes(x=area10,y=region,fill=region))+
  geom_density_ridges2(aes(fill=region),scale=2,color="white")+
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles=2,color="white",scale=2)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        panel.grid = element_blank())+
  scale_fill_manual(values=pal2,expand=c(0,0))+
  scale_y_discrete(expand=c(0,0.1),labels=c("Latitudinal \n generalists",
                                            "Extratropical", "Tropical"))+
  labs(x=expression(log["10"]~ "[Range size"~(km^2)~"]"),title="Bryophytes")  +
  scale_x_continuous(limits=c(2,10))+
  annotate("segment",x=8.6,xend=8.6, y=2,yend = 3)+
  annotate("segment",x=9.3,xend=9.3, y=2,yend = 1)+
  annotate("segment",x=10,xend=10, y=3,yend = 1)+
  
  annotate("segment",x=8.4,xend=8.6, y=3,yend = 3)+
  annotate("segment",x=8.4,xend=8.6, y=2,yend = 2)+
  
  annotate("segment",x=9.1,xend=9.3, y=2,yend = 2)+
  annotate("segment",x=9.1,xend=9.3, y=1,yend = 1)+
  
  annotate("segment",x=9.8,xend=10, y=3,yend = 3)+
  annotate("segment",x=9.8,xend=10, y=1,yend = 1)+
  
  annotate("text",label="ks ***", x=7.8, y=2.8,size=5)+
  annotate("text",label="ks ***", x=8.6, y=1.8,size=5)+
  annotate("text",label="ks ***", x=9.3, y=2.8  ,size=5)+
  
  annotate("text",label="kw ***", x=7.8, y=2.5,size=5)+
  annotate("text",label="kw ***", x=8.6, y=1.5,size=5)+
  annotate("text",label="kw ***", x=9.3, y=2.5  ,size=5)
adb

##### modelling ----------------------------------------------------------------------
#SMA
dsma <- db |> select(convex2,region,area10) |> drop_na() |> as.data.frame()
s1 <- sma(area10 ~ convex2+region,data=dsma, multcomp=TRUE,type="elevation")
multcompmatrix(s1)
summary(s1)
s2 <- sma(area10 ~ convex2+region,data=dsma, multcomp=TRUE,type="shift")
multcompmatrix(s2)
s3 <- sma(area10 ~ convex2*region,data=dsma, multcomp=TRUE)
multcompmatrix(s3)

#glmmm
m1b <- lmer(data=db, area10 ~ convex2:region +(1|region)) #random slope and intercept
r.squaredGLMM(m1b)
m1bcoef <- tidy(m1b,"ran_coefs") |> print()
coef(m1b)

tb <- pairs(emtrends(m1b,"region",var="convex2"),pbkrtest.limit = 6767) |>
  tidy() |> 
  mutate(group = "Bryophytes")

db <- db |> mutate(fit.m = predict(m1b, re.form=NA),
                   fit.c = predict(m1b, re.form=NULL))

bryopic <- get_uuid(name="Andreaeaceae", n = 1)
m1bg <-db |> 
  ggplot(aes(y=area10,x=convex2,color=region))+
  geom_point(color="gray70",alpha = 0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank())+
  geom_line(aes(y=fit.c,col=region),linewidth=1)+
  scale_color_manual(values =c("#56b4e9","#d55e00","#009e73"))+
  #geom_abline(intercept = 3.68, slope = 0.513,color="#56b4e9", linewidth=1)+ #Both
  #geom_abline(intercept = 2.70, slope = 0.745,color="#d55e00", linewidth=1)+ #Extratropical
  #geom_abline(intercept = 2.78, slope = 0.564,color="#009e73", linewidth=1)+ #Tropical
  labs(y=expression("log10 [Range size"~(km^2)~"]"),x="sqrt (Niche breadth)",title="Bryophytes")+
  scale_y_continuous(limits=c(1,8))+
  scale_x_continuous(limits=c(1,8))+
  #annotate("text",label =  "conditional~R^2==0.54",
             #expression(conditional~R^2 == 0.54),
  #         x=0,y=8,size=5,hjust = 0)+
  #legend
  annotate("segment", y = 2.5, yend = 2.5, x = 5.5, xend = 6.1, 
           color = "#56b4e9",linewidth = 1)+
  annotate("segment", y = 2, yend = 2, x = 5.5, xend = 6.1, 
           color = "#d55e00",linewidth = 1)+
  annotate("segment", y = 1.5, yend = 1.5, x = 5.5, xend = 6.1, 
           color = "#009e73",linewidth = 1)+
  
  annotate("text", label="slope = 0.49", x = 7.2, y =2.5,size= 5)+
  annotate("text", label="slope = 0.71", x = 7.2, y =2,size=5)+
  annotate("text", label="slope = 0.57", x = 7.2, y =1.5,size=5)+
  annotate("text",label = "conditional~R^2==0.51",
           parse = TRUE,
            # expression(conditional~R^2 == 0.51),
           x=1,y=8,size=5,hjust = 0)+
  add_phylopic(uuid = bryopic,  x = 8, y = 4,height= 2)

m1bg

'Range size * 10,000 '('km'^2)


### Ferns ----------------------------------------------------------------

df <- d |> 
  dplyr::filter(higher_plant_group=="Ferns and lycophytes")

#### distribution graph and analysis --------------------------------------------
##### NICHE ---------------------------------------------------------------------
t <- df |> dplyr::filter(region == "Tropical") |> 
  dplyr::select(convex) |> 
  pull()
e <- df |> dplyr::filter(region == "Extratropical") |> 
  dplyr::select(convex)|> 
  pull()
b <- df |> dplyr::filter(region == "Latitudinal generalists") |> 
  dplyr::select(convex)|> 
  pull()
ks.test(t,e)
ks.test(t,b)
ks.test(e,b)
df |> dunn_test(convex~region)

ndf <- df |> 
  ggplot(aes(x=convex2,y=region,fill=region))+
  geom_density_ridges2(aes(fill=region),scale=2,color="white")+
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles=2,color="white",scale=2)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())+
  scale_fill_manual(values=pal2,expand=c(0,0))+
  scale_y_discrete(expand=c(0,0.1),labels=c("Latitudinal \n generalists",
                                            "Extratropical", "Tropical"))+
  scale_x_continuous(limits=c(0,10), breaks=c(0,2,4,6,8,10))+
  labs(x="sqrt (Niche breadth)",title="Ferns and lycophytes")  +
  
  annotate("segment",x=8.6,xend=8.6, y=2,yend = 3)+
  annotate("segment",x=9.3,xend=9.3, y=2,yend = 1)+
  annotate("segment",x=10,xend=10, y=3,yend = 1)+
  
  annotate("segment",x=8.4,xend=8.6, y=3,yend = 3)+
  annotate("segment",x=8.4,xend=8.6, y=2,yend = 2)+
  
  annotate("segment",x=9.1,xend=9.3, y=2,yend = 2)+
  annotate("segment",x=9.1,xend=9.3, y=1,yend = 1)+
  
  annotate("segment",x=9.8,xend=10, y=3,yend = 3)+
  annotate("segment",x=9.8,xend=10, y=1,yend = 1)+
  
  annotate("text",label="ks *", x=7.8, y=2.8,size=5)+
  annotate("text",label="ks ***", x=8.6, y=1.8,size=5)+
  annotate("text",label="ks ***", x=9.3, y=2.8  ,size=5)+
  
  annotate("text",label="kw ***", x=7.8, y=2.5,size=5)+
  annotate("text",label="kw ***", x=8.6, y=1.5,size=5)+
  annotate("text",label="kw ***", x=9.3, y=2.5  ,size=5)
ndf

##### RANGE --------------------------------------------------------

t <- df |> dplyr::filter(region == "Tropical") |> 
  dplyr::filter(area>100) |> 
  dplyr::select(area10) |> 
  pull()
e <- df |> dplyr::filter(region == "Extratropical") |> 
  dplyr::filter(area>100) |> 
  dplyr::select(area10) |> 
  pull()
b <- df |> dplyr::filter(region == "Latitudinal generalists") |> 
  dplyr::filter(area>100) |> 
  dplyr::select(area10) |> 
  pull()
ks.test(t,e)
ks.test(t,b)
ks.test(e,b)
df |> 
  dplyr::filter(area>100) |> 
  dunn_test(area10~region) |> 
  mutate(p2 = round(p,5))


adf <- df |> 
  dplyr::filter(area>100) |> 
  ggplot(aes(x=area10,y=region,fill=region))+
  geom_density_ridges2(aes(fill=region),scale=2,color="white")+
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles=2,color="white",scale=2)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())+
  scale_fill_manual(values=pal2,expand=c(0,0))+
  #scale_y_continuous(limits=c(1,8))+
  scale_y_discrete(expand=c(0,0.1))+
  labs(x=expression(log["10"]~ "[Range size"~(km^2)~"]"),title="Ferns and lycophytes")+  
  scale_x_continuous(limits=c(2,10))+
  annotate("segment",x=8.6,xend=8.6, y=2,yend = 3)+
  annotate("segment",x=9.3,xend=9.3, y=2,yend = 1)+
  annotate("segment",x=10,xend=10, y=3,yend = 1)+
  
  annotate("segment",x=8.4,xend=8.6, y=3,yend = 3)+
  annotate("segment",x=8.4,xend=8.6, y=2,yend = 2)+
  
  annotate("segment",x=9.1,xend=9.3, y=2,yend = 2)+
  annotate("segment",x=9.1,xend=9.3, y=1,yend = 1)+
  
  annotate("segment",x=9.8,xend=10, y=3,yend = 3)+
  annotate("segment",x=9.8,xend=10, y=1,yend = 1)+
  
  annotate("text",label="ks ***", x=7.8, y=2.8,size=5)+
  annotate("text",label="ks ***", x=8.6, y=1.8,size=5)+
  annotate("text",label="ks ***", x=9.3, y=2.8  ,size=5)+
  
  annotate("text",label="kw **", x=7.8, y=2.5,size=5)+
  annotate("text",label="kw ***", x=8.6, y=1.5,size=5)+
  annotate("text",label="kw ***", x=9.3, y=2.5  ,size=5)
adf

#### modelling -----------------------------------------------------------
#SMA
dsma <- df |> select(convex2,region,area10) |> drop_na() |> as.data.frame()
s1 <- sma(area10 ~ convex2+region,data=dsma, multcomp=TRUE,type="elevation")
multcompmatrix(s1)
s2 <- sma(area10 ~ convex2+region,data=dsma, multcomp=TRUE,type="shift")
multcompmatrix(s2)
s3 <- sma(area10 ~ convex2*region,data=dsma, multcomp=TRUE)
multcompmatrix(s3)

#GLMM
m1f <- lmer(data=df, area10 ~ convex2:region + (1|region) )
r.squaredGLMM(m1f)
m1fcoef <- tidy(m1f,"ran_coefs") |> print()

tf <- pairs(emtrends(m1f,"region",var="convex2"),pbkrtest.limit = 6769) |> 
  tidy() |> 
  mutate(group="Ferns")

df <- df |> mutate(fit.m = predict(m1f, re.form=NA),
                   fit.c = predict(m1f, re.form=NULL))

fernpic <- get_uuid(name="Blechum", n = 1)
m1fg <-df |> 
  ggplot(aes(y=area10,x=convex2,color=region))+
  geom_point(color="gray70",alpha = 0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank())+
  geom_line(aes(y=fit.c,col=region),linewidth=1)+
  scale_color_manual(values =c("#56b4e9","#d55e00","#009e73"))+
  #geom_abline(intercept = 4.01, slope = 0.461,color="#56b4e9", linewidth=1)+ #Both
  #geom_abline(intercept = 3.24, slope = 0.542,color="#d55e00", linewidth=1)+ #Extratropical
  #geom_abline(intercept = 2.53, slope = 0.760,color="#009e73", linewidth=1)+ #Tropical
  labs(y=expression("log10 [Range size"~(km^2)~"]"),x="sqrt (Niche breadth)",title="Ferns")+
  scale_y_continuous(limits=c(1,8))+
  scale_x_continuous(limits=c(1,8))+
  #annotate("text",label = 
             #expression(conditional~R^2 == 0.57),
   #        x=0,y=8,size=5,hjust = 0)+
  #legend
  annotate("segment", y = 2.5, yend = 2.5, x = 5.5, xend = 6.1, 
           color = "#56b4e9",linewidth = 1)+
  annotate("segment", y = 2, yend = 2, x = 5.5, xend = 6.1, 
           color = "#d55e00",linewidth = 1)+
  annotate("segment", y = 1.5, yend = 1.5, x = 5.5, xend = 6.1, 
           color = "#009e73",linewidth = 1)+
  
  annotate("text", label="slope = 0.44", x = 7.2, y =2.5,size= 5)+
  annotate("text", label="slope = 0.52", x = 7.2, y =2,size=5)+
  annotate("text", label="slope = 0.76", x = 7.2, y =1.5,size=5)+
  annotate("text",label = "conditional~R^2==0.56",
           parse = TRUE,
             #expression(conditional~R^2 == 0.56),
           x=1,y=8,size=5,hjust = 0)+
  add_phylopic(uuid = fernpic, x = 7.4, y = 4, ysize = 1)
m1fg

###Gymnosperms ------------------------------------------------------------------

dg <- d |> 
  dplyr::filter(higher_plant_group=="Gymnosperms")

#### distribution graph and analysis --------------------------------------------
##### NICHE ---------------------------------------------------------------------
t <- dg |> dplyr::filter(region == "Tropical") |> 
  dplyr::select(convex) |> 
  pull()
e <- dg |> dplyr::filter(region == "Extratropical") |> 
  dplyr::select(convex)|> 
  pull()
b <- dg |> dplyr::filter(region == "Both") |> 
  dplyr::select(convex)|> 
  pull()
ks.test(t,e)
ks.test(t,b)
ks.test(e,b)
dg |> dunn_test(convex~region)

ndg <- dg |> 
  ggplot(aes(x=convex2,y=region,fill=region))+
  geom_density_ridges2(aes(fill=region),scale=2,color="white")+
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles=2,color="white",scale=2)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())+
  scale_fill_manual(values=pal2,expand=c(0,0))+
  #scale_y_continuous(limits=c(1,8))+
  scale_x_continuous(limits=c(0,10),breaks=c(0,2,4,6,8,10))+
  scale_y_discrete(expand=c(0,0.1))+
  labs(x="sqrt (Niche breadth)",title="Gymnosperms")+
  
  annotate("segment",x=8.6,xend=8.6, y=2,yend = 3)+
  annotate("segment",x=9.3,xend=9.3, y=2,yend = 1)+
  annotate("segment",x=10,xend=10, y=3,yend = 1)+
  
  annotate("segment",x=8.4,xend=8.6, y=3,yend = 3)+
  annotate("segment",x=8.4,xend=8.6, y=2,yend = 2)+
  
  annotate("segment",x=9.1,xend=9.3, y=2,yend = 2)+
  annotate("segment",x=9.1,xend=9.3, y=1,yend = 1)+
  
  annotate("segment",x=9.8,xend=10, y=3,yend = 3)+
  annotate("segment",x=9.8,xend=10, y=1,yend = 1)+
  
  annotate("text",label="ks n.s.", x=7.8, y=2.8,size=5)+
  annotate("text",label="ks ***", x=8.6, y=1.8,size=5)+
  annotate("text",label="ks ***", x=9.3, y=2.8  ,size=5)+
  
  annotate("text",label="kw n.s.", x=7.8, y=2.5,size=5)+
  annotate("text",label="kw ***", x=8.6, y=1.5,size=5)+
  annotate("text",label="kw ***", x=9.3, y=2.5  ,size=5)
ndg

##### RANGE --------------------------------------------------------------------------------

t <- dg |> dplyr::filter(region == "Tropical") |> 
  dplyr::filter(area>100) |> 
  dplyr::select(area10) |> 
  pull()
e <- dg |> dplyr::filter(region == "Extratropical") |> 
  dplyr::filter(area>100) |> 
  dplyr::select(area10) |> 
  pull()
b <- dg |> dplyr::filter(region == "Latitudinal generalists") |> 
  dplyr::filter(area>100) |> 
  dplyr::select(area10) |> 
  pull()
ks.test(t,e)
ks.test(t,b)
ks.test(e,b)
dg |> 
  dplyr::filter(area>100) |> 
  dunn_test(area10~region) |> 
  mutate(p2 = round(p,5))

adg <- dg |> 
  dplyr::filter(area>100) |> 
  ggplot(aes(x=area10,y=region,fill=region))+
  geom_density_ridges2(aes(fill=region),scale=2,color="white")+
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles=2,color="white",scale=2)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())+
  scale_fill_manual(values=pal2,expand=c(0,0))+
  labs(x=expression(log["10"]~ "[Range size"~(km^2)~"]"),title="Gymnosperms")+
  scale_x_continuous(limits=c(2,10))+
  scale_y_discrete(expand=c(0,0.1))+
  annotate("segment",x=8.6,xend=8.6, y=2,yend = 3)+
  annotate("segment",x=9.3,xend=9.3, y=2,yend = 1)+
  annotate("segment",x=10,xend=10, y=3,yend = 1)+
  
  annotate("segment",x=8.4,xend=8.6, y=3,yend = 3)+
  annotate("segment",x=8.4,xend=8.6, y=2,yend = 2)+
  
  annotate("segment",x=9.1,xend=9.3, y=2,yend = 2)+
  annotate("segment",x=9.1,xend=9.3, y=1,yend = 1)+
  
  annotate("segment",x=9.8,xend=10, y=3,yend = 3)+
  annotate("segment",x=9.8,xend=10, y=1,yend = 1)+
  
  annotate("text",label="ks **", x=7.8, y=2.8,size=5)+
  annotate("text",label="ks ***", x=8.6, y=1.8,size=5)+
  annotate("text",label="ks ***", x=9.3, y=2.8  ,size=5)+
  
  annotate("text",label="kw **", x=7.8, y=2.5,size=5)+
  annotate("text",label="kw ***", x=8.6, y=1.5,size=5)+
  annotate("text",label="kw ***", x=9.3, y=2.5  ,size=5)
adg

#### modelling -----------------------------------------------------------
#SMA
dsma <- dg|> select(convex2,region,area10) |> drop_na() |> as.data.frame()
s1 <- sma(area10 ~ convex2+region,data=dsma, multcomp=TRUE,type="elevation")
multcompmatrix(s1)
s2 <- sma(area10 ~ convex2+region,data=dsma, multcomp=TRUE,type="shift")
multcompmatrix(s2)
s3 <- sma(area10 ~ convex2*region,data=dsma, multcomp=TRUE)
multcompmatrix(s3)


#glmm
m1g <- lmer(data=dg, area10 ~ convex2:region + (1|region) )
r.squaredGLMM(m1g)
m1gcoef <- tidy(m1g,"ran_coefs") |> print()

tg <- pairs(emtrends(m1g,"region",var="convex2"),pbkrtest.limit = 6769) |> 
  tidy() |> 
  mutate(group="Gymnosperms")

dg <- dg |> mutate(fit.m = predict(m1g, re.form=NA),
                   fit.c = predict(m1g, re.form=NULL))

gypic <- get_uuid(name="Sequoia sempervirens", n = 1)
m1gg <-dg |> 
  #dplyr::filter(area >100) |> 
  ggplot(aes(y=area10,x=convex2,color=region))+
  geom_point(color="gray70",alpha = 0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank())+
  geom_line(aes(y=fit.c,col=region),linewidth=1)+
  scale_color_manual(values =c("#56b4e9","#d55e00","#009e73"))+
  #geom_abline(intercept = 4.44, slope = 0.389,color="#56b4e9", linewidth=1)+ #Both
  #geom_abline(intercept = 3.52, slope = 0.625,color="#d55e00", linewidth=1)+ #Extratropical
  #geom_abline(intercept = 3.11, slope = 0.687,color="#009e73", linewidth=1)+ #Tropical
  labs(y=expression("log10 [Range size"~(km^2)~"]"),x="sqrt (Niche breadth)",title="Gymnosperms")+
  scale_y_continuous(limits=c(1,8))+
  scale_x_continuous(limits=c(1,8))+
  #annotate("text",label = expression(conditional~R^2 == 0.49),
  #         x=0,y=8,size=5,hjust = 0)+
  #legend
  #annotate("segment", y = 2.5, yend = 2.5, x = 3.5, xend = 4.3, 
  #         color = "#56b4e9",linewidth = 1)+
  #annotate("segment", y = 2, yend = 2, x = 3.5, xend = 4.3, 
  #         color = "#d55e00",linewidth = 1)+
  #annotate("segment", y = 1.5, yend = 1.5, x = 3.5, xend = 4.3, 
  #         color = "#009e73",linewidth = 1)+
  
  #annotate("text", label="Both", x = 5.2, y =2.5,size= 5)+
  #annotate("text", label="Extratropical", x = 5.2, y =2,size=5)+
  #annotate("text", label="Tropical", x = 5.2, y =1.5,size=5)
  
  #legend
  annotate("segment", y = 2.5, yend = 2.5, x = 5.5, xend = 6.1, 
           color = "#56b4e9",linewidth = 1)+
  annotate("segment", y = 2, yend = 2, x = 5.5, xend = 6.1, 
           color = "#d55e00",linewidth = 1)+
  annotate("segment", y = 1.5, yend = 1.5, x = 5.5, xend = 6.1, 
           color = "#009e73",linewidth = 1)+
  
  annotate("text", label="slope = 0.39", x = 7.2, y =2.5,size= 5)+
  annotate("text", label="slope = 0.61", x = 7.2, y =2,size=5)+
  annotate("text", label="slope = 0.69", x = 7.2, y =1.5,size=5)+
  annotate("text",label = "conditional~R^2==0.58",
           parse = TRUE,
             #expression(conditional~R^2 == 0.48),
           x=1,y=8,size=5,hjust = 0)+
  add_phylopic(uuid = gypic, x = 7.5, y = 4, ysize =2)
m1gg



###Flowering plants ---------------------------------------------------------

dff <- d |> 
  dplyr::filter(higher_plant_group=="Flowering plants")

#### distribution graph and analysis --------------------------------------------
#####NICHE --------------------------------------------------------------------------
t <- dff |> dplyr::filter(region == "Tropical") |> 
  dplyr::select(convex) |> 
  pull()
e <- dff |> dplyr::filter(region == "Extratropical") |> 
  dplyr::select(convex)|> 
  pull()
b <- dff |> dplyr::filter(region == "Latitudinal generalists") |> 
  dplyr::select(convex)|> 
  pull()
ks.test(t,e)
ks.test(t,b)
ks.test(e,b)
d |> dunn_test(convex~region)

ndff <- dff |> 
  ggplot(aes(x=convex2,y=region,fill=region))+
  geom_density_ridges2(aes(fill=region),scale=2,color="white")+
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles=2,color="white",scale=2)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        panel.grid = element_blank())+
  scale_fill_manual(values=pal2,expand=c(0,0))+
  scale_y_discrete(expand=c(0,0.1))+
  scale_x_continuous(limits=c(0,10),breaks = c(0,2,4,6,8,10))+
  labs(x="sqrt (Niche breadth)",title="Flowering plants")+
  
  annotate("segment",x=8.6,xend=8.6, y=2,yend = 3)+
  annotate("segment",x=9.3,xend=9.3, y=2,yend = 1)+
  annotate("segment",x=10,xend=10, y=3,yend = 1)+
  
  annotate("segment",x=8.4,xend=8.6, y=3,yend = 3)+
  annotate("segment",x=8.4,xend=8.6, y=2,yend = 2)+
  
  annotate("segment",x=9.1,xend=9.3, y=2,yend = 2)+
  annotate("segment",x=9.1,xend=9.3, y=1,yend = 1)+
  
  annotate("segment",x=9.8,xend=10, y=3,yend = 3)+
  annotate("segment",x=9.8,xend=10, y=1,yend = 1)+
  
  annotate("text",label="ks ***", x=7.8, y=2.8,size=5)+
  annotate("text",label="ks ***", x=8.6, y=1.8,size=5)+
  annotate("text",label="ks ***", x=9.3, y=2.8  ,size=5)+
  
  annotate("text",label="kw ***", x=7.8, y=2.5,size=5)+
  annotate("text",label="kw ***", x=8.6, y=1.5,size=5)+
  annotate("text",label="kw ***", x=9.3, y=2.5  ,size=5)
ndff

##### RANGE ----------------------------------------------------------------------
t <- dff |> dplyr::filter(region == "Tropical") |> 
  dplyr::filter(area>100) |> 
  dplyr::select(area10) |> 
  pull()
e <- dff |> dplyr::filter(region == "Extratropical") |> 
  dplyr::filter(area>100) |> 
  dplyr::select(area10) |> 
  pull()
b <- dff |> dplyr::filter(region == "Both") |> 
  dplyr::filter(area>100) |> 
  dplyr::select(area10) |> 
  pull()
ks.test(t,e)
ks.test(t,b)
ks.test(e,b)
dff |> 
  dplyr::filter(area>100) |> 
  dunn_test(area10~region) |> 
  mutate(p2 = round(p,5))


adff <- dff |> 
  dplyr::filter(area>100) |> 
  ggplot(aes(x=area10,y=region,fill=region))+
  geom_density_ridges2(aes(fill=region),scale=2,color="white")+
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles=2,color="white",scale=2)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())+
  scale_fill_manual(values=pal2,expand=c(0,0))+
  labs(x=expression(log["10"]~ "[Range size"~(km^2)~"]"),title="Flowering plants")+
  scale_x_continuous(limits=c(2,10))+
  scale_y_discrete(expand=c(0,0.1))+
  annotate("segment",x=8.6,xend=8.6, y=2,yend = 3)+
  annotate("segment",x=9.3,xend=9.3, y=2,yend = 1)+
  annotate("segment",x=10,xend=10, y=3,yend = 1)+
  
  annotate("segment",x=8.4,xend=8.6, y=3,yend = 3)+
  annotate("segment",x=8.4,xend=8.6, y=2,yend = 2)+
  
  annotate("segment",x=9.1,xend=9.3, y=2,yend = 2)+
  annotate("segment",x=9.1,xend=9.3, y=1,yend = 1)+
  
  annotate("segment",x=9.8,xend=10, y=3,yend = 3)+
  annotate("segment",x=9.8,xend=10, y=1,yend = 1)+
  
  annotate("text",label="ks ***", x=7.8, y=2.8,size=5)+
  annotate("text",label="ks ***", x=8.6, y=1.8,size=5)+
  annotate("text",label="ks ***", x=9.3, y=2.8  ,size=5)+
  
  annotate("text",label="kw ***", x=7.8, y=2.5,size=5)+
  annotate("text",label="kw ***", x=8.6, y=1.5,size=5)+
  annotate("text",label="kw ***", x=9.3, y=2.5  ,size=5)
adff

#### modelling -----------------------------------------------------------
#SMA
dsma <- dff|> select(convex2,region,area10) |> drop_na() |> as.data.frame()
s1 <- sma(area10 ~ convex2+region,data=dsma, multcomp=TRUE,type="elevation")
multcompmatrix(s1)
s2 <- sma(area10 ~ convex2+region,data=dsma, multcomp=TRUE,type="shift")
multcompmatrix(s2)
s3 <- sma(area10 ~ convex2*region,data=dsma, multcomp=TRUE)
multcompmatrix(s3)


#glmm
m1ff <- lmer(data=dff, area10 ~ convex2:region + (1|region) )
r.squaredGLMM(m1ff)
m1ffcoef <- tidy(m1ff,"ran_coefs") |> print()

tff <- pairs(emtrends(m1ff,"region",var="convex2"),pbkrtest.limit = 6769) |> 
  tidy() |> 
  mutate(group="Flowering plants")

dff <- dff |> mutate(fit.m = predict(m1ff, re.form=NA),
                   fit.c = predict(m1ff, re.form=NULL))

treepic <- get_uuid(name="Melicope", n = 1)
m1ffg <-dff |> 
  #dplyr::filter(area >100) |> 
  ggplot(aes(y=area10,x=convex2,color=region))+
  geom_point(color="gray70",alpha = 0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank())+
  geom_line(aes(y=fit.c,col=region),linewidth=1)+
  scale_color_manual(values =c("#56b4e9","#d55e00","#009e73"))+
  #geom_abline(intercept = 4.42, slope = 0.404,color="#56b4e9", linewidth=1)+ #Both
  #geom_abline(intercept = 3.29, slope = 0.688,color="#d55e00", linewidth=1)+ #Extratropical
  #geom_abline(intercept = 3.02, slope = 0.698,color="#009e73", linewidth=1)+ #Tropical
  labs(y=expression("log10 [Range size"~(km^2)~"]"),x="sqrt (Niche breadth)",title="Flowering plants")+
  scale_y_continuous(limits=c(1,8))+
  scale_x_continuous(limits=c(1,8))+
  #annotate("text",label = expression(conditional~R^2 == 0.50),
  #         x=0,y=8,size=5,hjust = 0)+
  #legend
  annotate("segment", y = 2.5, yend = 2.5, x = 5.5, xend = 6.1, 
           color = "#56b4e9",linewidth = 1)+
  annotate("segment", y = 2, yend = 2, x = 5.5, xend = 6.1, 
           color = "#d55e00",linewidth = 1)+
  annotate("segment", y = 1.5, yend = 1.5, x = 5.5, xend = 6.1, 
           color = "#009e73",linewidth = 1)+
  
  annotate("text", label="slope = 0.40", x = 7.2, y =2.5,size= 5)+
  annotate("text", label="slope = 0.67", x = 7.2, y =2,size=5)+
  annotate("text", label="slope = 0.70", x = 7.2, y =1.5,size=5)+
  annotate("text",label ="conditional~R^2==0.49",
           parse = TRUE,
             #expression(conditional~R^2 == 0.49),
           x=1,y=8,size=5,hjust = 0)+
  add_phylopic(uuid = treepic, x = 7.5, y = 4, ysize = 1)
m1ffg


#### Save composite graph  --------------------------------------------------

cp=ndb + ndf + ndg + ndff+
adb + adf + adg + adff+
  plot_layout(guides = 'collect',ncol = 4,axis_titles = "collect")+
  plot_annotation(tag_levels = "a",tag_suffix = ")")
  ggsave(here("products","RESU_dist_region.jpeg"),
         units="cm",dpi = 300,height = 25,width = 50,plot=cp)

cp=m1bg + m1fg + m1ffg + m1gg +
  plot_layout(guides = 'collect')+
  plot_annotation(tag_levels = "a",tag_suffix = ")")&
  theme(legend.position="bottom")
ggsave(here("products","RESU_model1_v3.jpeg"),
       units="cm",dpi = 300,height = 30,width = 30,plot=cp)
ggsave(here("products","RESU_model1_v3.tiff"),
       units="cm",dpi = 300,height = 30,width = 30, plot=cp)

# make the results table

ctm1 <- rbind(tb,tg,tf,tff) |> 
  dplyr::mutate(p.value = round(adj.p.value,3)) |> 
  dplyr::select(contrast,p.value, group) |> 
  dplyr::relocate(group,.before = contrast) |> 
  dplyr::rename(region = contrast) |> 
  flextable::flextable() %>%                     
  flextable::autofit() %>%    
  flextable::save_as_docx(path = here("products","Table_model1.docx"))


## 2) model 2: niche breadth per elevation --------------------------------

pal <- c("#e69f00","#0072b2")

### All species ----------------------------------------------------------------------------
#### Elevation vs RANGE --------------------------------------------------------------------

# Piecewise regressions
fit <- lm(data=d, area~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
frs <-  as.data.frame(cbind(d$elevation, s$fitted.values))
colnames(frs) <-  c("elevation","area")

# plot
drs <- d |> 
  dplyr::filter(area>100) |> 
  dplyr::filter(convex>0) |> 
  ggplot(aes(x=elevation,y=area,color=ele2))+
  geom_point(alpha=0.5)+
  #geom_line(data=frs,color="black",linewidth=1)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma)+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(x="Elevation (m)",y=expression("Range size"~(km^2)))+
  annotate("text",label = "beta[1] == 2516 * ';' ~ beta[2] == -2668",
           parse = TRUE,
           #label = expression(beta ~"1" ~ "= 2516" ~";" ~ beta ~"2 = -2668"),
          x=100,y=30000000,size=5,hjust=0)+
  annotate("text", label = "psi = 235",
           x=100,y=29000000,size=5,hjust=0)           
drs  

t1 <- d |>   dplyr::filter(ele2 =="Low elevation") |> 
  dplyr::select(area10) |> 
  pull()
t2 <- d |>   dplyr::filter(ele2 =="Mountain") |> 
  dplyr::select(area10) |> 
  pull()
ks.test(t1,t2)
d |>   dunn_test(area10~ele2,detailed=T)


adele <- d |> 
  filter(area>100) |> 
  ggplot(aes(x=area10,y=e))+
  geom_density_ridges(aes(scale=5,color=paste(l,m)),fill="transparent",linewidth=1)+
  #stat_density_ridges(quantile_lines = TRUE,
  #                    quantiles=2,color="white",scale=5,alpha=0.8)+
  theme_bw()+
  theme(text=element_text(size=12,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        rect = element_rect(fill = "transparent"),
        plot.background = element_rect(fill="transparent"),
        panel.background = element_rect(fill="transparent"),
        panel.grid = element_blank())+
  #scale_fill_manual(values=pal,expand=c(0,0))+
  scale_color_manual(values=pal,expand=c(0,0))+
  labs(x=expression(log["10"]~"[Range size"~(km^2)~"]"))+
  scale_x_continuous(limits=c(0,10),expand=c(0,0.1))+
  scale_y_discrete(expand=c(0,0))+
  annotate("text",label=
             "ks ***
kw ***", 
           x=8, y=1.5  ,size=3,hjust=0)
adele


#### Elevation vs. NICHE ---------------------------------------------------------
# Piecewise regressions
fit <- lm(data=d, convex~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(d$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","convex")

# plot
dnb <- d |> 
  dplyr::filter(area>100) |> 
  dplyr::filter(convex>0) |> 
  ggplot(aes(x=elevation,y=convex,color=ele2))+
  geom_point(alpha=0.5)+
  geom_line(data=fnb,color="black",linewidth=1)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma)+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(x="Elevation (m)",y="Niche breadth")+
  annotate("text", label = "beta[1] == 0.01 * ';' ~ beta[2] == -0.01",
           parse = TRUE,
           #label = expression(beta ~"1" ~ "= 0.01" ~";" ~ beta ~"2 = -0.01"),
           x=100,y=70,size=5,hjust = 0)+
  annotate("text", label = "psi = 235",
           x=100,y=67,size=5,hjust=0)     
dnb  


ndele <- d |> 
  filter(area>100) |> 
  ggplot(aes(x=convex2,y=e))+
  geom_density_ridges(aes(scale=5,color=paste(l,m)),fill="transparent",linewidth=1)+
  #stat_density_ridges(quantile_lines = TRUE,
  #                    quantiles=2,color="white",scale=5,alpha=0.8)+
  theme_bw()+
  theme(text=element_text(size=12,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        rect = element_rect(fill = "transparent"),
        plot.background = element_rect(fill="transparent"),
        panel.background = element_rect(fill="transparent"),
        panel.grid = element_blank())+
  #scale_fill_manual(values=pal,expand=c(0,0))+
  scale_color_manual(values=pal,expand=c(0,0))+
  labs(x="sqrt(Niche breadth)")+
  scale_x_continuous(limits=c(0,10),expand=c(0,0.1))+
  scale_y_discrete(expand=c(0,0))+
  annotate("text",label=
  "ks ***
kw ***", 
           x=8, y=1.5  ,size=3,hjust=0)
ndele


####Make composite graph -------------------------------------------------------
cp=(dnb + inset_element(ndele, left = 0.6, bottom = 0.6, right = 0.98, top = 0.9, align_to = 'full',on_top = T))+
  (drs + inset_element(adele, left = 0.6, bottom = 0.6, right = 0.98, top = 0.9, align_to = 'full',on_top = T))+
  plot_layout(guides = 'collect',ncol = 2,axis_titles = "collect",tag_level = "new")+
  plot_annotation(tag_levels = list(c("a)","","b)")))
ggsave(here("products","RESU_dist_ele_allPlants_v3.jpeg"),
       units="cm",dpi = 300,height = 15,width = 32,plot=cp)
ggsave(here("products","RESU_dist_ele_allPlants_v3.tiff"),
       units="cm",dpi = 300,height = 15,width = 32,plot=cp)

### All species per region-------------------------------------------
#### NICHE ---------------------------------------------------------
t1 <- d |> dplyr::filter(region == "Tropical") |> 
  dplyr::filter(ele2 =="Low elevation") |> 
  dplyr::select(convex) |> 
  pull()
t2 <- d |> dplyr::filter(region == "Tropical") |> 
  dplyr::filter(ele2 =="Mountain") |> 
  dplyr::select(convex) |> 
  pull()
ks.test(t1,t2)
d |> dplyr::filter(region == "Tropical") |> 
  dunn_test(convex~ele2)

ndt <- d |> 
  dplyr::filter(region == "Tropical") |> 
  ggplot(aes(x=convex,y=ele2,fill=ele2))+
  geom_density_ridges2(aes(fill=ele2),scale=2,color="white")+
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles=2,color="white",scale=2)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        panel.grid = element_blank())+
  scale_fill_manual(values=pal,expand=c(0,0))+
  scale_y_discrete(expand=c(0,0.1))+
  #scale_y_continuous(limits=c(1,8))+
  scale_x_continuous(limits=c(0,60))+
  labs(x="Niche breadth",title = "Tropical region")+
  annotate("segment",x=60,xend=60, y=2,yend = 1)+
  
  annotate("segment",x=58,xend=60, y=2,yend = 2)+
  annotate("segment",x=58,xend=60, y=1,yend = 1)+
  
  annotate("text",label="ks ***", x=55, y=1.8  ,size=5)+
  annotate("text",label="kw ***", x=55, y=1.5  ,size=5)
ndt

# inset graph
indt <- d |>  filter(convex !=0 ) |> 
  dplyr::filter(region == "Tropical") |> 
  group_by(cEle,ele2) |> 
  summarise(maxConvex = max(convex),mConvex=median(convex)) |> 
  filter(maxConvex !=0) |> 
  ggplot(aes(x=cEle,y=mConvex,fill=ele2))+
  geom_col(alpha=0.6)+
  #geom_line(aes(x=cLat,y=mConvex,color=region),linewidth=1)+
  #facet_wrap(~fct_rev(ele2))+
  theme_bw()+
  theme(text=element_text(size=12,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        rect = element_rect(fill = "transparent"),
        axis.title.x = element_blank())+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000),expand = c(0,0))+
  scale_y_continuous(limits=c(0,8),breaks = c(0,4,8))+
  #scale_y_log10(labels=trans_format("log10", math_format(10^.x)))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(x="Elevation",y="Niche breadth")
  coord_flip()
indt

d |>  filter(convex !=0 ) |> 
  dplyr::filter(region == "Tropical") |> 
  group_by(cEle,ele2) |> 
  summarise(maxConvex = median(convex),mConvex=median(convex)) 

e1 <- d |> dplyr::filter(region == "Extratropical") |> 
  dplyr::filter(ele2 =="Low elevation") |> 
  dplyr::select(convex)|>
  pull()
e2 <- d |> dplyr::filter(region == "Extratropical") |> 
  dplyr::filter(ele2 =="Mountain") |> 
  dplyr::select(convex)|>
  pull()
ks.test(e1,e2)
d |> dplyr::filter(region == "Extratropical") |> 
  dunn_test(convex~ele2)

nde <- d |> 
  dplyr::filter(region == "Extratropical") |> 
  ggplot(aes(x=convex,y=ele2,fill=ele2))+
  geom_density_ridges2(aes(fill=ele2),scale=2,color="white")+
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles=2,color="white",scale=2)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())+
  scale_fill_manual(values=pal,expand=c(0,0))+
  #scale_y_continuous(limits=c(1,8))+
  scale_y_discrete(expand=c(0,0.1))+
  scale_x_continuous(limits=c(0,60))+
  labs(x="Niche breadth",title = "Extratropical region")+
  annotate("segment",x=60,xend=60, y=2,yend = 1)+
  
  annotate("segment",x=58,xend=60, y=2,yend = 2)+
  annotate("segment",x=58,xend=60, y=1,yend = 1)+
  
  annotate("text",label="ks ***", x=55, y=1.8  ,size=5)+
  annotate("text",label="kw ***", x=55, y=1.5  ,size=5)
nde

# inset graph
inde <- d |>  filter(convex !=0 ) |> 
  dplyr::filter(region == "Extratropical") |> 
  group_by(cEle,ele2) |> 
  summarise(maxConvex = max(convex),mConvex=median(convex)) |> 
  filter(maxConvex !=0) |> 
  ggplot(aes(x=cEle,y=mConvex,fill=ele2))+
  geom_col(alpha=0.6)+
  #geom_line(aes(x=cLat,y=mConvex,color=region),linewidth=1)+
  #facet_wrap(~fct_rev(ele2))+
  theme_bw()+
  theme(text=element_text(size=12,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        rect = element_rect(fill = "transparent"),
        axis.title.x = element_blank())+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000),expand = c(0,0))+
  scale_y_continuous(limits=c(0,8),breaks = c(0,4,8))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(x="Elevation",y="Niche breadth")+
  coord_flip()
inde

b1 <- d |> dplyr::filter(region == "Both") |> 
  dplyr::filter(ele2 =="Low elevation") |> 
  dplyr::select(convex)|> 
  pull()
b2 <- d |> dplyr::filter(region == "Both") |> 
  dplyr::filter(ele2 =="Mountain") |> 
  dplyr::select(convex)|> 
  pull()
ks.test(b1,b2)
d |> dplyr::filter(region == "Both") |> 
  dunn_test(convex~ele2)

ndb <- d|> 
  dplyr::filter(region == "Both") |> 
  ggplot(aes(x=convex,y=ele2,fill=ele2))+
  geom_density_ridges2(aes(fill=ele2),scale=2,color="white")+
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles=2,color="white",scale=2)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())+
  scale_fill_manual(values=pal,expand=c(0,0))+
  scale_y_discrete(expand=c(0,0.1))+
  scale_x_continuous(limits=c(0,60))+
  labs(x="Niche breadth",title= "Both regions")+
  annotate("segment",x=60,xend=60, y=2,yend = 1)+
  
  annotate("segment",x=58,xend=60, y=2,yend = 2)+
  annotate("segment",x=58,xend=60, y=1,yend = 1)+
  
  annotate("text",label="ks ***", x=55, y=1.8  ,size=5)+
  annotate("text",label="kw ***", x=55, y=1.5  ,size=5)
ndb

# inset graph
indb <- d |>  filter(convex !=0 ) |> 
  dplyr::filter(region == "Both") |> 
  group_by(cEle,ele2) |> 
  summarise(maxConvex = max(convex),mConvex=median(convex)) |> 
  filter(maxConvex !=0) |> 
  ggplot(aes(x=cEle,y=mConvex,fill=ele2))+
  geom_col(alpha=0.6)+
  theme_bw()+
  theme(text=element_text(size=12,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        rect = element_rect(fill = "transparent"),
        axis.title.x = element_blank())+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000),expand = c(0,0))+
  scale_y_continuous(limits=c(0,8),breaks = c(0,4,8))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(x="Elevation",y="Niche breadth")+
  coord_flip()
indb

####RANGE -------------------------------------------------------------------------
t1 <- d |> dplyr::filter(region == "Tropical") |> 
  dplyr::filter(ele2 =="Low elevation") |> 
  dplyr::select(area10) |> 
  pull()
t2 <- d |> dplyr::filter(region == "Tropical") |> 
  dplyr::filter(ele2 =="Mountain") |> 
  dplyr::select(area10) |> 
  pull()
ks.test(t1,t2)
d |> dplyr::filter(region == "Tropical") |> 
  dunn_test(area10~ele2,detailed=T)

adt <- d |> 
  dplyr::filter(region == "Tropical") |> 
  ggplot(aes(x=area10,y=ele2,fill=ele2))+
  geom_density_ridges2(aes(fill=ele2),scale=2,color="white")+
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles=2,color="white",scale=2)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        rect = element_rect(fill = "transparent"),
        plot.background = element_rect(fill="transparent"),
        panel.background = element_rect(fill="transparent"),
        panel.grid = element_blank())+
  scale_fill_manual(values=pal,expand=c(0,0))+
  labs(x=expression("log10 [Range size"~(km^2)~"]"),title="Tropical region")+
  scale_x_continuous(limits=c(2,10))+
  scale_y_discrete(expand=c(0,0.1))+
  annotate("segment",x=10,xend=10, y=2,yend = 1)+
  
  annotate("segment",x=9.8,xend=10, y=2,yend = 2)+
  annotate("segment",x=9.8,xend=10, y=1,yend = 1)+
  
  annotate("text",label="ks ***", x=9.3, y=1.8  ,size=5)+
  
  annotate("text",label="kw ***", x=9.3, y=1.5  ,size=5)
adt

# inset graph
iadt <- d |>  filter(convex !=0 ) |> 
  dplyr::filter(region == "Tropical") |> 
  group_by(cEle,ele2) |> 
  summarise(maxConvex = max(area10),mConvex=median(area10)) |> 
  filter(maxConvex !=0) |> 
  ggplot(aes(x=cEle,y=mConvex,fill=ele2))+
  geom_col(alpha=0.6)+
  #geom_line(aes(x=cLat,y=mConvex,color=region),linewidth=1)+
  #facet_wrap(~fct_rev(ele2))+
  theme_bw()+
  theme(text=element_text(size=12,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"),
        axis.title.x = element_blank())+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000),expand = c(0,0))+
  scale_y_continuous(limits=c(0,6),breaks = c(0,3,6))+
  #scale_y_log10(labels=trans_format("log10", math_format(10^.x)))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(x="Elevation")+
  coord_flip()
iadt


t1 <- d |> dplyr::filter(region == "Extratropical") |> 
  dplyr::filter(ele2 =="Low elevation") |> 
  dplyr::select(area10) |> 
  pull()
t2 <- d |> dplyr::filter(region == "Extratropical") |> 
  dplyr::filter(ele2 =="Mountain") |> 
  dplyr::select(area10) |> 
  pull()
ks.test(t1,t2)
d |> dplyr::filter(region == "Extratropical") |> 
  dunn_test(area10~ele2)

ade <- d |> 
  dplyr::filter(region == "Extratropical") |> 
  ggplot(aes(x=area10,y=ele2,fill=ele2))+
  geom_density_ridges2(aes(fill=ele2),scale=2,color="white")+
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles=2,color="white",scale=2)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())+
  scale_fill_manual(values=pal,expand=c(0,0))+
  labs(x=expression("log10 [Range size"~(km^2)~"]"),title="Extratropical region")+
  scale_x_continuous(limits=c(2,10))+
  scale_y_discrete(expand=c(0,0.1))+
  annotate("segment",x=10,xend=10, y=2,yend = 1)+
  
  annotate("segment",x=9.8,xend=10, y=2,yend = 2)+
  annotate("segment",x=9.8,xend=10, y=1,yend = 1)+
  
  annotate("text",label="ks ***", x=9.3, y=1.8  ,size=5)+
  
  annotate("text",label="kw ***", x=9.3, y=1.5  ,size=5)
ade

# inset graph
iade <- d |>  filter(convex !=0 ) |> 
  dplyr::filter(region == "Extratropical") |> 
  group_by(cEle,ele2) |> 
  summarise(maxConvex = max(area10),mConvex=median(area10)) |> 
  filter(maxConvex !=0) |> 
  ggplot(aes(x=cEle,y=mConvex,fill=ele2))+
  geom_col(alpha=0.6)+
  #geom_line(aes(x=cLat,y=mConvex,color=region),linewidth=1)+
  #facet_wrap(~fct_rev(ele2))+
  theme_bw()+
  theme(text=element_text(size=12,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        rect = element_rect(fill = "transparent"),
        axis.title.x = element_blank())+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000),expand = c(0,0))+
  scale_y_continuous(limits=c(0,6),breaks = c(0,3,6))+
  #scale_y_log10(labels=trans_format("log10", math_format(10^.x)))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(x="Elevation")+
  coord_flip()
iade

t1 <- d |> dplyr::filter(region == "Both") |> 
  dplyr::filter(ele2 =="Low elevation") |> 
  dplyr::select(area10) |> 
  pull()
t2 <- d |> dplyr::filter(region == "Both") |> 
  dplyr::filter(ele2 =="Mountain") |> 
  dplyr::select(area10) |> 
  pull()
ks.test(t1,t2)
d |> dplyr::filter(region == "Both") |> 
  dunn_test(area10~ele2)

adb <- d |> 
  dplyr::filter(region == "Both") |> 
  ggplot(aes(x=area10,y=ele2,fill=ele2))+
  geom_density_ridges2(aes(fill=ele2),scale=2,color="white")+
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles=2,color="white",scale=2)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())+
  scale_fill_manual(values=pal,expand=c(0,0))+
  labs(x=expression("log10 [Range size"~(km^2)~"]"),title="Both regions")+
  scale_x_continuous(limits=c(2,10))+
  scale_y_discrete(expand=c(0,0.1))+
  annotate("segment",x=10,xend=10, y=2,yend = 1)+
  
  annotate("segment",x=9.8,xend=10, y=2,yend = 2)+
  annotate("segment",x=9.8,xend=10, y=1,yend = 1)+
  
  annotate("text",label="ks ***", x=9.3, y=1.8  ,size=5)+
  
  annotate("text",label="kw ***", x=9.3, y=1.5  ,size=5)
adb

# inset graph
iadb <- d |>  filter(convex !=0 ) |> 
  dplyr::filter(region == "Both") |> 
  group_by(cEle,ele2) |> 
  summarise(maxConvex = max(area10),mConvex=median(area10)) |> 
  filter(maxConvex !=0) |> 
  ggplot(aes(x=cEle,y=mConvex,fill=ele2))+
  geom_col(alpha=0.6)+
  #geom_line(aes(x=cLat,y=mConvex,color=region),linewidth=1)+
  #facet_wrap(~fct_rev(ele2))+
  theme_bw()+
  theme(text=element_text(size=12,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill="transparent"),
        axis.title.x = element_blank())+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000),expand = c(0,0))+
  scale_y_continuous(limits=c(0,6),breaks = c(0,3,6))+
  #scale_y_log10(labels=trans_format("log10", math_format(10^.x)))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(x="Elevation")+
  coord_flip()
iadb

###### Make composite graph --------------------------------------------------------------------

(ndt + inset_element(indt, left = 0.6, bottom = 0.6, right = 0.98, top = 0.88, align_to = 'full',on_top = T))+
(nde+inset_element(inde, left = 0.6, bottom = 0.6, right = 0.98, top = 0.88, align_to = 'full',on_top = T))+
(ndb+inset_element(indb, left = 0.6, bottom = 0.6, right = 0.98, top = 0.88, align_to = 'full',on_top = T))+
  (adt+inset_element(iadt, left = 0.6, bottom = 0.6, right = 0.98, top = 0.88, align_to = 'full',on_top = F))+
  (ade+inset_element(iade, left = 0.6, bottom = 0.6, right = 0.98, top = 0.88, align_to = 'full',on_top = T))+
  (adb+inset_element(iadb, left = 0.6, bottom = 0.6, right = 0.98, top = 0.88, align_to = 'full',on_top = T))+
  plot_layout(guides = 'collect',ncol = 3,axis_titles = "collect",tag_level = "new")+
  plot_annotation(tag_levels = list(c("a)","","b)","","c)",
                                      "","d)","","e)","","f)")))
ggsave(here("products","RESU_dist_ele_allPlants_v2.jpeg"),
       units="cm",dpi = 300,height = 25,width = 40)



ndt + nde + ndb+
adt + ade + adb +
  plot_layout(guides = 'collect',ncol = 3,axis_titles = "collect")+
  plot_annotation(tag_levels = "a",tag_suffix = ")")
ggsave(here("products","RESU_dist_ele_allPlants.jpeg"),
       units="cm",dpi = 300,height = 25,width = 35)

#### SMA -------------------------------------------------------------------

dsma <- d |> select(convex2,region,area10,ele2,region) |>
  unite(geo2, c("region","ele2")) |> 
  drop_na() |> as.data.frame()
s1 <- sma(area10 ~ convex2+geo2,data=dsma,multcomp=TRUE,type="elevation")
multcompmatrix(s1)
summary(s1)
s2 <- sma(area10 ~ convex2+geo2,data=dsma, multcomp=TRUE,type="shift")
multcompmatrix(s2)
s3 <- sma(area10 ~ convex2*geo2,data=dsma, multcomp=TRUE)
multcompmatrix(s3)

### Bryophytes --------------------------------------------------------------

db <- d |> 
  dplyr::filter(higher_plant_group=="Bryophytes")

#### distribution graph and analysis --------------------------------------------
#####NICHE -------------------------------------------------------------------------

# Piecewise regressions
db1 <- db |> 
  dplyr::filter(region == "Tropical")
fit <- lm(data=db1, convex~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","convex")


ndbt <- db |> 
  dplyr::filter(region == "Tropical") |> 
  ggplot(aes(x=elevation,y=convex,color=ele2))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
 scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma)+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(y="Niche breadth",x="Elevation (m)",title="Bryophytes",
       subtitle = "Tropical region")+
  geom_line(data=fnb,color="black",linewidth=1)+
  annotate("text",label = expression(beta ~"1" ~ "= 0.01" ~";" ~ beta ~"2 = -0.01"),
           x=100,y=70,size=5,hjust = 0)+
  annotate("text", label = "psi = 434",
           x=100,y=65,size=5,hjust=0)  
ndbt

# Piecewise regressions
db1 <- db |> 
  dplyr::filter(region == "Extratropical")
fit <- lm(data=db1, convex~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","convex")


ndbe <- db |> 
  dplyr::filter(region == "Extratropical") |> 
  ggplot(aes(x=elevation,y=convex,color=ele2))+
  geom_point(alpha=0.5)+
  geom_line(data=fnb,color="black",linewidth=1)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma)+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(y="Niche breadth",x="Elevation (m)",title="Bryophytes",
       subtitle = "Extratropical region")+
  annotate("text",label = expression(beta ~"1" ~ "= 0.02" ~";" ~ beta ~"2 = -0.02"),
           x=100,y=70,size=5,hjust = 0)+
  annotate("text", label = "psi = 439",
           x=100,y=65,size=5,hjust=0)  
ndbe

# Piecewise regressions
db1 <- db |> 
  dplyr::filter(region == "Latitudinal generalists")
fit <- lm(data=db1, convex~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","convex")


ndbb <- db |> 
  dplyr::filter(region == "Latitudinal generalists") |> 
  ggplot(aes(x=elevation,y=convex,color=ele2))+
  geom_point(alpha=0.5)+
  geom_line(data=fnb,color="black",linewidth=1)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma)+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(y="Niche breadth",x="Elevation (m)",title="Bryophytes",
       subtitle = "Tropical and Extratropical")+
  annotate("text",label = expression(beta ~"1" ~ "= 0.02" ~";" ~ beta ~"2 = -0.02"),
           x=100,y=70,size=5,hjust = 0)+
  annotate("text", label = "psi = 390",
           x=100,y=65,size=5,hjust=0)  
ndbb

#####RANGE -------------------------------------------------------------------------
db1 <- db |> 
  dplyr::filter(region == "Tropical")
fit <- lm(data=db1, area~elevation )
summary(fit)
s <- segmented(fit)
print(s)

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","area")

adbt <- db |> 
  dplyr::filter(region == "Tropical") |> 
  ggplot(aes(x=elevation,y=area,color=ele2))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma, limits=c(0,20000000))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(x="Elevation (m)",y=expression("Range size"~(km^2)),
       title="Bryophytes",subtitle = "Tropical region")+
  geom_line(data=fnb,color="black",linewidth=1)+
  annotate("text",label = expression(beta ~"1" ~ "= -496" ~";" ~ beta ~"2 = 484"),
           x=100,y=20000000,size=5,hjust=0)+
  annotate("text", label = "psi = 333",
           x=100,y=18000000,size=5,hjust=0)   
adbt

db1 <- db |> 
  dplyr::filter(region == "Extratropical")
fit <- lm(data=db1, area~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","area")

adbe <- db |> 
  dplyr::filter(region == "Extratropical") |> 
  ggplot(aes(x=elevation,y=area,color=ele2))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma, limits=c(0,20000000))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(x="Elevation (m)",y=expression("Range size"~(km^2)),
       title="Bryophytes",subtitle = "Extratropical region")+
geom_line(data=fnb,color="black",linewidth=1)+
  annotate("text",label = expression(beta ~"1" ~ "= 5096" ~";" ~ beta ~"2 = -5458"),
           x=100,y=20000000,size=5,hjust=0)+
  annotate("text", label = "psi = 260",
           x=100,y=18000000,size=5,hjust=0)   
adbe


db1 <- db |> 
  dplyr::filter(region == "Latitudinal generalists")
fit <- lm(data=db1, area~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","area")


adbb <- db |> 
  dplyr::filter(region == "Latitudinal generalists") |> 
  ggplot(aes(x=elevation,y=area,color=ele2))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma, limits=c(0,20000000))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(x="Elevation (m)",y=expression("Range size"~(km^2)),
       title="Bryophytes",subtitle = "Tropical and Extratropical")+
  geom_line(data=fnb,color="black",linewidth=1)+
  annotate("text",label = expression(beta ~"1" ~ "= 4153" ~";" ~ beta ~"2 = -4528"),
           x=100,y=20000000,size=5,hjust=0)+
  annotate("text", label = "psi = 302",
           x=100,y=18000000,size=5,hjust=0)   
adbb

ndbt + ndbe + ndbb
adbt + adbe + adbb
#### modelling ----------------------------------------------------------------------

db <- d |> 
  dplyr::filter(higher_plant_group=="Bryophytes")
m2b <- lmer(data=db, area10 ~ convex2:ele2:region+(1|region))
r.squaredGLMM(m2b)
m2bcoef <- tidy(m2b,"ran_coefs") |> print()
coef(m2b)
tidy(m2b)

bp2 <- pairs(emtrends(m2b,c("ele2","region"),var="convex2"),pbkrtest.limit = 6768)|>
  tidy() |> 
  mutate(group = "Bryophytes")

db <- db |> mutate(fit.m = predict(m2b, re.form=NA),
                     fit.c = predict(m2b, re.form=NULL))

m1bgET <-db |> 
  dplyr::filter(region=="Tropical") |> 
  ggplot(aes(y=area10,x=convex2,color=ele2))+
  geom_point(color="gray70",alpha = 0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank())+
  geom_line(aes(y=fit.c,col=ele2),linewidth=1)+
  scale_color_manual(values =c("#56b4e9","#d55e00"))+
  #geom_abline(intercept = 3.68, slope = 0.569,color="#56b4e9", linewidth=1)+ #Both
  #geom_abline(intercept = 2.69, slope = 0.844,color="#d55e00", linewidth=1)+ #Extratropical
  #geom_abline(intercept = 2.72, slope = 0.675,color="#009e73", linewidth=1)+ #Tropical
  labs(y=expression("log10 [Range size"~(km^2)~"]"),x="Niche breadth",title="Bryophytes",
       subtitle = "Tropical")+
  scale_y_continuous(limits=c(1,8))+
  scale_x_continuous(limits=c(1,8))+
#legend
  annotate("segment", y = 2.5, yend = 2.5, x = 5.5, xend = 6.1, 
         color = "#56b4e9",linewidth = 1)+
  annotate("segment", y = 2, yend = 2, x = 5.5, xend = 6.1, 
           color = "#d55e00",linewidth = 1)+
  
  annotate("text", label="slope = 0.60", x = 7.2, y =2.5,size= 5)+
  annotate("text", label="slope = 0.54", x = 7.2, y =2,size=5)
  #annotate("text",label = expression(conditional~R^2 == 0.44),
  #         x=1,y=8,size=5,hjust = 0)
m1bgET

m1bgEE <-db |> 
  dplyr::filter(region=="Extratropical") |> 
  ggplot(aes(y=area10,x=convex2,color=ele2))+
  geom_point(color="gray70",alpha = 0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank())+
  geom_line(aes(y=fit.c,col=ele2),linewidth=1)+
  scale_color_manual(values =c("#56b4e9","#d55e00"))+
  #geom_abline(intercept = 3.68, slope = 0.569,color="#56b4e9", linewidth=1)+ #Both
  #geom_abline(intercept = 2.69, slope = 0.844,color="#d55e00", linewidth=1)+ #Extratropical
  #geom_abline(intercept = 2.72, slope = 0.675,color="#009e73", linewidth=1)+ #Tropical
  labs(y=expression("log10 [Range size"~(km^2)~"]"),x="Niche breadth",title="Bryophytes",
       subtitle = "Extratropical")+
  scale_y_continuous(limits=c(1,8))+
  scale_x_continuous(limits=c(1,8))+
  #legend
  annotate("segment", y = 2.5, yend = 2.5, x = 5.5, xend = 6.1, 
           color = "#56b4e9",linewidth = 1)+
  annotate("segment", y = 2, yend = 2, x = 5.5, xend = 6.1, 
           color = "#d55e00",linewidth = 1)+
  
  annotate("text", label="slope = 0.73", x = 7.2, y =2.5,size= 5)+
  annotate("text", label="slope = 0.60", x = 7.2, y =2,size=5)
#annotate("text",label = expression(conditional~R^2 == 0.44),
#         x=1,y=8,size=5,hjust = 0)
m1bgEE

m1bgEB <-db |> 
  dplyr::filter(region=="Both") |> 
  ggplot(aes(y=area10,x=convex2,color=ele2))+
  geom_point(color="gray70",alpha = 0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank())+
  geom_line(aes(y=fit.c,col=ele2),linewidth=1)+
  scale_color_manual(values =c("#56b4e9","#d55e00"))+
  #geom_abline(intercept = 3.68, slope = 0.569,color="#56b4e9", linewidth=1)+ #Both
  #geom_abline(intercept = 2.69, slope = 0.844,color="#d55e00", linewidth=1)+ #Extratropical
  #geom_abline(intercept = 2.72, slope = 0.675,color="#009e73", linewidth=1)+ #Tropical
  labs(y=expression("log10 [Range size"~(km^2)~"]"),x="Niche breadth",title="Bryophytes",
       subtitle = "Both")+
  scale_y_continuous(limits=c(1,8))+
  scale_x_continuous(limits=c(1,8))+
  #legend
  annotate("segment", y = 2.5, yend = 2.5, x = 5.5, xend = 6.1, 
           color = "#56b4e9",linewidth = 1)+
  annotate("segment", y = 2, yend = 2, x = 5.5, xend = 6.1, 
           color = "#d55e00",linewidth = 1)+
  
  annotate("text", label="slope = 0.51", x = 7.2, y =2.5,size= 5)+
  annotate("text", label="slope = 0.45", x = 7.2, y =2,size=5)
#annotate("text",label = expression(conditional~R^2 == 0.44),
#         x=1,y=8,size=5,hjust = 0)
m1bgEB

dbl <- db |> dplyr::filter(ele2 == "Low elevation") |> 
  dplyr::select(convex2) |> 
  pull()
dbm <- db |> dplyr::filter(ele2 == "Mountain") |> 
  dplyr::select(convex2)|> 
  pull()
ks.test(dbl,dbm)

vtbe1 <- db |>
  ggplot(aes(x=ele2,y=convex2,fill=ele2))+
  geom_flat_violin(position = position_nudge(x = .12, y = 0), alpha = 0.6)+
  geom_jitter(aes(color = ele2),
              width=0.08, alpha = 0.6)+
  #facet_wrap(.~region)+
  #scale_y_continuous(limits=c(1,8))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none",
        #axis.text.y = element_blank(),
        strip.background = element_rect(colour="white", fill="white"))+
  #labs(y=expression("log10 [Range size"~(km^2)~"]"), x="")+
  labs(y="Niche breadth", x="",title = "Bryophytes")+
  scale_color_grey()+
  scale_fill_grey()+
  scale_y_continuous(limits=c(1,8))+
  coord_flip()+
  annotate("text", label="p-value < 0.001", x = 0.7, y =6,size=5)
vtbe1 

dbl <- db |> dplyr::filter(ele2 == "Low elevation") |> 
  dplyr::select(area10) |> 
  pull()
dbm <- db |> dplyr::filter(ele2 == "Mountain") |> 
  dplyr::select(area10)|> 
  pull()
ks.test(dbl,dbm)

vtbe2 <- db |>
  ggplot(aes(x=ele2,y=area10,fill=ele2))+
  geom_flat_violin(position = position_nudge(x = .12, y = 0), alpha = 0.6)+
  geom_jitter(aes(color = ele2),
              width=0.08, alpha = 0.6)+
 #facet_wrap(.~region)+
  #scale_y_continuous(limits=c(1,8))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.text.y = element_blank(),
        strip.background = element_rect(colour="white", fill="white"))+
  #labs(y=expression("log10 [Range size"~(km^2)~"]"), x="")+
  labs(y=expression("log10 [Range size"~(km^2)~"]"), x="",title="Bryophytes")+
  scale_color_grey()+
  scale_fill_grey()+
  scale_y_continuous(limits=c(2,8))+
  coord_flip()+
  annotate("text", label="p-value < 0.001", x = 0.7, y =6,size=5)
vtbe2 

### Ferns ------------------------------------------------
df <- d |> 
  dplyr::filter(higher_plant_group=="Ferns and lycophytes")

#### distribution graph and analysis --------------------------------------------
#####NICHE -------------------------------------------------------------------------
# Piecewise regressions
db1 <- df |> 
  dplyr::filter(region == "Tropical")
fit <- lm(data=db1, convex~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","convex")

ndft <- df |> 
  dplyr::filter(region == "Tropical") |> 
  ggplot(aes(x=elevation,y=convex,color=ele2))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma)+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(y="Niche breadth",x="Elevation",title="Ferns and lycophytes",
       subtitle = "Tropical region")+
  geom_line(data=fnb,color="black",linewidth=1)+
  annotate("text",label = expression(beta ~"1" ~ "= 0.01" ~";" ~ beta ~"2 = -0.01"),
           x=100,y=70,size=5,hjust = 0)+
  annotate("text", label = "psi = 418",
           x=100,y=65,size=5,hjust=0)  
ndft

# Piecewise regressions
db1 <- df |> 
  dplyr::filter(region == "Extratropical")
fit <- lm(data=db1, convex~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","convex")


ndfe <- df |> 
  dplyr::filter(region == "Extratropical") |> 
  ggplot(aes(x=elevation,y=convex,color=ele2))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma)+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(y="Niche breadth",x="Elevation",title="Ferns and lycophytes",
       subtitle = "Extratropical region")+
  geom_line(data=fnb,color="black",linewidth=1)+
  annotate("text",label = expression(beta ~"1" ~ "= 0.007" ~";" ~ beta ~"2 = -0.008"),
           x=100,y=70,size=5,hjust = 0)+
  annotate("text", label = "psi = 521",
           x=100,y=65,size=5,hjust=0)  
ndfe

# Piecewise regressions
db1 <- df |> 
  dplyr::filter(region == "Latitudinal generalists")
fit <- lm(data=db1, convex~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","convex")

ndfb <- df |> 
  dplyr::filter(region == "Latitudinal generalists") |> 
  ggplot(aes(x=elevation,y=convex,color=ele2))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma)+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(y="Niche breadth",x="Elevation",title="Ferns and lycophytes",
       subtitle = "Tropical and Extratropical")+
  geom_line(data=fnb,color="black",linewidth=1)+
  annotate("text",label = expression(beta ~"1" ~ "= 0.02" ~";" ~ beta ~"2 = -0.02"),
           x=100,y=70,size=5,hjust = 0)+
  annotate("text", label = "psi = 392",
           x=100,y=65,size=5,hjust=0)  
ndfb

#####RANGE -------------------------------------------------------------------------
# Piecewise regression
db1 <- df |> 
  dplyr::filter(region == "Tropical")
fit <- lm(data=db1, area~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","area")

adft <- df |> 
  dplyr::filter(region == "Tropical") |> 
  ggplot(aes(x=elevation,y=area,color=ele2))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma, limits=c(0,15000000))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(x="Elevation (m)",y=expression("Range size"~(km^2)),
       title="Ferns and lycophytes",subtitle = "Tropical region")+
  geom_line(data=fnb,color="black",linewidth=1)+
  annotate("text",label = expression(beta ~"1" ~ "= -224" ~";" ~ beta ~"2 = 199"),
           x=100,y=15000000,size=5,hjust=0)+
  annotate("text", label = "psi = 838",
           x=100,y=13000000,size=5,hjust=0)
adft

# Piecewise regression
db1 <- df |> 
  dplyr::filter(region == "Extratropical")
fit <- lm(data=db1, area~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","area")

adfe <- df |> 
  dplyr::filter(region == "Extratropical") |> 
  ggplot(aes(x=elevation,y=area,color=ele2))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma, limits=c(0,15000000))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(x="Elevation (m)",y=expression("Range size"~(km^2)),
       title="Ferns and lycophytes",subtitle = "Extratropical region")+
  geom_line(data=fnb,color="black",linewidth=1)+
  annotate("text",label = expression(beta ~"1" ~ "= -125" ~";" ~ beta ~"2 = 192"),
           x=100,y=15000000,size=5,hjust=0)+
  annotate("text", label = "psi = 2116",
           x=100,y=13000000,size=5,hjust=0)
adfe


# Piecewise regression
db1 <- df |> 
  dplyr::filter(region == "Latitudinal generalists")
fit <- lm(data=db1, area~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","area")


adfb <- df |> 
  dplyr::filter(region == "Latitudinal generalists") |> 
  ggplot(aes(x=elevation,y=area,color=ele2))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma, limits=c(0,15000000))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(x="Elevation (m)",y=expression("Range size"~(km^2)),
       title="Ferns and lycophytes",subtitle = "Tropical and Extratropical")+
  geom_line(data=fnb,color="black",linewidth=1)+
  annotate("text",label = expression(beta ~"1" ~ "= 5371" ~";" ~ beta ~"2 = -5583"),
           x=100,y=15000000,size=5,hjust=0)+
  annotate("text", label = "psi = 255",
           x=100,y=13000000,size=5,hjust=0)
adfb

ndft + ndfe + ndfb
adft + adfe + adfb

#### modelling ------------------------------------------------------
m2f <- lmer(data=df, area10 ~ convex2:ele2:region+(1|region))
r.squaredGLMM(m2f)
m2fcoef <- tidy(m2f,"ran_coefs") |> print()
tidy(m2f)

df <- df |> mutate(fit.m = predict(m2f, re.form=NA),
                   fit.c = predict(m2f, re.form=NULL))

fp2 <- pairs(emtrends(m2f,c("ele2","region"),var="convex2"),pbkrtest.limit = 6767)|>
  tidy() |> 
  mutate(group = "Ferns")

m1fgEB <-df |> 
  dplyr::filter(region=="Both") |> 
  ggplot(aes(y=area10,x=convex2,color=ele2))+
  geom_point(color="gray70",alpha = 0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank())+
  geom_line(aes(y=fit.c,col=ele2),linewidth=1)+
  scale_color_manual(values =c("#56b4e9","#d55e00"))+
    labs(y=expression("log10 [Range size"~(km^2)~"]"),x="Niche breadth",title="Ferns",
       subtitle = "Both")+
  scale_y_continuous(limits=c(1,8))+
  scale_x_continuous(limits=c(1,8))+
  #legend
  annotate("segment", y = 2.5, yend = 2.5, x = 5.5, xend = 6.1, 
           color = "#56b4e9",linewidth = 1)+
  annotate("segment", y = 2, yend = 2, x = 5.5, xend = 6.1, 
           color = "#d55e00",linewidth = 1)+
  
  annotate("text", label="slope = 0.46", x = 7.2, y =2.5,size= 5)+
  annotate("text", label="slope = 0.43", x = 7.2, y =2,size=5)
#annotate("text",label = expression(conditional~R^2 == 0.44),
#         x=1,y=8,size=5,hjust = 0)
m1fgEB

m1fgET <-df |> 
  dplyr::filter(region=="Tropical") |> 
  ggplot(aes(y=area10,x=convex2,color=ele2))+
  geom_point(color="gray70",alpha = 0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank())+
  geom_line(aes(y=fit.c,col=ele2),linewidth=1)+
  scale_color_manual(values =c("#56b4e9","#d55e00"))+
    labs(y=expression("log10 [Range size"~(km^2)~"]"),x="Niche breadth",title="Ferns",
       subtitle = "Tropical")+
  scale_y_continuous(limits=c(1,8))+
  scale_x_continuous(limits=c(1,8))+
  #legend
  annotate("segment", y = 2.5, yend = 2.5, x = 5.5, xend = 6.1, 
           color = "#56b4e9",linewidth = 1)+
  annotate("segment", y = 2, yend = 2, x = 5.5, xend = 6.1, 
           color = "#d55e00",linewidth = 1)+
  
  annotate("text", label="slope = 0.80", x = 7.2, y =2.5,size= 5)+
  annotate("text", label="slope = 0.71", x = 7.2, y =2,size=5)
#annotate("text",label = expression(conditional~R^2 == 0.44),
#         x=1,y=8,size=5,hjust = 0)
m1fgET

m1fgEE <-df |> 
  dplyr::filter(region=="Extratropical") |> 
  ggplot(aes(y=area10,x=convex2,color=ele2))+
  geom_point(color="gray70",alpha = 0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank())+
  geom_line(aes(y=fit.c,col=ele2),linewidth=1)+
  scale_color_manual(values =c("#56b4e9","#d55e00"))+
    labs(y=expression("log10 [Range size"~(km^2)~"]"),x="Niche breadth",title="Ferns",
       subtitle = "Extratropical")+
  scale_y_continuous(limits=c(1,8))+
  scale_x_continuous(limits=c(1,8))+
  #legend
  annotate("segment", y = 2.5, yend = 2.5, x = 5.5, xend = 6.1, 
           color = "#56b4e9",linewidth = 1)+
  annotate("segment", y = 2, yend = 2, x = 5.5, xend = 6.1, 
           color = "#d55e00",linewidth = 1)+
  
  annotate("text", label="slope = 0.56", x = 7.2, y =2.5,size= 5)+
  annotate("text", label="slope = 0.46", x = 7.2, y =2,size=5)
#annotate("text",label = expression(conditional~R^2 == 0.44),
#         x=1,y=8,size=5,hjust = 0)
m1fgEE


dfl <- df |> dplyr::filter(ele2 == "Low elevation") |> 
  dplyr::select(convex2) |> 
  pull()
dfm <- df |> dplyr::filter(ele2 == "Mountain") |> 
  dplyr::select(convex2)|> 
  pull()
ks.test(dfl,dfm)

vtfe1 <- df |>
  ggplot(aes(x=ele2,y=convex2,fill=ele2))+
  geom_flat_violin(position = position_nudge(x = .12, y = 0), alpha = 0.6)+
  geom_jitter(aes(color = ele2),
              width=0.08, alpha = 0.6)+
  # facet_wrap(.~band)+
  #scale_y_continuous(limits=c(1,8))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none",
        #axis.text.y = element_blank(),
        strip.background = element_rect(colour="white", fill="white"))+
  #labs(y=expression("log10 [Range size"~(km^2)~"]"), x="")+
  labs(y="Niche breadth", x="",title="Ferns")+
  scale_color_grey()+
  scale_fill_grey()+
  scale_y_continuous(limits=c(1,8))+
  coord_flip()+
  annotate("text", label="p-value = 0.02", x = 0.7, y =6,size=5)
vtfe1 

dfl <- df |> dplyr::filter(ele2 == "Low elevation") |> 
  dplyr::select(area10) |> 
  pull()
dfm <- df |> dplyr::filter(ele2 == "Mountain") |> 
  dplyr::select(area10)|> 
  pull()
ks.test(dfl,dfm)

vtfe2 <- df |>
  ggplot(aes(x=ele2,y=area10,fill=ele2))+
  geom_flat_violin(position = position_nudge(x = .12, y = 0), alpha = 0.6)+
  geom_jitter(aes(color = ele2),
              width=0.08, alpha = 0.6)+
  # facet_wrap(.~band)+
  #scale_y_continuous(limits=c(1,8))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_rect(colour="white", fill="white"))+
  #labs(y=expression("log10 [Range size"~(km^2)~"]"), x="")+
  labs(y=expression("log10 [Range size"~(km^2)~"]"), x="",title = "Ferns")+
  scale_color_grey()+
  scale_fill_grey()+
  scale_y_continuous(limits=c(2,8))+
  coord_flip()+
  annotate("text", label="p-value < 0.001", x = 0.7, y =6,size=5)
vtfe2 


### Gymnosperms ---------------------------------------------------------------------
dg <- d |> 
  dplyr::filter(higher_plant_group=="Gymnosperms")

#### distribution graph and analysis --------------------------------------------
#####NICHE -------------------------------------------------------------------------
# Piecewise regressions
db1 <- dg |> 
  dplyr::filter(region == "Tropical")
fit <- lm(data=db1, convex~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","convex")


ndgt <- dg |> 
  ggplot(aes(x=elevation,y=convex,color=ele2))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma)+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(y="Niche breadth",x="Elevation",title="Gymnosperms",
       subtitle = "Tropical region")+
  geom_line(data=fnb,color="black",linewidth=1)+
  annotate("text",label = expression(beta ~"1" ~ "= 0.06" ~";" ~ beta ~"2 = -0.06"),
           x=100,y=70,size=5,hjust = 0)+
  annotate("text", label = "psi = 108",
           x=100,y=65,size=5,hjust=0)  
ndgt

# Piecewise regressions
db1 <- dg |> 
  dplyr::filter(region == "Extratropical")
fit <- lm(data=db1, convex~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","convex")


ndge <- dg |> 
  dplyr::filter(region == "Extratropical") |> 
  ggplot(aes(x=elevation,y=convex,color=ele2))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma)+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(y="Niche breadth",x="Elevation",title="Gymnosperms",
       subtitle = "Extratropical region")+
  geom_line(data=fnb,color="black",linewidth=1)+
  annotate("text",label = expression(beta ~"1" ~ "= 0.008" ~";" ~ beta ~"2 = -0.008"),
           x=100,y=70,size=5,hjust = 0)+
  annotate("text", label = "psi = 454",
           x=100,y=65,size=5,hjust=0)  
ndge


# Piecewise regressions
db1 <- dg |> 
  dplyr::filter(region == "Latitudinal generalists")
fit <- lm(data=db1, convex~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","convex")

ndgb <- dg |> 
  dplyr::filter(region == "Latitudinal generalists") |> 
  ggplot(aes(x=elevation,y=convex,color=ele2))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma)+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(y="Niche breadth",x="Elevation",title="Gymnosperms",
       subtitle = "Tropical and Extratropical")+
  geom_line(data=fnb,color="black",linewidth=1)+
  annotate("text",label = expression(beta ~"1" ~ "= 0.01" ~";" ~ beta ~"2 = -0.01"),
           x=100,y=70,size=5,hjust = 0)+
  annotate("text", label = "psi = 440",
           x=100,y=65,size=5,hjust=0)  
ndgb

#####RANGE -------------------------------------------------------------------------
# Piecewise regression
db1 <- dg |> 
  dplyr::filter(region == "Tropical")
fit <- lm(data=db1, area~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","area")

adgt <- dg |> 
  dplyr::filter(region == "Tropical") |> 
  ggplot(aes(x=elevation,y=area,color=ele2))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma, limits=c(0,10000000))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(x="Elevation (m)",y=expression("Range size"~(km^2)),
       title="Gymnosperms",subtitle = "Tropical region")+
  geom_line(data=fnb,color="black",linewidth=1)+
  annotate("text",label = expression(beta ~"1" ~ "= -24275" ~";" ~ beta ~"2 = 24269"),
           x=100,y=10000000,size=5,hjust=0)+
  annotate("text", label = "psi = 106",
           x=100,y=9000000,size=5,hjust=0)
adgt

# Piecewise regression
db1 <- dg |> 
  dplyr::filter(region == "Extratropical")
fit <- lm(data=db1, area~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","area")

adge <- dg |> 
  dplyr::filter(region == "Extratropical") |> 
  ggplot(aes(x=elevation,y=area,color=ele2))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma, limits=c(0,10000000))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(x="Elevation (m)",y=expression("Range size"~(km^2)),
       title="Gymnosperms",subtitle = "Extratropical region")+
  geom_line(data=fnb,color="black",linewidth=1)+
  annotate("text",label = expression(beta ~"1" ~ "= 2033" ~";" ~ beta ~"2 = -2199"),
           x=100,y=10000000,size=5,hjust=0)+
  annotate("text", label = "psi = 402",
           x=100,y=9000000,size=5,hjust=0)
adge

# Piecewise regression
db1 <- dg |> 
  dplyr::filter(region == "Latitudinal generalists")
fit <- lm(data=db1, area~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","area")


adgb <- dg |> 
  dplyr::filter(region == "Latitudinal generalists") |> 
  ggplot(aes(x=elevation,y=area,color=ele2))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma, limits=c(0,10000000))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(x="Elevation (m)",y=expression("Range size"~(km^2)),
       title="Gymnosperms",subtitle = "Tropical and Extratropical")+
  geom_line(data=fnb,color="black",linewidth=1)+
  annotate("text",label = expression(beta ~"1" ~ "= 880" ~";" ~ beta ~"2 = -999"),
           x=100,y=10000000,size=5,hjust=0)+
  annotate("text", label = "psi = 707",
           x=100,y=9000000,size=5,hjust=0)
adgb

ndgt + ndge + ndgb
adgt + adge + adgb

#### modelling -------------------------------------------------------------
m2g <- lmer(data=dg, area10 ~ convex2:ele2:region+(1|region))
r.squaredGLMM(m2g)
m2gcoef <- tidy(m2g,"ran_coefs") |> print()
tidy(m2g)

dg <- dg |> mutate(fit.m = predict(m2g, re.form=NA),
                   fit.c = predict(m2g, re.form=NULL))

gp2 <- pairs(emtrends(m2g,c("ele2","region"),var="convex2"),pbkrtest.limit = 6767)|>
  tidy() |> 
  mutate(group = "Gymnosperms")

m1ggEB <-dg |> 
  dplyr::filter(region=="Both") |> 
  ggplot(aes(y=area10,x=convex2,color=ele2))+
  geom_point(color="gray70",alpha = 0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank())+
  geom_line(aes(y=fit.c,col=ele2),linewidth=1)+
  scale_color_manual(values =c("#56b4e9","#d55e00"))+
  #geom_abline(intercept = 3.68, slope = 0.569,color="#56b4e9", linewidth=1)+ #Both
  #geom_abline(intercept = 2.69, slope = 0.844,color="#d55e00", linewidth=1)+ #Extratropical
  #geom_abline(intercept = 2.72, slope = 0.675,color="#009e73", linewidth=1)+ #Tropical
  labs(y=expression("log10 [Range size"~(km^2)~"]"),x="Niche breadth",title="Gymnosperms",
       subtitle = "Both")+
  scale_y_continuous(limits=c(1,8))+
  scale_x_continuous(limits=c(1,8))+
  #legend
  annotate("segment", y = 2.5, yend = 2.5, x = 5.5, xend = 6.1, 
           color = "#56b4e9",linewidth = 1)+
  annotate("segment", y = 2, yend = 2, x = 5.5, xend = 6.1, 
           color = "#d55e00",linewidth = 1)+
  
  annotate("text", label="slope = 0.37", x = 7.2, y =2.5,size= 5)+
  annotate("text", label="slope = 0.43", x = 7.2, y =2,size=5)
#annotate("text",label = expression(conditional~R^2 == 0.44),
#         x=1,y=8,size=5,hjust = 0)
m1ggEB


m1ggET <-dg |> 
  dplyr::filter(region=="Tropical") |> 
  ggplot(aes(y=area10,x=convex2,color=ele2))+
  geom_point(color="gray70",alpha = 0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank())+
  geom_line(aes(y=fit.c,col=ele2),linewidth=1)+
  scale_color_manual(values =c("#56b4e9","#d55e00"))+
  #geom_abline(intercept = 3.68, slope = 0.569,color="#56b4e9", linewidth=1)+ #Both
  #geom_abline(intercept = 2.69, slope = 0.844,color="#d55e00", linewidth=1)+ #Extratropical
  #geom_abline(intercept = 2.72, slope = 0.675,color="#009e73", linewidth=1)+ #Tropical
  labs(y=expression("log10 [Range size"~(km^2)~"]"),x="Niche breadth",title="Gymnosperms",
       subtitle = "Tropical")+
  scale_y_continuous(limits=c(1,8))+
  scale_x_continuous(limits=c(1,8))+
  #legend
  annotate("segment", y = 2.5, yend = 2.5, x = 5.5, xend = 6.1, 
           color = "#56b4e9",linewidth = 1)+
  annotate("segment", y = 2, yend = 2, x = 5.5, xend = 6.1, 
           color = "#d55e00",linewidth = 1)+
  
  annotate("text", label="slope = 0.69", x = 7.2, y =2.5,size= 5)+
  annotate("text", label="slope = 0.67", x = 7.2, y =2,size=5)
#annotate("text",label = expression(conditional~R^2 == 0.44),
#         x=1,y=8,size=5,hjust = 0)
m1ggET


m1ggEE <-dg |> 
  dplyr::filter(region=="Extratropical") |> 
  ggplot(aes(y=area10,x=convex2,color=ele2))+
  geom_point(color="gray70",alpha = 0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank())+
  geom_line(aes(y=fit.c,col=ele2),linewidth=1)+
  scale_color_manual(values =c("#56b4e9","#d55e00"))+
  #geom_abline(intercept = 3.68, slope = 0.569,color="#56b4e9", linewidth=1)+ #Both
  #geom_abline(intercept = 2.69, slope = 0.844,color="#d55e00", linewidth=1)+ #Extratropical
  #geom_abline(intercept = 2.72, slope = 0.675,color="#009e73", linewidth=1)+ #Tropical
  labs(y=expression("log10 [Range size"~(km^2)~"]"),x="Niche breadth",title="Gymnosperms",
       subtitle = "Extratropical")+
  scale_y_continuous(limits=c(1,8))+
  scale_x_continuous(limits=c(1,8))+
  #legend
  annotate("segment", y = 2.5, yend = 2.5, x = 5.5, xend = 6.1, 
           color = "#56b4e9",linewidth = 1)+
  annotate("segment", y = 2, yend = 2, x = 5.5, xend = 6.1, 
           color = "#d55e00",linewidth = 1)+
  
  annotate("text", label="slope = 0.63", x = 7.2, y =2.5,size= 5)+
  annotate("text", label="slope = 0.57", x = 7.2, y =2,size=5)
#annotate("text",label = expression(conditional~R^2 == 0.44),
#         x=1,y=8,size=5,hjust = 0)
m1ggEE

dgl <- dg |> dplyr::filter(ele2 == "Low elevation") |> 
  dplyr::select(convex2) |> 
  pull()
dgm <- dg |> dplyr::filter(ele2 == "Mountain") |> 
  dplyr::select(convex2)|> 
  pull()
ks.test(dgl,dgm)

vtge1 <- dg |>
  ggplot(aes(x=ele2,y=convex2,fill=ele2))+
  geom_flat_violin(position = position_nudge(x = .12, y = 0), alpha = 0.6)+
  geom_jitter(aes(color = ele2),
              width=0.08, alpha = 0.6)+
  # facet_wrap(.~band)+
  #scale_y_continuous(limits=c(1,8))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "bottom",
        #axis.text.y = element_blank(),
        strip.background = element_rect(colour="white", fill="white"))+
  #labs(y=expression("log10 [Range size"~(km^2)~"]"), x="")+
  labs(y="Niche breadth", x="",title="Gymnosperms")+
  scale_color_grey()+
  scale_fill_grey()+
  scale_y_continuous(limits=c(1,8))+
  coord_flip()+
  annotate("text", label="p-value = 0.02", x = 0.7, y =6,size=5)
vtge1 

dgl <- dg |> dplyr::filter(ele2 == "Low elevation") |> 
  dplyr::select(area10) |> 
  pull()
dgm <- dg |> dplyr::filter(ele2 == "Mountain") |> 
  dplyr::select(area10)|> 
  pull()
ks.test(dgl,dgm)

vtge2 <- dg |>
  ggplot(aes(x=ele2,y=area10,fill=ele2))+
  geom_flat_violin(position = position_nudge(x = .12, y = 0), alpha = 0.6)+
  geom_jitter(aes(color = ele2),
              width=0.08, alpha = 0.6)+
  # facet_wrap(.~band)+
  #scale_y_continuous(limits=c(1,8))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none",
       axis.text.y = element_blank(),
        strip.background = element_rect(colour="white", fill="white"))+
  #labs(y=expression("log10 [Range size"~(km^2)~"]"), x="")+
  labs(y=expression("log10 [Range size"~(km^2)~"]"), x="",title="Gymnosperms")+
  scale_color_grey()+
  scale_fill_grey()+
  scale_y_continuous(limits=c(2,8))+
  coord_flip()+
  annotate("text", label="p-value = 0.10", x = 0.7, y =6,size=5)
vtge2 


### Flowering plants -------------------------------------------------
dff <- d |> 
  dplyr::filter(higher_plant_group=="Flowering plants")

#### distribution graph and analysis --------------------------------------------
#####NICHE -------------------------------------------------------------------------
# Piecewise regressions
db1 <- dff |> 
  dplyr::filter(region == "Tropical")
fit <- lm(data=db1, convex~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","convex")


ndfft <- dff |> 
  dplyr::filter(region == "Tropical") |> 
  ggplot(aes(x=elevation,y=convex,color=ele2))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma)+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(y="Niche breadth",x="Elevation",title="Flowering plants",
       subtitle = "Tropical region")+
  geom_line(data=fnb,color="black",linewidth=1)+
  annotate("text",label = expression(beta ~"1" ~ "= 0.02" ~";" ~ beta ~"2 = -0.02"),
           x=100,y=70,size=5,hjust = 0)+
  annotate("text", label = "psi = 193",
           x=100,y=65,size=5,hjust=0) 
ndfft

#Piecewise regressions
db1 <- dff |> 
  dplyr::filter(region == "Extratropical")
fit <- lm(data=db1, convex~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","convex")

ndffe <- dff |> 
  dplyr::filter(region == "Extratropical") |> 
  ggplot(aes(x=elevation,y=convex,color=ele2))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma)+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(y="Niche breadth",x="Elevation",title="Flowering plants",
       subtitle = "Extratropical region")+
  geom_line(data=fnb,color="black",linewidth=1)+
  annotate("text",label = expression(beta ~"1" ~ "= 0.009" ~";" ~ beta ~"2 = -0.009"),
           x=100,y=70,size=5,hjust = 0)+
  annotate("text", label = "psi = 357",
           x=100,y=65,size=5,hjust=0) 
ndffe

#Piecewise regressions
db1 <- dff |> 
  dplyr::filter(region == "Latitudinal generalists")
fit <- lm(data=db1, convex~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","convex")

ndffb <- dff |> 
  dplyr::filter(region == "Latitudinal generalists") |> 
  ggplot(aes(x=elevation,y=convex,color=ele2))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma)+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(y="Niche breadth",x="Elevation",title="Flowering plants",
       subtitle = "Tropical and Extratropical")+
  geom_line(data=fnb,color="black",linewidth=1)+
  annotate("text",label = expression(beta ~"1" ~ "= 0.02" ~";" ~ beta ~"2 = -0.02"),
           x=100,y=70,size=5,hjust = 0)+
  annotate("text", label = "psi = 266",
           x=100,y=65,size=5,hjust=0) 
ndffb

#####RANGE -------------------------------------------------------------------------
# Piecewise regression
db1 <- dff |> 
  dplyr::filter(region == "Tropical")
fit <- lm(data=db1, area~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","area")

adfft <- dff |> 
  dplyr::filter(region == "Tropical") |> 
  ggplot(aes(x=elevation,y=area,color=ele2))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma, limits=c(0,20000000))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(x="Elevation (m)",y=expression("Range size"~(km^2)),
       title="Flowering plants",subtitle = "Tropical region")+
  geom_line(data=fnb,color="black",linewidth=1)+
  annotate("text",label = expression(beta ~"1" ~ "= -254" ~";" ~ beta ~"2 = 245"),
           x=100,y=20000000,size=5,hjust=0)+
  annotate("text", label = "psi = 1095",
           x=100,y=18000000,size=5,hjust=0)
adfft

# Piecewise regression
db1 <- dff |> 
  dplyr::filter(region == "Extratropical")
fit <- lm(data=db1, area~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","area")


adffe <- dff |> 
  dplyr::filter(region == "Extratropical") |> 
  ggplot(aes(x=elevation,y=area,color=ele2))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma, limits=c(0,20000000))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(x="Elevation (m)",y=expression("Range size"~(km^2)),
       title="Flowering plants",subtitle = "Extratropical region")+
  geom_line(data=fnb,color="black",linewidth=1)+
  annotate("text",label = expression(beta ~"1" ~ "= 1367" ~";" ~ beta ~"2 = -1404"),
           x=100,y=20000000,size=5,hjust=0)+
  annotate("text", label = "psi = 225",
           x=100,y=18000000,size=5,hjust=0)
adffe


# Piecewise regression
db1 <- dff |> 
  dplyr::filter(region == "Latitudinal generalists")
fit <- lm(data=db1, area~elevation )
summary(fit)
s <- segmented(fit)
print(s)
confint(s, "elevation")

p <- as.data.frame(predict(s,interval="confidence"))
fnb <-  as.data.frame(cbind(db1$elevation, s$fitted.values))
colnames(fnb) <-  c("elevation","area")

adffb <- dff |> 
  dplyr::filter(region == "Latitudinal generalists") |> 
  ggplot(aes(x=elevation,y=area,color=ele2))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        rect = element_rect(fill = "transparent",color="transparent"))+
  scale_x_continuous(limits = c(0,5000),breaks = c(0,1000,2000,3000,4000,5000))+
  scale_y_continuous(labels=comma, limits=c(0,20000000))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  labs(x="Elevation (m)",y=expression("Range size"~(km^2)),
       title="Flowering plants",subtitle = "Tropical and Extratropical")+
  geom_line(data=fnb,color="black",linewidth=1)+
  annotate("text",label = expression(beta ~"1" ~ "= 4663" ~";" ~ beta ~"2 = -5217"),
           x=100,y=20000000,size=5,hjust=0)+
  annotate("text", label = "psi = 267",
           x=100,y=18000000,size=5,hjust=0)
adffb

ndfft + ndffe + ndffb
adfft + adffe + adffb


#### modelling -------------------------------------------------------
m2ff <- lmer(data=dff, area10 ~ convex2:ele2:region+(1|region) )
r.squaredGLMM(m2ff)
m2ffcoef <- tidy(m2ff,"ran_coefs") |> print()
tidy(m2ff)

dff <- dff |> mutate(fit.m = predict(m2ff, re.form=NA),
                   fit.c = predict(m2ff, re.form=NULL))

ffp2 <- pairs(emtrends(m2ff,c("ele2","region"),var="convex2"),pbkrtest.limit = 141939)|>
  tidy() |> 
  mutate(group = "Flowering plants")



m1ffgEB <-dff |> 
  dplyr::filter(region=="Both") |> 
  ggplot(aes(y=area10,x=convex2,color=ele2))+
  geom_point(color="gray70",alpha = 0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank())+
  geom_line(aes(y=fit.c,col=ele2),linewidth=1)+
  scale_color_manual(values =c("#56b4e9","#d55e00"))+
  #geom_abline(intercept = 3.68, slope = 0.569,color="#56b4e9", linewidth=1)+ #Both
  #geom_abline(intercept = 2.69, slope = 0.844,color="#d55e00", linewidth=1)+ #Extratropical
  #geom_abline(intercept = 2.72, slope = 0.675,color="#009e73", linewidth=1)+ #Tropical
  labs(y=expression("log10 [Range size"~(km^2)~"]"),x="Niche breadth",title="Flowering plants",
       subtitle = "Both")+
  scale_y_continuous(limits=c(1,8))+
  scale_x_continuous(limits=c(1,8))+
  #legend
  annotate("segment", y = 2.5, yend = 2.5, x = 5.5, xend = 6.1, 
           color = "#56b4e9",linewidth = 1)+
  annotate("segment", y = 2, yend = 2, x = 5.5, xend = 6.1, 
           color = "#d55e00",linewidth = 1)+
  
  annotate("text", label="slope = 0.41", x = 7.2, y =2.5,size= 5)+
  annotate("text", label="slope = 0.34", x = 7.2, y =2,size=5)
#annotate("text",label = expression(conditional~R^2 == 0.44),
#         x=1,y=8,size=5,hjust = 0)
m1ffgEB


m1ffgET <-dff |> 
  dplyr::filter(region=="Tropical") |> 
  ggplot(aes(y=area10,x=convex2,color=ele2))+
  geom_point(color="gray70",alpha = 0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank())+
  geom_line(aes(y=fit.c,col=ele2),linewidth=1)+
  scale_color_manual(values =c("#56b4e9","#d55e00"))+
  #geom_abline(intercept = 3.68, slope = 0.569,color="#56b4e9", linewidth=1)+ #Both
  #geom_abline(intercept = 2.69, slope = 0.844,color="#d55e00", linewidth=1)+ #Extratropical
  #geom_abline(intercept = 2.72, slope = 0.675,color="#009e73", linewidth=1)+ #Tropical
  labs(y=expression("log10 [Range size"~(km^2)~"]"),x="Niche breadth",title="Flowering plants",
       subtitle = "Tropical")+
  scale_y_continuous(limits=c(1,8))+
  scale_x_continuous(limits=c(1,8))+
  #legend
  annotate("segment", y = 2.5, yend = 2.5, x = 5.5, xend = 6.1, 
           color = "#56b4e9",linewidth = 1)+
  annotate("segment", y = 2, yend = 2, x = 5.5, xend = 6.1, 
           color = "#d55e00",linewidth = 1)+
  
  annotate("text", label="slope = 0.74", x = 7.2, y =2.5,size= 5)+
  annotate("text", label="slope = 0.61", x = 7.2, y =2,size=5)
#annotate("text",label = expression(conditional~R^2 == 0.44),
#         x=1,y=8,size=5,hjust = 0)
m1ffgET


m1ffgEE <-dff |> 
  dplyr::filter(region=="Extratropical") |> 
  ggplot(aes(y=area10,x=convex2,color=ele2))+
  geom_point(color="gray70",alpha = 0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank())+
  geom_line(aes(y=fit.c,col=ele2),linewidth=1)+
  scale_color_manual(values =c("#56b4e9","#d55e00"))+
  #geom_abline(intercept = 3.68, slope = 0.569,color="#56b4e9", linewidth=1)+ #Both
  #geom_abline(intercept = 2.69, slope = 0.844,color="#d55e00", linewidth=1)+ #Extratropical
  #geom_abline(intercept = 2.72, slope = 0.675,color="#009e73", linewidth=1)+ #Tropical
  labs(y=expression("log10 [Range size"~(km^2)~"]"),x="Niche breadth",title="Flowering plants",
       subtitle = "Extratropical")+
  scale_y_continuous(limits=c(1,8))+
  scale_x_continuous(limits=c(1,8))+
  #legend
  annotate("segment", y = 2.5, yend = 2.5, x = 5.5, xend = 6.1, 
           color = "#56b4e9",linewidth = 1)+
  annotate("segment", y = 2, yend = 2, x = 5.5, xend = 6.1, 
           color = "#d55e00",linewidth = 1)+
  
  annotate("text", label="slope = 0.71", x = 7.2, y =2.5,size= 5)+
  annotate("text", label="slope = 0.61", x = 7.2, y =2,size=5)
#annotate("text",label = expression(conditional~R^2 == 0.44),
#         x=1,y=8,size=5,hjust = 0)
m1ffgEE

dffl <- dff |> dplyr::filter(ele2 == "Low elevation") |> 
  dplyr::select(convex2) |> 
  pull()
dffm <- dff|> dplyr::filter(ele2 == "Mountain") |> 
  dplyr::select(convex2)|> 
  pull()
ks.test(dffl,dffm)

vtffe1 <- dff |>
  ggplot(aes(x=ele2,y=convex2,fill=ele2))+
  geom_flat_violin(position = position_nudge(x = .12, y = 0), alpha = 0.6)+
  geom_jitter(aes(color = ele2),
              width=0.08, alpha = 0.6)+
  # facet_wrap(.~band)+
  #scale_y_continuous(limits=c(1,8))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none",
        #axis.text.y = element_blank(),
        strip.background = element_rect(colour="white", fill="white"))+
  #labs(y=expression("log10 [Range size"~(km^2)~"]"), x="")+
  labs(y="Niche breadth", x="",title="Flowering plants")+
  scale_color_grey()+
  scale_fill_grey()+
  scale_y_continuous(limits=c(1,8))+
  coord_flip()+
  annotate("text", label="p-value < 0.001", x = 0.7, y =6,size=5)
vtffe1 

dffl <- dff |> dplyr::filter(ele2 == "Low elevation") |> 
  dplyr::select(area10) |> 
  pull()
dffm <- dff|> dplyr::filter(ele2 == "Mountain") |> 
  dplyr::select(area10)|> 
  pull()
ks.test(dffl,dffm)

vtffe2 <- dff |>
  ggplot(aes(x=ele2,y=area10,fill=ele2))+
  geom_flat_violin(position = position_nudge(x = .12, y = 0), alpha = 0.6)+
  geom_jitter(aes(color = ele2),
              width=0.08, alpha = 0.6)+
  # facet_wrap(.~band)+
  #scale_y_continuous(limits=c(1,8))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none",
       axis.text.y = element_blank(),
        strip.background = element_rect(colour="white", fill="white"))+
  #labs(y=expression("log10 [Range size"~(km^2)~"]"), x="")+
  labs(y=expression("log10 [Range size"~(km^2)~"]"), x="",title="Flowering plants")+
  scale_color_grey()+
  scale_fill_grey()+
  coord_flip()+
  scale_y_continuous(limits=c(2,8))+
  annotate("text", label="p-value < 0.001", x = 0.7, y =6,size=5)
vtffe2 


### Save composite plots and tables ----------------------------
ncp = ndbt + ndbe + ndbb+
  ndft + ndfe + ndfb+
  ndgt + ndge + ndgb+
  ndfft + ndffe + ndffb+
plot_layout(guides = 'collect',ncol = 3)+
plot_annotation(tag_levels = "a",tag_suffix = ")")&
theme(legend.position="none")
ggsave(here("products","RESU_dist_ele_niche_v2.jpeg"),
       units="cm",dpi = 300,height = 50,width = 50,plot=ncp)

rcp = adbt + adbe + adbb+
  adft + adfe + adfb+
  adgt + adge + adgb+
  adfft + adffe + adffb+
plot_layout(guides = 'collect',ncol = 3)+
  plot_annotation(tag_levels = "a",tag_suffix = ")")&
  theme(legend.position="none")
ggsave(here("products","RESU_dist_ele_range_v2.jpeg"),
       units="cm",dpi = 300,height = 50,width = 50,plot=rcp)

  

cp=m1bgEl + m1bgEm +
m1fgEl + m1fgEm +
m1ffgEl + m1ffgEm +
m1ggEl + m1ggEm +
  plot_layout(guides = 'collect',ncol = 2)+
  plot_annotation(tag_levels = "a",tag_suffix = ")")&
  theme(legend.position="bottom")
ggsave(here("products","RESU_model2_v5.jpeg"),
       units="cm",dpi = 300,
       width=28, height = 50,
       plot=cp)

cp=m1bgET + m1bgEE + m1bgEB+
m1fgET + m1fgEE + m1fgEB+
m1ffgET + m1ffgEE + m1ffgEB+
m1ggET + m1ggEE + m1ggEB+
  plot_layout(guides = 'collect',ncol = 3)+
  plot_annotation(tag_levels = "a",tag_suffix = ")")&
  theme(legend.position="bottom")
ggsave(here("products","RESU_model2_v6.jpeg"),
       units="cm",dpi = 300,
       width=38, height = 45,
       plot=cp)

# make the results table

bp3 <- bp2 |> slice(c(1,10,15))
fp3 <- fp2 |> slice(c(1,10,15))
gp3 <- gp2 |> slice(c(1,10,15))
ffp3 <- ffp2 |> slice(c(1,10,15))

ctm2 <- rbind(bp3,fp3,gp3,ffp3) |> 
  dplyr::mutate(p.value = round(adj.p.value,3)) |> 
  dplyr::select(contrast,p.value, group) |> 
  dplyr::relocate(group,.before = contrast) |> 
  dplyr::rename(Interaction = contrast) |> 
  flextable::flextable() %>%                     
  flextable::autofit() %>%    
  flextable::save_as_docx(path = here("products","Table_model2.docx"))

# make the distribution figure

cp=vtbe1 + vtbe2 +
vtfe1 + vtfe2 + 
vtge1 + vtge2 +
vtffe1 + vtffe2 +
  plot_layout(guides = 'collect',ncol = 2,
              axis_titles = "collect")+
  plot_annotation(tag_levels = "a",tag_suffix = ")")&
  theme(legend.position="bottom")
  ggsave(here("products","RESU_model2_v5_dist.jpeg"),
         units="cm",dpi = 300,
         width=30, height = 50,
         plot=cp)

## 3) model 3: dominance -----------------------------------------------
pal3 <- c("#cc79a7","#f0e442")

# Load Hyperdominant species
data_hy <- dataGlobal

### 3.1 Distribution graphs------------------------------------------------------------
###All species ------------------------------------------------------------------------
#### NICHE ---------------------------------------------------------------------------
ndtH <- data_hy |> 
  dplyr::filter(area>100) |> 
  dplyr::filter(convex>0) |> 
  ggplot(aes(x=convex2,y=d))+
  geom_density_ridges(aes(fill=source,scale=5,color=paste(h,nd)),linewidth=1,alpha=0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.grid = element_blank())+
  annotate("segment", y = 2.7, yend = 2.7, x = 5.5, xend = 6.1, 
           color = "#f0e442",linewidth = 1)+
  annotate("text", label="Non-dominants", x = 8, y =2.7,size= 5)+
  
  annotate("segment", y = 2.5, yend = 2.5, x = 5.5, xend = 6.1, 
           color ="#cc79a7",linewidth = 1)+
  annotate("text", label="Hyperdominants", x = 8, y =2.5,size=5)+
  scale_color_manual(values=pal3,expand=c(0,0),
                     label=c("Hyperdominants","Non-dominants"),
                     guide="legend")+
  scale_fill_manual(values=pal3,label=c("Hyperdominants","Non-dominants"),
                    guide="legend")+
  scale_y_discrete(expand=c(0,0))+
  scale_x_continuous(limits=c(0,10),breaks=c(0,2,4,6,8,10))+
  labs(x="sqrt (Niche breadth)")+
  annotate("text",label=
             "ks ***
kw ***", 
           x=8, y=1.5  ,size=5,hjust=0)
ndtH

xx <- data_hy |> filter(source=="Hyperdominants")
summary(lm(log10(area)~log10(n),data=xx))
xx1 <- data_hy |> filter(source=="Non-dominants")
summary(lm(log10(area)~log10(n),data=xx1))

occplotNB<- data_hy |> 
  dplyr::filter(area>100) |> 
  dplyr::filter(convex>0) |> 
  ggplot(aes(y=convex2,x=n,color=source))+
  geom_point(alpha = 0.6,color="grey60", size=1)+
  #facet_wrap(.~higher_plant_group)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=10),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"))+
  scale_x_log10(labels=comma)+
  scale_y_continuous(limits=c(0,8))+
  # scale_y_log10()+
  geom_smooth(method="lm")+
  scale_color_manual(values = pal3)+
  labs(x=expression(log["10"] ~ "(Occurrences)"),y="sqrt (Niche breadth)")+
  annotate("text", label = "R^2 == 0.27 ~ '***'",
           parse = TRUE,
             #expression(R^2 == 0.27 ~  "***" ),
           x=10000,y=7.9,size=4,hjust = 0)+
  annotate("text", label = "R^2 == 0.31 ~ '***'",
           parse = TRUE,
             #expression(R^2 == 0.31 ~ "***" ),
           x=10000,y=6.9,size=4,hjust = 0)+
  annotate("segment", y = 7.9, yend = 7.9, x = 8000, xend = 5000, 
         color = "#cc79a7",linewidth = 1)+
  annotate("segment", y = 6.9, yend = 6.9, x = 8000, xend = 5000, 
           color = "#f0e442",linewidth = 1)
occplotNB

#### RANGE ---------------------------------------------------------------------------
adtH <- data_hy |> 
  dplyr::filter(area>100) |> 
  dplyr::filter(convex>0) |> 
  ggplot(aes(x=area10,y=d))+
  geom_density_ridges(aes(fill=source,scale=5,color=paste(h,nd)),linewidth=1,alpha=0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.grid = element_blank())+
  scale_color_manual(values=pal3)+
  scale_fill_manual(values=pal3)+
  scale_y_discrete(expand=c(0,0))+
  scale_x_continuous(limits=c(2,10),breaks=c(2,4,6,8,10))+
  labs(x=expression("log10 [Range size"~(km^2)~"]"))+
  annotate("text",label=
             "ks ***
kw ***", 
           x=8, y=2  ,size=5,hjust=0)+
  annotate("segment", y = 4.3, yend = 4.3, x = 2, xend = 2.6, 
           color = "#f0e442",linewidth = 1)+
  annotate("text", label="Non-dominants", x = 4, y =4.3,size= 5)+
  annotate("segment", y = 3.9, yend = 3.9, x = 2, xend = 2.6, 
           color ="#cc79a7",linewidth = 1)+
  annotate("text", label="Hyperdominants", x = 4, y =3.9,size=5)
adtH

xx <- data_hy |> filter(source=="Hyperdominants")
summary(lm(sqrt(convex)~log10(n),data=xx))
xx1 <- data_hy |> filter(source=="Non-dominants")
summary(lm(sqrt(convex)~log10(n),data=xx1))

occplotRS<- data_hy |> 
  dplyr::filter(area>100) |> 
  dplyr::filter(convex>0) |> 
  ggplot(aes(y=area10,x=n,color=source))+
  geom_point(alpha = 0.6,color="grey60", size=1)+
  #facet_wrap(.~higher_plant_group)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=10),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"))+
  scale_x_log10(labels=comma)+
  scale_y_continuous(limits=c(2,8))+
  # scale_y_log10()+
  geom_smooth(method="lm")+
  scale_color_manual(values = pal3)+
  labs(x=expression(log["10"] ~ "(Occurrences)"),y=expression("log10 [Range size"~(km^2)~"]"))+
  annotate("text", label = "R^2 == 0.25 ~ '***'",
           parse = TRUE, 
            # expression(R^2 == 0.25 ~  "***" ),
           x=10000,y=2.5,size=4,hjust = 0)+
  annotate("text", label = "R^2 == 0.21 ~ '***'",
           parse = TRUE,
           #label = expression(R^2 == 0.21 ~ "***" ),
           x=10000,y=3.5,size=4,hjust = 0)+
  annotate("segment", y = 2.5, yend = 2.5, x = 8000, xend = 5000, 
           color = "#cc79a7",linewidth = 1)+
  annotate("segment", y = 3.5, yend = 3.5, x = 8000, xend = 5000, 
           color = "#f0e442",linewidth = 1)
occplotRS

if(F){
(ndtH + 
    inset_element(occplotNB, left = 0.6, bottom = 0.58, right = 0.98, top = 0.93, align_to = 'full'))+
(adtH + 
     inset_element(occplotRS, left = 0.1, bottom = 0.58, right = 0.5, top = 0.93, align_to = 'full'))+
  plot_layout(guides = 'collect',ncol = 2,axis_titles = "collect",tag_level = "new")+
  plot_annotation(tag_levels = list(c("a)","","b)","")))
ggsave(here("products","RESU_dist_region_allPlants_v3_inset.jpeg"),
       units="cm",dpi = 300,height = 15,width = 35)
}

ndtH + adtH +
occplotNB + occplotRS +
  plot_annotation(tag_levels = "a", tag_suffix = ")")&
  theme(text = element_text(size=20))
ggsave(here("products","RESU_dist_region_allPlants_v3_fourPanel.jpeg"),
       units="cm",dpi = 300,height = 20,width = 30)
ggsave(here("products","RESU_dist_region_allPlants_v3_fourPanel.tiff"),
       units="cm",dpi = 300,height = 20,width = 30)



#### Tropical America -------------------------------------------------------------
##### NICHE ---------------------------------------------------------------------
ham <- data_hy |> dplyr::filter(band == "Tropical America") |> 
  dplyr:: filter(source == "Hyperdominants") |> 
  dplyr::select(convex) |> 
  pull()
nam <- data_hy |> dplyr::filter(band == "Tropical America") |> 
  dplyr:: filter(source == "Non-dominants") |> 
  dplyr::select(convex) |> 
  pull()
ks.test(ham,nam)
data_hy |> dplyr::filter(band == "Tropical America") |> 
  dunn_test(convex~source)

ndtam <- data_hy |> 
  dplyr::filter(band == "Tropical America") |> 
  ggplot(aes(x=convex2,y=d))+
  geom_density_ridges(aes(fill=source,scale=5,color=paste(h,nd)),linewidth=1,alpha=0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.grid = element_blank())+
  scale_color_manual(values=pal3)+
  scale_fill_manual(values=pal3)+
  scale_y_discrete(expand=c(0,0))+
  scale_x_continuous(limits=c(0,10),breaks=c(0,2,4,6,8,10))+
  labs(x="sqrt (Niche breadth)")+
  annotate("text",label=
             "ks ***
kw ***", 
           x=8, y=1.5  ,size=5,hjust=0)
ndtam

##### RANGE ---------------------------------------------------------------------
ham <- data_hy |> dplyr::filter(band == "Tropical America") |> 
  dplyr:: filter(source == "Hyperdominants") |> 
  dplyr::filter(area>1) |> 
  dplyr::select(area10) |> 
  pull()
nam <- data_hy |> dplyr::filter(band == "Tropical America") |> 
  dplyr:: filter(source == "Non-dominants") |> 
  dplyr::filter(area>1) |> 
  dplyr::select(area10) |> 
  pull()
ks.test(ham,nam)
data_hy |> dplyr::filter(band == "Tropical America") |> 
  dplyr::filter(area>1) |> 
  dunn_test(area10~source)

adtam <- data_hy |> 
  dplyr::filter(area>1) |> 
  dplyr::filter(band == "Tropical America") |> 
  ggplot(aes(x=area10,y=d))+
  geom_density_ridges(aes(fill=source,scale=5,color=paste(h,nd)),linewidth=1,alpha=0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.grid = element_blank())+
  scale_color_manual(values=pal3)+
  scale_fill_manual(values=pal3)+
  scale_y_discrete(expand=c(0,0))+
  scale_x_continuous(limits=c(2,10),breaks=c(2,4,6,8,10))+
  labs(x=expression("log10 [Range size"~(km^2)~"]"))+
  annotate("text",label=
             "ks ***
kw ***", 
           x=8, y=2  ,size=5,hjust=0)
adtam

#### Tropical Africa -----------------------------------------------------------

ham <- data_hy |> dplyr::filter(band == "Tropical Africa") |> 
  dplyr:: filter(source == "Hyperdominants") |> 
  dplyr::select(convex) |> 
  pull()
nam <- data_hy |> dplyr::filter(band == "Tropical Africa") |> 
  dplyr:: filter(source == "Non-dominants") |> 
  dplyr::select(convex) |> 
  pull()
ks.test(ham,nam)
data_hy |> dplyr::filter(band == "Tropical Africa") |> 
  dunn_test(convex~source)

ndtaf <- data_hy |> 
  dplyr::filter(band == "Tropical Africa") |> 
  ggplot(aes(x=convex2,y=d))+
  geom_density_ridges(aes(fill=source,scale=5,color=paste(h,nd)),linewidth=1,alpha=0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.grid = element_blank())+
  scale_color_manual(values=pal3)+
  scale_fill_manual(values=pal3)+
  scale_y_discrete(expand=c(0,0))+
  scale_x_continuous(limits=c(0,10),breaks=c(0,2,4,6,8,10))+
  labs(x="sqrt (Niche breadth)")+
  annotate("text",label=
             "ks ***
kw ***", 
           x=8, y=1.5  ,size=5,hjust=0)
ndtaf


ham <- data_hy |> dplyr::filter(band == "Tropical Africa") |> 
  dplyr:: filter(source == "Hyperdominants") |> 
  dplyr::filter(area>1) |> 
  dplyr::select(area10) |> 
  pull()
nam <- data_hy |> dplyr::filter(band == "Tropical Africa") |> 
  dplyr:: filter(source == "Non-dominants") |> 
  dplyr::filter(area>1) |> 
  dplyr::select(area10) |> 
  pull()
ks.test(ham,nam)
data_hy |> dplyr::filter(band == "Tropical Africa") |> 
  dplyr::filter(area>1) |> 
  dunn_test(area10~source)

adtaf <- data_hy |> 
  dplyr::filter(area>1) |> 
  dplyr::filter(band == "Tropical Africa") |> 
  ggplot(aes(x=area10,y=d))+
  geom_density_ridges(aes(fill=source,scale=5,color=paste(h,nd)),linewidth=1,alpha=0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.grid = element_blank())+
  scale_color_manual(values=pal3)+
  scale_fill_manual(values=pal3)+
  scale_y_discrete(expand=c(0,0))+
  scale_x_continuous(limits=c(2,10),breaks=c(2,4,6,8,10))+
  labs(x=expression("log10 [Range size"~(km^2)~"]"))+
  annotate("text",label=
             "ks ***
kw ***", 
           x=8, y=2  ,size=5,hjust=0)
adtaf

#### Tropical Asia -------------------------------------------------------------

ham <- data_hy |> dplyr::filter(band == "Tropical Asia") |> 
  dplyr:: filter(source == "Hyperdominants") |> 
  dplyr::select(convex) |> 
  pull()
nam <- data_hy |> dplyr::filter(band == "Tropical Asia") |> 
  dplyr:: filter(source == "Non-dominants") |> 
  dplyr::select(convex) |> 
  pull()
ks.test(ham,nam)
data_hy |> dplyr::filter(band == "Tropical Asia") |> 
  dunn_test(convex~source)

ndtas <- data_hy |> 
  dplyr::filter(band == "Tropical Asia") |> 
  ggplot(aes(x=convex2,y=d))+
  geom_density_ridges(aes(fill=source,scale=5,color=paste(h,nd)),linewidth=1,alpha=0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.grid = element_blank())+
  scale_color_manual(values=pal3)+
  scale_fill_manual(values=pal3)+
  scale_y_discrete(expand=c(0,0))+
  scale_x_continuous(limits=c(0,10),breaks=c(0,2,4,6,8,10))+
  labs(x="sqrt (Niche breadth)")+
  annotate("text",label=
             "ks ***
kw ***", 
           x=8, y=1.5  ,size=5,hjust=0)
ndtas

ham <- data_hy |> dplyr::filter(band == "Tropical Asia") |> 
  dplyr:: filter(source == "Hyperdominants") |> 
  dplyr::filter(area>1) |> 
  dplyr::select(area10) |> 
  pull()
nam <- data_hy |> dplyr::filter(band == "Tropical Asia") |> 
  dplyr:: filter(source == "Non-dominants") |> 
  dplyr::filter(area>1) |> 
  dplyr::select(area10) |> 
  pull()
ks.test(ham,nam)
data_hy |> dplyr::filter(band == "Tropical Asia") |> 
  dplyr::filter(area>1) |> 
  dunn_test(area10~source)

adtas <- data_hy |> 
  dplyr::filter(area>1) |> 
  dplyr::filter(band == "Tropical Asia") |> 
  ggplot(aes(x=area10,y=d))+
  geom_density_ridges(aes(fill=source,scale=5,color=paste(h,nd)),linewidth=1,alpha=0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.grid = element_blank())+
  scale_color_manual(values=pal3)+
  scale_fill_manual(values=pal3)+
  scale_y_discrete(expand=c(0,0))+
  scale_x_continuous(limits=c(2,10),breaks=c(2,4,6,8,10))+
  labs(x=expression("log10 [Range size"~(km^2)~"]"))+
  annotate("text",label=
             "ks ***
kw ***", 
           x=8, y=2  ,size=5,hjust=0)
adtas

### Save composite graph ---------------------------------------------------------

ndtam + ndtaf + ndtas+
  adtam + adtaf + adtas+
  plot_layout(guides = 'collect',ncol = 3,axis_titles = "collect")+
  plot_annotation(tag_levels = "a",tag_suffix = ")")
ggsave(here("products","RESU_dist_region_dominance_v2.jpeg"),
       units="cm",dpi = 300,height = 25,width = 35)

### 3.2 modelling -------------------------------------------------------------------
#### SMA -------------------------------------------------------------------

dsma <- data_hy |> select(convex2,region,area10,band,source) |>
  #unite(geo2, c("band","source")) |> 
  drop_na() |> as.data.frame()
s1 <- sma(area10 ~ convex2+source,data=dsma,type="elevation")
multcompmatrix(s1)
summary(s1)
s2 <- sma(area10 ~ convex2+source,data=dsma, multcomp=TRUE,type="shift")
multcompmatrix(s2)
s3 <- sma(area10 ~ convex2*source,data=dsma, multcomp=TRUE)
multcompmatrix(s3)


mh <- lmer(data=data_hy,area10~convex2:source:band+(1|band) )
r.squaredGLMM(mh)
tidy(mh)
mhcoef <- tidy(mh,"ran_coefs") |> print()
coef(mh)

data_hy <- data_hy |> mutate(fit.m = predict(mh, re.form=NA),
                     fit.c = predict(mh, re.form=NULL))

# make the results table

pairs(emtrends(mh,c("band","source"),var = "convex2")) |> 
  tidy() |> 
  dplyr::slice(c(3,8,12)) |>  
  dplyr::mutate(p.value = round(adj.p.value,3)) |> 
  dplyr::select(contrast,p.value) |> 
  flextable::flextable() %>%                     
  flextable::autofit() %>%    
  flextable::save_as_docx(path = here("products","Table_model3.docx"))



###Tropical Africa --------------------------------------------------
tAf <-data_hy |> 
  dplyr::filter(band=="Tropical Africa") |> 
  ggplot(aes(y=area10,x=convex2,color=source))+
  geom_point(color="gray70",alpha = 0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())+
  geom_line(aes(y=fit.c,col=source),linewidth=1)+
  scale_color_manual(values =c("#cc79a7","#f0e442"))+
  labs(y=expression("log10 [Range size"~(km^2)~"]"),x="Niche breadth",
       title="Tropical Africa")+
  scale_y_continuous(limits=c(1,8))+
  scale_x_continuous(limits=c(1,8))+
  
  annotate("segment", y = 1.5, yend = 1.5, x = 5.5, xend = 6.1, 
           color = "#f0e442",linewidth = 1)+
  annotate("text", label="slope = 0.65", x = 7.2, y =1.5,size= 4)+
  
  annotate("segment", y = 1, yend = 1, x = 5.5, xend = 6.1, 
           color ="#cc79a7",linewidth = 1)+
  annotate("text", label="slope = 0.72", x = 7.2, y =1,size=4)
tAf

a <- data_hy |>
  dplyr::filter(band=="Tropical Africa") |> 
  dplyr::filter(source == "Non-dominants") |> 
  dplyr::select(area10) |> 
  pull()
b <- data_hy |>
  dplyr::filter(band=="Tropical Africa") |> 
  dplyr::filter(source == "Hyperdominants") |> 
  dplyr::select(area10) |> 
  pull()
ks.test(a,b)

vtAf <- data_hy |>
  dplyr::filter(band=="Tropical Africa") |> 
  ggplot(aes(x=source,y=area10,fill=source))+
  geom_flat_violin(position = position_nudge(x = .12, y = 0), alpha = 0.6)+
  geom_jitter(aes(color = source),
             width=0.08, alpha = 0.6)+
 # facet_wrap(.~band)+
  scale_y_continuous(limits=c(2,8))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none",
        axis.text.y  = element_blank(),
        strip.background = element_rect(colour="white", fill="white"))+
  labs(y=expression("log10 [Range size"~(km^2)~"]"), x="",title="Tropical Africa")+
  scale_color_manual(values=pal3)+
  scale_fill_manual(values=pal3)+
  coord_flip()+
  annotate("text", label="p-value < 0.001", x = 0.7, y =6,size=5)
vtAf 

a <- data_hy |>
  dplyr::filter(band=="Tropical Africa") |> 
  dplyr::filter(source == "Non-dominants") |> 
  dplyr::select(convex2) |> 
  pull()
b <- data_hy |>
  dplyr::filter(band=="Tropical Africa") |> 
  dplyr::filter(source == "Hyperdominants") |> 
  dplyr::select(convex2) |> 
  pull()
ks.test(a,b)

vtAfN <- data_hy |>
  dplyr::filter(band=="Tropical Africa") |> 
  ggplot(aes(x=source,y=convex2,fill=source))+
  geom_flat_violin(position = position_nudge(x = .12, y = 0), alpha = 0.6)+
  geom_jitter(aes(color = source),
              width=0.08, alpha = 0.6)+
  # facet_wrap(.~band)+
  scale_y_continuous(limits=c(1,8))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_rect(colour="white", fill="white"))+
  labs(y="Niche breadth", x="",title="Tropical Africa")+
  scale_color_manual(values=pal3)+
  scale_fill_manual(values=pal3)+
  coord_flip()+
  annotate("text", label="p-value < 0.001", x = 0.7, y =6,size=5)
vtAfN 


###Tropical America --------------------------------------------------
tAm <-data_hy |> 
  dplyr::filter(band=="Tropical America") |> 
  ggplot(aes(y=area10,x=convex2,color=source))+
  geom_point(color="gray70",alpha = 0.6)+
  theme_bw()+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())+
  geom_line(aes(y=fit.c,col=source),linewidth=1)+
  scale_color_manual(values =c("#cc79a7","#f0e442"))+
  labs(y=expression("log10 [Range size"~(km^2)~"]"),x="Niche breadth",title="Tropical America")+
  scale_y_continuous(limits=c(1,8))+
  scale_x_continuous(limits=c(1,8))+
  annotate("segment", y = 1.5, yend = 1.5, x = 5.5, xend = 6.1, 
           color = "#f0e442",linewidth = 1)+
  annotate("text", label="slope = 0.66", x = 7.2, y =1.5,size= 4)+
  
  annotate("segment", y = 1, yend = 1, x = 5.5, xend = 6.1, 
           color ="#cc79a7",linewidth = 1)+
  annotate("text", label="slope = 0.72", x = 7.2, y =1,size=4)
tAm

a <- data_hy |>
  dplyr::filter(band=="Tropical America") |> 
  dplyr::filter(source == "Non-dominants") |> 
  dplyr::select(area10) |> 
  pull()
b <- data_hy |>
  dplyr::filter(band=="Tropical America") |> 
  dplyr::filter(source == "Hyperdominants") |> 
  dplyr::select(area10) |> 
  pull()
ks.test(a,b)

vtAm <- data_hy |>
  dplyr::filter(band=="Tropical America") |> 
  ggplot(aes(x=source,y=area10,fill=source))+
  geom_flat_violin(position = position_nudge(x = .12, y = 0), alpha = 0.6)+
  geom_jitter(aes(color = source),
              width=0.08, alpha = 0.6)+
  # facet_wrap(.~band)+
  scale_y_continuous(limits=c(2,8))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none",
        #axis.text.y = element_blank(),
        strip.background = element_rect(colour="white", fill="white"))+
  labs(y=expression("log10 [Range size"~(km^2)~"]"), x="",
       title="Tropical America")+
  scale_color_manual(values=pal3)+
  scale_fill_manual(values=pal3)+
  coord_flip()+
  annotate("text", label="p-value < 0.001", x = 0.7, y =6,size=5)
vtAm

a <- data_hy |>
  dplyr::filter(band=="Tropical America") |> 
  dplyr::filter(source == "Non-dominants") |> 
  dplyr::select(convex2) |> 
  pull()
b <- data_hy |>
  dplyr::filter(band=="Tropical America") |> 
  dplyr::filter(source == "Hyperdominants") |> 
  dplyr::select(convex2) |> 
  pull()
ks.test(a,b)

vtAmN <- data_hy |>
  dplyr::filter(band=="Tropical America") |> 
  ggplot(aes(x=source,y=convex2,fill=source))+
  geom_flat_violin(position = position_nudge(x = .12, y = 0), alpha = 0.6)+
  geom_jitter(aes(color = source),
              width=0.08, alpha = 0.6)+
  # facet_wrap(.~band)+
  scale_y_continuous(limits=c(1,8))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none",
        #axis.text.y = element_blank(),
        strip.background = element_rect(colour="white", fill="white"))+
  labs(y='Niche breadth', x="",title="Tropical America")+
  scale_color_manual(values=pal3)+
  scale_fill_manual(values=pal3)+
  coord_flip()+
  annotate("text", label="p-value < 0.001", x = 0.7, y =6,size=5)
vtAmN


###Tropical Asia ------------------------------------------------------------
tAs <-data_hy |> 
  dplyr::filter(band=="Tropical Asia") |> 
  ggplot(aes(y=area10,x=convex2,color=source))+
  geom_point(color="gray70",alpha = 0.6)+
  theme_bw()+
  #facet_wrap(~band)+
  theme(text=element_text(size=20,color="black"),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())+
  geom_line(aes(y=fit.c,col=source),linewidth=1)+
  scale_color_manual(values =c("#cc79a7","#f0e442"))+
  labs(y=expression("log10 [Range size"~(km^2)~"]"),x="Niche breadth",title="Tropical Asia")+
  scale_y_continuous(limits=c(1,8)) +
  scale_x_continuous(limits=c(1,8))+
  annotate("segment", y = 1.5, yend = 1.5, x = 5.5, xend = 6.1, 
           color = "#f0e442",linewidth = 1)+
  annotate("text", label="slope = 0.61", x = 7.2, y =1.5,size= 4)+
  
  annotate("segment", y = 1, yend = 1, x = 5.5, xend = 6.1, 
           color ="#cc79a7",linewidth = 1)+
  annotate("text", label="slope = 0.62", x = 7.2, y =1,size=4)
tAs

a <- data_hy |>
  dplyr::filter(band=="Tropical Asia") |> 
  dplyr::filter(source == "Non-dominants") |> 
  dplyr::select(area10) |> 
  pull()
b <- data_hy |>
  dplyr::filter(band=="Tropical Asia") |> 
  dplyr::filter(source == "Hyperdominants") |> 
  dplyr::select(area10) |> 
  pull()
ks.test(a,b)

vtAs <- data_hy |>
  dplyr::filter(band=="Tropical Asia") |> 
  ggplot(aes(x=source,y=area10,fill=source))+
  geom_flat_violin(position = position_nudge(x = .12, y = 0), alpha = 0.6)+
  geom_jitter(aes(color = source),
              width=0.08, alpha = 0.6)+
  # facet_wrap(.~band)+
  scale_y_continuous(limits=c(2,8))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_rect(colour="white", fill="white"))+
  labs(y=expression("log10 [Range size"~(km^2)~"]"), x="",
       title="Tropical Asia")+
  scale_color_manual(values=pal3)+
  scale_fill_manual(values=pal3)+
  coord_flip()+
  annotate("text", label="p-value < 0.001", x = 0.7, y =6,size=5)
vtAs

a <- data_hy |>
  dplyr::filter(band=="Tropical Asia") |> 
  dplyr::filter(source == "Non-dominants") |> 
  dplyr::select(convex2) |> 
  pull()
b <- data_hy |>
  dplyr::filter(band=="Tropical Asia") |> 
  dplyr::filter(source == "Hyperdominants") |> 
  dplyr::select(convex2) |> 
  pull()
ks.test(a,b)

vtAsN <- data_hy |>
  dplyr::filter(band=="Tropical Asia") |> 
  ggplot(aes(x=source,y=convex2,fill=source))+
  geom_flat_violin(position = position_nudge(x = .12, y = 0), alpha = 0.6)+
  geom_jitter(aes(color = source),
              width=0.08, alpha = 0.6)+
  # facet_wrap(.~band)+
  scale_y_continuous(limits=c(1,8))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.text.y = element_blank(),
        strip.background = element_rect(colour="white", fill="white"))+
  labs(y="Niche breadth", x="",title="Tropical Asia")+
  scale_color_manual(values=pal3)+
  scale_fill_manual(values=pal3)+
  coord_flip()+
  annotate("text", label="p-value < 0.001", x = 0.7, y =6,size=5)
vtAsN

#### make composite graphic ------------------------------------------------------
(tAm|(vtAm/vtAmN))/
(tAf|(vtAf/vtAfN))/
(tAs|(vtAs/vtAsN)) +
  plot_annotation(tag_levels = "a",tag_suffix = ")")
  ggsave(here("products","RESU_hyperdominants_CV_sqrt_v2.jpeg"),units="cm",dpi=300,
         width=30, height = 45)

  tAm + tAf + tAs+
    plot_annotation(tag_levels = "a",tag_suffix = ")")
  ggsave(here("products","RESU_hyperdominants_CV_sqrt_v3.jpeg"),units="cm",dpi=300,
         width=35, height = 12)

  vtAm+vtAmN+
    vtAf+vtAfN+
    vtAs+vtAsN+
  plot_layout(guides = 'collect',ncol = 3,byrow = F,
              axis_titles = 'collect')+
  plot_annotation(tag_levels = "a",tag_suffix = ")")&
  theme(legend.position="bottom")
  ggsave(here("products","RESU_hyperdominants_CV_sqrt_DIST_v2.jpeg"),
         units="cm",dpi = 300,
         width=40, height = 25)
                    
  
## 4) Occurrence analysis ------------------------------------------------------

if(F){
a <- archive(here("temp","rawdata","rm_data_20230524.tar.gz"),)

strings <- sort(c(dataGlobal$species))
sps <- c(a$path)
sps2 <- str_extract(sps, "(?<=/)[^/]+(?=\\.csv)")

b <- which(sps2 %in% strings)
a <- a[b,"path"]
a <- unlist(a, use.names = FALSE)

ex <- list()
#for (i  in 1:length(a)){ #this loop needs to be run again starting from the 94457
  for (i  in 76129:length(a)){ 
    print(i)
ex[[i]] <-  read_csv(archive_read(here("temp","rawdata","rm_data_20230524.tar.gz"),
                                  a[i]),show_col_types = F) |> 
  dplyr::select(species) |> 
  count(species)
}
ex2 <- map_df(ex,~.x)
write_csv(ex2,here("products","Occ_per_species_94457.csv"))

a <- archive(here("temp","rawdata","rm_data_20230524_missing_spp.tar.gz"))

strings <- sort(c(hdataGlobal$species))
sps <- c(a$path)
sps2 <- str_extract(sps, "(?<=/)[^/]+(?=\\.csv)")

b <- which(sps2 %in% strings)
a <- a[b,"path"]
a <- unlist(a, use.names = FALSE)

ex <- list()
for (i  in 1:length(a)){
  ex[[i]] <-  read_csv(archive_read(here("temp","rawdata","rm_data_20230524_missing_spp.tar.gz"),a[i])) |> dplyr::select(species) |> count(species)
}
ex2 <- map_df(ex,~.x)
write_csv(ex2,here("products","Occ_missing_per_species.csv"))

}

summary(lm(docc$convex2~log10(docc$n)))


occplotNB<- docc |> 
  ggplot(aes(y=convex2,x=n))+
  geom_point(color="gray70",alpha = 0.6)+
  #facet_wrap(.~higher_plant_group)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"))+
  scale_x_log10(labels=comma)+
  scale_y_continuous(limits=c(0,8))+
  # scale_y_log10()+
  geom_smooth(method="lm",col="black")+
  #scale_color_manual(values = "#56b4e9")+
  labs(x=expression(log["10"] ~ "(Occurrences)"),y="Niche breadth")+
   annotate("text",label = expression(R^2 == 0.21 ~ "," ~ p < 0.001 ),
         x=1,y=8,size=5,hjust = 0)
occplotNB

summary(lm(docc$area10~log10(docc$n)))

occplotRS<- docc|> 
  #filter( band=="Tropical Asia") |> 
  ggplot(aes(y=area10,x=n))+
  geom_point(color="gray70",alpha = 0.4)+
  #facet_wrap(.~higher_plant_group)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"))+
  scale_x_log10(labels=comma)+
 scale_y_continuous(limits=c(2,8))+
  # scale_y_log10()+
  geom_smooth(method="lm",col="black",se=T)+
  #scale_color_manual(values = "#56b4e9")+
  labs(x=expression(log["10"] ~ "(Occurrences)"),
       y=expression(log["10"] ~ "[Range size"~(km^2)~"]"))+
  annotate("text",label = expression(R^2 == 0.30 ~ "," ~ p < 0.001 ),
         x=1,y=8,size=5,hjust = 0)
occplotRS

# make composite graphic
  cp=occplotRS + occplotNB +
  plot_annotation(tag_levels = "a",tag_suffix = ")")
ggsave(here("products","RESU_OCCR_all.jpeg"),units="cm",dpi=300,
       width=30, height = 15, plot=cp)
  
#join with occ points for hyperdominant species
  
occ <- read_csv(here("products","Occ_per_species_hyper.csv"))
occ2 <- read_csv(here("products","Occ_missing_per_species_hyper.csv"))
occ3 <- read_csv(here("products","Occ_per_species.csv"),show_col_types = F)
occ <- rbind(occ,occ2)
data_hy <- left_join(data_hy,occ, by="species")
data_hy <- data_hy |> drop_na()

# Tropica Africa

daf <- data_hy |> 
  filter(band=="Tropical Africa" & source=="Hyperdominants")

summary(lm(daf$convex2~log10(daf$n)))


occplotNBaf<- daf |> 
  ggplot(aes(y=convex2,x=n))+
  geom_point(color="gray70",alpha = 0.6)+
  #facet_wrap(.~higher_plant_group)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"))+
  scale_x_log10(labels=comma)+
  scale_y_continuous(limits=c(0,8))+
  # scale_y_log10()+
  geom_smooth(method="lm",col="black")+
  #scale_color_manual(values = "#56b4e9")+
  labs(x=expression(log["10"] ~ "(Occurrences)"),y="Niche breadth",title="Tropical Africa")+
  annotate("text",label = expression(R^2 == 0.04 ~ "," ~ p == 0.01 ),
           x=1,y=8,size=5,hjust = 0)
occplotNBaf

summary(lm(daf$area10~log10(daf$n)))

occplotRSaf<- daf|> 
  #filter( band=="Tropical Asia") |> 
  ggplot(aes(y=area10,x=n))+
  geom_point(color="gray70",alpha = 0.6)+
  #facet_wrap(.~higher_plant_group)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"))+
  scale_x_log10(labels=comma)+
  scale_y_continuous(limits=c(2,8))+
  # scale_y_log10()+
  geom_smooth(method="lm",col="black")+
  #scale_color_manual(values = "#56b4e9")+
  labs(x=expression(log["10"] ~ "(Occurrences)"),y=expression(log["10"] ~ "[Range size"~(km^2)~"]"),
       ,title="Tropical Africa")+
  annotate("text",label = expression(R^2 == 0.13 ~ "," ~ p < 0.001 ),
           x=1,y=8,size=5,hjust = 0)
occplotRSaf

# Tropica America

dam <- data_hy |> 
  filter(band=="Tropical America" & source=="Hyperdominants")

summary(lm(dam$convex2~log10(dam$n)))


occplotNBam<- dam |> 
  ggplot(aes(y=convex2,x=n))+
  geom_point(color="gray70",alpha = 0.6)+
  #facet_wrap(.~higher_plant_group)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"))+
  scale_x_log10(labels=comma)+
  scale_y_continuous(limits=c(0,8))+
  # scale_y_log10()+
  geom_smooth(method="lm",col="black")+
  #scale_color_manual(values = "#56b4e9")+
  labs(x=expression(log["10"] ~ "(Occurrences)"),y="Niche breadth",title="Tropical America")+
  annotate("text",label = expression(R^2 == 0.36 ~ "," ~ p < 0.001 ),
           x=1,y=8,size=5,hjust = 0)
occplotNBam

summary(lm(dam$area10~log10(dam$n)))

occplotRSam<- dam|> 
  #filter( band=="Tropical Asia") |> 
  ggplot(aes(y=area10,x=n))+
  geom_point(color="gray70",alpha = 0.6)+
  #facet_wrap(.~higher_plant_group)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"))+
  scale_x_log10(labels=comma)+
  scale_y_continuous(limits=c(2,8))+
  # scale_y_log10()+
  geom_smooth(method="lm",col="black")+
  #scale_color_manual(values = "#56b4e9")+
  labs(x=expression(log["10"] ~ "(Occurrences)"),y=expression(log["10"]~" [Range size"~(km^2)~"]"),
       ,title="Tropical America")+
  annotate("text",label = expression(R^2 == 0.37 ~ "," ~ p < 0.001 ),
           x=1,y=8,size=5,hjust = 0)
occplotRSam

# Tropica Asia

das <- data_hy |> 
  filter(band=="Tropical Asia" & source=="Hyperdominants")

summary(lm(das$convex2~log10(das$n)))

occplotNBas<- das |> 
  ggplot(aes(y=convex2,x=n))+
  geom_point(color="gray70",alpha = 0.6)+
  #facet_wrap(.~higher_plant_group)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"))+
  scale_x_log10(labels=comma)+
  scale_y_continuous(limits=c(0,8))+
  # scale_y_log10()+
  geom_smooth(method="lm",col="black")+
  #scale_color_manual(values = "#56b4e9")+
  labs(x=expression(log["10"] ~ "(Occurrences)"),y="Niche breadth",title="Tropical Asia")+
  annotate("text",label = expression(R^2 == 0.23 ~ "," ~ p < 0.001 ),
           x=1,y=8,size=5,hjust = 0)
occplotNBas

summary(lm(das$area10~log10(das$n)))

occplotRSas<- das|> 
  #filter( band=="Tropical Asia") |> 
  ggplot(aes(y=area10,x=n))+
  geom_point(color="gray70",alpha = 0.6)+
  #facet_wrap(.~higher_plant_group)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour="white", fill="white"))+
  scale_x_log10(labels=comma)+
  scale_y_continuous(limits=c(2,8))+
  # scale_y_log10()+
  geom_smooth(method="lm",col="black")+
  #scale_color_manual(values = "#56b4e9")+
  labs(x=expression(log["10"] ~ "(Occurrences)"),
       y=expression(log["10"] ~"[Range size"~(km^2)~"]"),
       ,title="Tropical Asia")+
  annotate("text",label = expression(R^2 == 0.02 ~ "," ~ p == 0.06 ),
           x=1,y=8,size=5,hjust = 0)
occplotRSas

#make composite graph
occplotNBaf+ occplotNBam+ occplotNBas +
  occplotRSaf+ occplotRSam + occplotRSas+
  plot_annotation(tag_levels = "a",tag_suffix = ")")
ggsave(here("products","RESU_hyperdominants_OCCR_BAND_v2.jpeg"),units="cm",dpi=300,
       width=40, height = 25)



#some stats
data_hy |> dplyr::filter(source=="Hyperdominants" & band=="Tropical Africa") |> View()


## 5) NULL models -------------------------------------------
nf4p <- read_csv("dataNUll.csv")

nf4p <- left_join(nf4p,d)

df <- nf4p |> dplyr::select(cvn,cv,band, higher_plant_group,type) |> 
  dplyr::filter(cv!=0) |> 
  group_by(cvn,band, higher_plant_group,type) |> summarise(mean=mean(cv), 
                                                              median=median(cv))
ns <- ggplot(df,aes(x=mean, group=type))+
  #geom_density(fill="skyblue3")+
  geom_histogram(fill="gray68",alpha=0.6,color="black")+
  geom_vline(data=subset(df, df$cvn=="convex"),aes(xintercept=mean),
             color="red3",linetype="dashed")+
  facet_grid(band~higher_plant_group)+
  labs(x="Mean niche size", y="Frequency")+
  #ggtitle("Coefficients log10(hypervolume) ~ log10(Geographic area)")+
  theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=20),
          legend.title = element_blank(),
          legend.position = "bottom",
          strip.background = element_rect(colour="white", fill="white"))+
  scale_y_continuous(expand=c(0,0))
ns
                             
#### get coefs ------------------------------------------------------------

coefs = as.data.frame(matrix(NA, ncol=4, nrow = length(grep("convex",colnames(nf4)))))
colnames(coefs) = c("cv","a","b","R2")

for(j in 1:length(grep("convex",colnames(nf4)))){
  coefs[j,1] <- names(nf4)[grep("convex",colnames(nf4))[j]]
  td <- nf4 |> dplyr::select(area,coefs[j,1])
  td <- td[which(td[,2]!=0),]
  td <- as.data.frame(td)
  #model <- lm(td$convex~log10(td$area))
  model <- lm(log10(td$area)~sqrt(td$convex))
  s <- summary(model)
  coefs[j,2] <- model$coefficients[1]
  coefs[j,3] <- model$coefficients[2]
  coefs[j,4] <- round(s$r.squared,2)
}


#get the quantiles
a <- coefs|> dplyr::select(b)
ci <-  a |> summarize(lower = quantile(b, probs = .05),
                      upper = quantile(b, probs = .95))

##### Make the graphic -----------------------------------------------
cplot <- ggplot(coefs,aes(b))+
  #geom_density(fill="skyblue3")+
  geom_histogram(bins=30,fill="grey80",color="gray23")+
 geom_vline(xintercept = 0.59,color="red",linetype="solid",linewidth=1)+
 # geom_vline(xintercept = 1.5,color="red",linetype="solid",linewidth=1)+
  geom_vline(data=ci, aes(xintercept = lower),color="grey20",linetype="dashed")+
  geom_vline(data=ci, aes(xintercept = upper),color="grey20",linetype="dashed")+
 # labs(x="Slope values", y="Frequency", title="antes")+
  labs(x="Slope values", y="Frequency")+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,12))+
  #ggtitle("Coefficients log10(hypervolume) ~ log10(Geographic area)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(colour="white", fill="white"))+
  annotate("text",label="z = 32, p < 0.0001",x=0.8,y=12,size=5)
cplot


ggsave(here("products","RESU_general_nullmodels_sqrt.jpeg"),units = "cm",
       height = 15,width = 20)

co <- coefs[-1,]
mean(co$b)

# Known mean and standard deviation of the population
mu <- mean(co$b)  # Population mean
sigma <- sd(co$b)  # Population standard deviation
n <- length(co$b)

# Value to test
x <- coefs[1,"b"]

z <- (x - mu) / (sigma / sqrt(n));print(z)
p_value <- 2 * pnorm(-abs(z)); print(p_value)


nf5 <- nf4 |> unite(class, c(band,higher_plant_group))
length(unique(nf5$class))

clist <- list()
for (i in unique(nf5$class)){
  temp <- subset(nf5, nf5$class==i)
  
  temp2 <-  as.data.frame(matrix(NA, ncol=7, nrow = length(grep("convex",colnames(temp)))))
  colnames(temp2) = c("cv","a","b","R2","class","ciL","ciH")
  
  for(j in 1:length(grep("convex",colnames(temp)))){
    temp2[j,1] <- names(temp)[grep("convex",colnames(temp))[j]]
    td <- temp |> dplyr::select(area,temp2[j,1])
    td <- td[which(td[,2]!=0),]
    td <- as.data.frame(td)
    model <- lm(td[,2]~log10(td[,1]))
    s <- summary(model)
    temp2[j,2] <- model$coefficients[1]
    temp2[j,3] <- model$coefficients[2]
    temp2[j,4] <- round(s$r.squared,2)
    temp2[j,5] <- i
  }
  
  a <- temp2|> dplyr::select(b)
  
  temp2[,c("ciL","ciH")]<-  a |> summarize(lower = quantile(b, probs = .05),
                        upper = quantile(b, probs = .95))
  clist[[i]] <- temp2
  
}
coefs_c <- map_df(clist,~.x)
coefs_c <- coefs_c |> separate_wider_delim(class, delim="_",names=c("band","HPG"))

cplot2 <- ggplot(coefs_c,aes(b))+
  #geom_density(fill="skyblue3")+
  geom_histogram(bins=20,fill="grey",color="gray23")+
  facet_grid(band~HPG)+
  geom_vline(data=subset(coefs_c, coefs_c$cv=="convex"),aes(xintercept = b),color="red",linetype="solid")+
  geom_vline(data=coefs_c, aes(xintercept = ciL),color="grey20",linetype="dashed")+
  geom_vline(data=coefs_c, aes(xintercept = ciH),color="grey20",linetype="dashed")+
  
  labs(x="Slope of null model (beta coefficient)", y="Frequency")+
  #ggtitle("Coefficients log10(hypervolume) ~ log10(Geographic area)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(colour="white", fill="white"))

cplot2
ggsave(here("products","RESU_species_nullmodels.jpeg"),units = "cm",
       height = 25,width = 30)


