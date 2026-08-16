#!/usr/bin/env Rscript

# ============================================================
# 06_plot_coverage.R
#
# Generate coverage plots from:
#
#   results/summary/all_depth.tsv
#   metadata/experimental_groups.tsv
#
# Outputs:
#
#   results/figures/coverage_by_group.png
#   results/figures/median_depth_heatmap.png
#
# ============================================================


library(ggplot2)
library(dplyr)
library(ggpubr)



# ============================================================
# Paths
# ============================================================

args <- commandArgs(trailingOnly = FALSE)

file_arg <- grep("^--file=", args, value = TRUE)

if(length(file_arg) == 1){
  
  script_path <- normalizePath(
    sub("^--file=", "", file_arg)
  )
  
  repo_root <- dirname(dirname(script_path))
  
}else{
  
  repo_root <- normalizePath(".")
  
}


coverage_file <- file.path(
  repo_root,
  "results",
  "summary",
  "all_depth.tsv"
)

metadata_file <- file.path(
  repo_root,
  "metadata",
  "experimental_groups.tsv"
)

output_dir <- file.path(
  repo_root,
  "results",
  "figures"
)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)



# ============================================================
# Check files
# ============================================================

if(!file.exists(coverage_file)){
  
  stop(
    paste(
      "[ERROR] Coverage file not found:",
      coverage_file
    )
  )
  
}


if(!file.exists(metadata_file)){
  
  stop(
    paste(
      "[ERROR] Metadata file not found:",
      metadata_file
    )
  )
  
}



# ============================================================
# Read data
# ============================================================

message("[INFO] Reading coverage data")


coverage <- read.delim(
  coverage_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  na.strings = ""
)


metadata <- read.delim(
  metadata_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  na.strings = ""
)


coverage$position <- as.numeric(coverage$position)

coverage$coverage <- as.numeric(coverage$coverage)



# ============================================================
# Fix metadata IDs
#
# R tends to convert:
#
#   01 -> 1
#   02 -> 2
#
# We restore two-digit IDs before joining.
# ============================================================

metadata$ID <- sprintf(
  "%02d",
  as.integer(metadata$ID)
)



# ============================================================
# Extract animal ID from sample names
# ============================================================

coverage$ID <- NA_character_



# BALF samples:
#
# 01_BALF -> 01
# 12_BALF -> 12

balf <- grepl(
  "^[0-9]{2}_BALF$",
  coverage$sample
)


coverage$ID[balf] <- substr(
  coverage$sample[balf],
  1,
  2
)



# Nasal swabs:
#
# animal5_3_NS  -> 05
# animal16_9_NS -> 16

ns <- grepl(
  "^animal[0-9]+_",
  coverage$sample
)


animal_number <- sub(
  "^animal([0-9]+)_.*$",
  "\\1",
  coverage$sample[ns]
)


coverage$ID[ns] <- sprintf(
  "%02d",
  as.integer(animal_number)
)



# ============================================================
# Add experimental groups
# ============================================================

coverage <- left_join(
  coverage,
  metadata[, c("ID", "Group")],
  by = "ID"
)



# Inoculum has no animal ID

coverage$Group[
  coverage$sample == "H1N2"
] <- "Inoculum"



coverage$Group <- factor(
  coverage$Group,
  levels = c(
    "Inoculum",
    "Vaccinated",
    "Nonvaccinated"
  ),
  ordered = TRUE
)



# ============================================================
# Segment order
# ============================================================

coverage$segment <- factor(
  coverage$segment,
  levels = c(
    "PB2",
    "PB1",
    "PA",
    "HA",
    "NP",
    "NA",
    "M",
    "NS"
  ),
  ordered = TRUE
)



# ============================================================
# Pretty sample names
# ============================================================

coverage$sample_label <- coverage$sample



coverage$sample_label[
  coverage$sample == "H1N2"
] <- "Inoculum"



# BALF:
#
# 01_BALF -> 1 BALF

idx <- grepl(
  "^[0-9]{2}_BALF$",
  coverage$sample
)


coverage$sample_label[idx] <- paste0(
  as.integer(
    substr(
      coverage$sample[idx],
      1,
      2
    )
  ),
  " BALF"
)



# Nasal swabs:
#
# animal5_3_NS -> 5 3 dpi

idx <- grepl(
  "^animal[0-9]+_[0-9]+_NS$",
  coverage$sample
)


tmp <- strsplit(
  coverage$sample[idx],
  "_"
)


coverage$sample_label[idx] <- vapply(
  tmp,
  function(x){
    
    animal <- sub(
      "^animal",
      "",
      x[1]
    )
    
    day <- x[2]
    
    paste(
      animal,
      paste0(
        day,
        " dpi"
      )
    )
    
  },
  character(1)
)



# ============================================================
# Remove zero coverage for log-scale plots
# ============================================================

coverage_log <- coverage[
  coverage$coverage > 0,
]



# ============================================================
# Colours
# ============================================================

textcol <- "gray40"



# Inoculum

colour_inoculum <- c(
  "Inoculum" = "#D35400"
)



# Vaccinated animals

vaccinated_samples <- sort(
  unique(
    coverage_log$sample_label[
      coverage_log$Group == "Vaccinated"
    ]
  )
)


blue_palette <- colorRampPalette(
  c(
    "#AED6F1",
    "#1B4F72"
  )
)


vaccinated_colours <- blue_palette(
  length(vaccinated_samples)
)


names(vaccinated_colours) <- vaccinated_samples



# Non-vaccinated animals

nonvaccinated_samples <- sort(
  unique(
    coverage_log$sample_label[
      coverage_log$Group == "Nonvaccinated"
    ]
  )
)


green_palette <- colorRampPalette(
  c(
    "#A2D9CE",
    "#0B5345"
  )
)


nonvaccinated_colours <- green_palette(
  length(nonvaccinated_samples)
)


names(nonvaccinated_colours) <- nonvaccinated_samples



# ============================================================
# Inoculum coverage
# ============================================================

message("[INFO] Plotting inoculum")


coverage_inoculum <- coverage_log[
  coverage_log$Group == "Inoculum",
]


COVI <- ggplot(
  coverage_inoculum,
  aes(
    x = position,
    y = coverage,
    colour = sample_label,
    group = sample_label
  )
) +
  geom_point(size=0.1) +
  geom_line(linewidth=0.25) +
  scale_y_log10(
    breaks=c(1,10,100,1000,10000),
    labels=c("1","10","100","1k","10k")
  ) +
  coord_cartesian(
    ylim=c(1,10000)
  ) +
  scale_colour_manual(
    values=colour_inoculum
  ) +
  geom_hline(
    yintercept=100,
    color="gray30",
    linetype="dashed",
    linewidth=0.4
  ) +
  labs(
    x="Genome position (nt)",
    y="Read depth (x)",
    title="A. Inoculum",
    colour="Sample"
  ) +
  facet_grid(
    ~segment,
    scales="free_x",
    space="free",
    switch="x"
  ) +
  theme_classic() +
  theme(
    panel.spacing=unit(0,"lines"),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_text(
      size=10,
      colour=textcol
    ),
    axis.title.x=element_text(
      size=10,
      face="bold",
      colour=textcol
    ),
    axis.title.y=element_text(
      size=10,
      face="bold",
      colour=textcol
    ),
    axis.ticks.x=element_blank(),
    legend.title=element_text(
      size=10,
      face="bold",
      colour=textcol
    ),
    legend.position="right",
    legend.justification="left",
    legend.direction="vertical",
    legend.margin=margin(
      grid::unit(0,"cm")
    ),
    legend.text=element_text(
      colour=textcol,
      size=10
    ),
    legend.key.height=grid::unit(
      0.6,
      "cm"
    ),
    legend.key.width=grid::unit(
      0.3,
      "cm"
    ),
    plot.title=element_text(
      hjust=0,
      size=12,
      face="bold",
      colour=textcol
    )
  ) +
  guides(
    colour=guide_legend(
      override.aes=list(
        size=2.5
      )
    )
  ) +
  theme(
    strip.text.x=element_text(
      size=10,
      face="bold",
      colour=textcol
    ),
    strip.placement="outside",
    strip.background.x=element_rect(
      color=NA,
      fill=NA
    ),
    strip.background.y=element_rect(
      color=NA,
      fill=NA
    )
  )



# ============================================================
# Vaccinated animals coverage
# ============================================================

message("[INFO] Plotting vaccinated animals")


coverage_vaccinated <- coverage_log[
  coverage_log$Group == "Vaccinated",
]


COVV <- ggplot(
  coverage_vaccinated,
  aes(
    x = position,
    y = coverage,
    colour = sample_label,
    group = sample_label
  )
) +
  geom_point(size=0.1) +
  geom_line(linewidth=0.25) +
  scale_y_log10(
    breaks=c(1,10,100,1000,10000),
    labels=c("1","10","100","1k","10k")
  ) +
  coord_cartesian(
    ylim=c(1,10000)
  ) +
  scale_colour_manual(
    values=vaccinated_colours
  ) +
  geom_hline(
    yintercept=100,
    color="gray30",
    linetype="dashed",
    linewidth=0.4
  ) +
  labs(
    x="Genome position (nt)",
    y="Read depth (x)",
    title="B. Vaccinated animals",
    colour="Sample"
  ) +
  facet_grid(
    ~segment,
    scales="free_x",
    space="free",
    switch="x"
  ) +
  theme_classic() +
  theme(
    panel.spacing=unit(0,"lines"),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_text(
      size=10,
      colour=textcol
    ),
    axis.title.x=element_text(
      size=10,
      face="bold",
      colour=textcol
    ),
    axis.title.y=element_text(
      size=10,
      face="bold",
      colour=textcol
    ),
    axis.ticks.x=element_blank(),
    legend.title=element_text(
      size=10,
      face="bold",
      colour=textcol
    ),
    legend.position="right",
    legend.justification="left",
    legend.direction="vertical",
    legend.margin=margin(
      grid::unit(0,"cm")
    ),
    legend.text=element_text(
      colour=textcol,
      size=10
    ),
    legend.key.height=grid::unit(
      0.6,
      "cm"
    ),
    legend.key.width=grid::unit(
      0.3,
      "cm"
    ),
    plot.title=element_text(
      hjust=0,
      size=12,
      face="bold",
      colour=textcol
    )
  ) +
  guides(
    colour=guide_legend(
      override.aes=list(
        size=2.5
      )
    )
  ) +
  theme(
    strip.text.x=element_text(
      size=10,
      face="bold",
      colour=textcol
    ),
    strip.placement="outside",
    strip.background.x=element_rect(
      color=NA,
      fill=NA
    ),
    strip.background.y=element_rect(
      color=NA,
      fill=NA
    )
  )



# ============================================================
# Non-vaccinated animals coverage
# ============================================================

message("[INFO] Plotting non-vaccinated animals")


coverage_nonvaccinated <- coverage_log[
  coverage_log$Group == "Nonvaccinated",
]


COVN <- ggplot(
  coverage_nonvaccinated,
  aes(
    x = position,
    y = coverage,
    colour = sample_label,
    group = sample_label
  )
) +
  geom_point(size=0.1) +
  geom_line(linewidth=0.25) +
  scale_y_log10(
    breaks=c(1,10,100,1000,10000),
    labels=c("1","10","100","1k","10k")
  ) +
  coord_cartesian(
    ylim=c(1,10000)
  ) +
  scale_colour_manual(
    values=nonvaccinated_colours
  ) +
  geom_hline(
    yintercept=100,
    color="gray30",
    linetype="dashed",
    linewidth=0.4
  ) +
  labs(
    x="Genome position (nt)",
    y="Read depth (x)",
    title="C. Non-vaccinated animals",
    colour="Sample"
  ) +
  facet_grid(
    ~segment,
    scales="free_x",
    space="free",
    switch="x"
  ) +
  theme_classic() +
  theme(
    panel.spacing=unit(0,"lines"),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_text(
      size=10,
      colour=textcol
    ),
    axis.title.x=element_text(
      size=10,
      face="bold",
      colour=textcol
    ),
    axis.title.y=element_text(
      size=10,
      face="bold",
      colour=textcol
    ),
    axis.ticks.x=element_blank(),
    legend.title=element_text(
      size=10,
      face="bold",
      colour=textcol
    ),
    legend.position="right",
    legend.justification="left",
    legend.direction="vertical",
    legend.margin=margin(
      grid::unit(0,"cm")
    ),
    legend.text=element_text(
      colour=textcol,
      size=10
    ),
    legend.key.height=grid::unit(
      0.6,
      "cm"
    ),
    legend.key.width=grid::unit(
      0.3,
      "cm"
    ),
    plot.title=element_text(
      hjust=0,
      size=12,
      face="bold",
      colour=textcol
    )
  ) +
  guides(
    colour=guide_legend(
      override.aes=list(
        size=2.5
      )
    )
  ) +
  theme(
    strip.text.x=element_text(
      size=10,
      face="bold",
      colour=textcol
    ),
    strip.placement="outside",
    strip.background.x=element_rect(
      color=NA,
      fill=NA
    ),
    strip.background.y=element_rect(
      color=NA,
      fill=NA
    )
  )



# ============================================================
# Combine coverage plots
# ============================================================

message("[INFO] Combining coverage plots")


COV_ALL <- ggarrange(
  COVI,
  COVV,
  COVN,
  ncol=1,
  nrow=3,
  align="v",
  heights=c(
    0.75,
    1,
    1
  )
)


ggsave(
  file.path(
    output_dir,
    "coverage_by_group.png"
  ),
  COV_ALL,
  width=12,
  height=7,
  dpi=300
)



# ============================================================
# Median depth heatmap
# ============================================================

message("[INFO] Building median depth heatmap")


median_depth <- coverage %>%
  group_by(
    sample_label,
    segment,
    Group
  ) %>%
  summarise(
    median_depth=median(
      coverage,
      na.rm=TRUE
    ),
    .groups="drop"
  )



median_depth$Group <- factor(
  median_depth$Group,
  levels=c(
    "Inoculum",
    "Vaccinated",
    "Nonvaccinated"
  ),
  ordered=TRUE
)



HEAT <- ggplot(
  median_depth,
  aes(
    x=segment,
    y=sample_label,
    fill=median_depth
  )
) +
  geom_tile(
    colour="white",
    linewidth=0.2
  ) +
  scale_fill_viridis_c(
    trans="log10",
    name="Median depth"
  ) +
  facet_grid(
    Group~.,
    scales="free_y",
    space="free_y",
    switch="y"
  ) +
  labs(
    x="Genome segment",
    y=NULL,
    title="Median read depth per genome segment"
  ) +
  theme_classic() +
  theme(
    axis.text.x=element_text(
      size=10,
      face="bold",
      colour=textcol
    ),
    axis.text.y=element_text(
      size=10,
      colour=textcol
    ),
    axis.title.x=element_text(
      size=10,
      face="bold",
      colour=textcol
    ),
    axis.ticks=element_blank(),
    plot.title=element_text(
      hjust=0,
      size=12,
      face="bold",
      colour=textcol
    ),
    strip.placement="outside",
    strip.background=element_blank(),
    strip.text.y=element_text(
      size=10,
      face="bold",
      colour=textcol
    ),
    legend.title=element_text(
      size=10,
      face="bold",
      colour=textcol
    ),
    legend.text=element_text(
      size=10,
      colour=textcol
    )
  )



ggsave(
  file.path(
    output_dir,
    "median_depth_heatmap.png"
  ),
  HEAT,
  width=7,
  height=7,
  dpi=300
)



# ============================================================
# Done
# ============================================================

message("")
message("============================================================")
message("[SUCCESS] Coverage figures generated")
message("============================================================")
message("")
message("[INFO] Coverage plot:")
message(
  file.path(
    output_dir,
    "coverage_by_group.png"
  )
)
message("")
message("[INFO] Median-depth heatmap:")
message(
  file.path(
    output_dir,
    "median_depth_heatmap.png"
  )
)