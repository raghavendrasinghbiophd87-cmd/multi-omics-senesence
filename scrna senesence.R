library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dplyr)
library(patchwork)

graphics.off()

# 1. Build Multi-Pathway Phospho-Site Activation Matrix (Proliferating -> Early Senescent -> SASP Senescent)
phospho_mat = matrix(
  c(
    # DNA Damage & p53 Axis (Activated in Early & SASP Senescence)
    -1.15,  0.85,  0.95,   # ATM (Ser1981)
    -1.20,  0.90,  1.05,   # H2AX (Ser139 / yH2AX)
    -1.10,  0.95,  0.90,   # CHK2 (Thr68)
    -1.25,  1.10,  0.95,   # TP53 (Ser15 / Ser20)
    # Cell Cycle Pocket Proteins (Hypo-phosphorylated / Deactivated in Senescence)
    1.35, -0.85, -1.10,   # RB1 (Ser780 / Ser807/811) [Loss of phospho = Active Arrest]
    1.30, -0.90, -1.05,   # CDK1/2 (Thr160/161)
    # mTORC1 & Translational SASP Axis
    -0.95,  0.25,  1.25,   # mTOR (Ser2448)
    -0.90,  0.20,  1.20,   # 4E-BP1 (Thr37/46)
    -0.85,  0.15,  1.15,   # RPS6KB1 (Thr389)
    # cGAS-STING & Interferon Axis
    -0.80, -0.10,  1.30,   # STING1 (Ser366)
    -0.85, -0.05,  1.32,   # TBK1 (Ser172)
    -0.80, -0.10,  1.28,   # IRF3 (Ser386 / Ser396)
    # p38 MAPK / NF-kB / IL-6 STAT3 SASP Axis
    -0.95,  0.30,  1.25,   # p38 MAPK (Thr180/Tyr182)
    -0.90,  0.15,  1.30,   # IKKB (Ser177/181)
    -0.92,  0.20,  1.32,   # RELA/p65 (Ser536)
    -0.95,  0.10,  1.35,   # JAK2 (Tyr1007/1008)
    -0.98,  0.05,  1.38,   # STAT3 (Tyr705)
    -0.95,  0.25,  1.35    # STAT3 (Ser727)
  ),
  nrow = 18, byrow = TRUE,
  dimnames = list(
    c("ATM (Ser1981)", "H2AX (Ser139)", "CHK2 (Thr68)", "p53 (Ser15/Ser20)",
      "RB1 (Ser780/807) [Hypo-P]", "CDK1/2 (Thr160/161)",
      "mTOR (Ser2448)", "4E-BP1 (Thr37/46)", "S6K1 (Thr389)",
      "STING1 (Ser366)", "TBK1 (Ser172)", "IRF3 (Ser386/396)",
      "p38 MAPK (Thr180/Tyr182)", "IKKb (Ser177/181)", "RelA/p65 (Ser536)",
      "JAK2 (Tyr1007/1008)", "STAT3 (Tyr705)", "STAT3 (Ser727)"),
    c("Proliferating", "Early_Senescent", "SASP_Senescent")
  )
)

row_annot = data.frame(
  Pathway = factor(c(rep("p53 / DDR", 4), rep("Cell Cycle", 2), rep("mTORC1", 3),
                     rep("cGAS-STING", 3), rep("p38 / NF-kB", 3), rep("IL-6 / STAT3", 3)),
                   levels = c("p53 / DDR", "Cell Cycle", "mTORC1", "cGAS-STING", "p38 / NF-kB", "IL-6 / STAT3")),
  row.names = rownames(phospho_mat)
)

col_annot = data.frame(
  CellState = factor(c("Proliferating", "Early_Senescent", "SASP_Senescent"),
                     levels = c("Proliferating", "Early_Senescent", "SASP_Senescent")),
  row.names = colnames(phospho_mat)
)

ann_colors = list(
  CellState = c(Proliferating = "#2563EB", Early_Senescent = "#F59E0B", SASP_Senescent = "#DC2626"),
  Pathway   = c("p53 / DDR" = "#10B981", "Cell Cycle" = "#6366F1", "mTORC1" = "#D97706",
                "cGAS-STING" = "#0D9488", "p38 / NF-kB" = "#E11D48", "IL-6 / STAT3" = "#7C3AED")
)

out_dir = file.path(Sys.getenv("USERPROFILE"), "Desktop", "Senescence_Results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

hm_path = file.path(out_dir, "16_MultiPathway_Phosphorylation_Heatmap.png")

# Save high-res multi-pathway phospho heatmap
pheatmap(phospho_mat,
         annotation_col    = col_annot,
         annotation_row    = row_annot,
         annotation_colors = ann_colors,
         cluster_cols      = FALSE,
         cluster_rows      = FALSE,
         gaps_row          = c(4, 6, 9, 12, 15),
         display_numbers   = TRUE,
         number_format     = "%.2f",
         fontsize_number   = 8.5,
         fontsize_row      = 10,
         fontsize_col      = 11,
         angle_col         = 0,
         color             = colorRampPalette(c("#1D4ED8", "white", "#DC2626"))(100),
         main              = "Multi-Pathway Phosphorylation Dynamics in Cancer Senescence",
         filename          = hm_path,
         width             = 9.5,
         height            = 8.5)

# Also draw to RStudio Plots pane
grid::grid.newpage()
pheatmap(phospho_mat,
         annotation_col    = col_annot,
         annotation_row    = row_annot,
         annotation_colors = ann_colors,
         cluster_cols      = FALSE,
         cluster_rows      = FALSE,
         gaps_row          = c(4, 6, 9, 12, 15),
         display_numbers   = TRUE,
         number_format     = "%.2f",
         fontsize_number   = 8.5,
         fontsize_row      = 10,
         fontsize_col      = 11,
         angle_col         = 0,
         color             = colorRampPalette(c("#1D4ED8", "white", "#DC2626"))(100),
         main              = "Multi-Pathway Phosphorylation Dynamics in Cancer Senescence")

browseURL(hm_path)
library(ggplot2)
library(ggrepel)
library(dplyr)

graphics.off()

# 1. Protein Backbone Lengths
backbone_df = data.frame(
  Protein = factor(c("gp130 (IL6ST)", "TYK2 / JAK1/2", "SHP2 (PTPN11)", "SOCS3", "STAT1 / STAT3"),
                   levels = c("STAT1 / STAT3", "SOCS3", "SHP2 (PTPN11)", "TYK2 / JAK1/2", "gp130 (IL6ST)")),
  y_pos   = c(5, 4, 3, 2, 1),
  aa_len  = c(918, 1187, 593, 225, 770)
)

# 2. Extended Phosphorylation Sites & Their Impact on NF-kB & Senescence
ext_phospho = data.frame(
  Protein  = factor(c(
    "gp130 (IL6ST)", "gp130 (IL6ST)", "gp130 (IL6ST)",
    "TYK2 / JAK1/2", "TYK2 / JAK1/2", "TYK2 / JAK1/2",
    "SHP2 (PTPN11)", "SHP2 (PTPN11)",
    "SOCS3", "SOCS3",
    "STAT1 / STAT3", "STAT1 / STAT3", "STAT1 / STAT3"
  ), levels = c("STAT1 / STAT3", "SOCS3", "SHP2 (PTPN11)", "TYK2 / JAK1/2", "gp130 (IL6ST)")),
  y_base   = c(5, 5, 5, 4, 4, 4, 3, 3, 2, 2, 1, 1, 1),
  Position = c(759, 782, 905, 570, 1007, 1054, 542, 580, 204, 221, 701, 705, 727),
  Label    = c(
    "Y759: Recruits SHP2 -> IKK/NF-kB",
    "S782 (MK2): Endocytosis Rheostat",
    "Y905/Y915: STAT1/3 Docking",
    "Y570 (JAK2): Inhibitory Brake",
    "Y1007/1034: Kinase Activation",
    "Y1054 (TYK2): IFN/SASP Synergy",
    "Y542: Activates PI3K/AKT -> IKKb",
    "Y580: Boosts ERK -> p65 S536",
    "Y204: Degrades SOCS3 (Locks SASP ON)",
    "Y221: Blocks Elongin B/C Binding",
    "Y701 (STAT1): CXCL10/Immune SASP",
    "Y705 (STAT3): Nuclear Dimerization",
    "S727 (STAT3): p300/Ac-p65 & mtROS"
  ),
  Effect   = c(
    "Activates NF-kB / SASP", "Inhibitory / Rheostat", "Activates NF-kB / SASP",
    "Inhibitory / Rheostat", "Activates NF-kB / SASP", "Activates NF-kB / SASP",
    "Activates NF-kB / SASP", "Activates NF-kB / SASP",
    "Disables Brake (Locks SASP)", "Disables Brake (Locks SASP)",
    "Immune SASP Shift", "Activates NF-kB / SASP", "Activates NF-kB / SASP"
  )
)

# 3. Build Zero-Overlap Phospho-Site & NF-kB Impact Chart
p_ext_il6 = ggplot() +
  geom_segment(data = backbone_df, aes(x = 0, xend = aa_len, y = y_pos, yend = y_pos),
               color = "#CBD5E1", linewidth = 3.5, lineend = "round") +
  geom_segment(data = ext_phospho, aes(x = Position, xend = Position, y = y_base, yend = y_base + 0.22),
               color = "#475569", linewidth = 0.8) +
  geom_point(data = ext_phospho, aes(x = Position, y = y_base + 0.22, color = Effect), size = 4.2) +
  geom_text_repel(
    data               = ext_phospho,
    aes(x = Position, y = y_base + 0.24, label = Label, color = Effect),
    direction          = "both",
    nudge_y            = 0.26,
    box.padding        = 0.45,
    point.padding      = 0.25,
    force              = 7,
    force_pull         = 0.35,
    min.segment.length = 0,
    segment.color      = "#64748B",
    segment.size       = 0.45,
    size               = 3.1,
    fontface           = "bold",
    show.legend        = FALSE
  ) +
  scale_color_manual(values = c(
    "Activates NF-kB / SASP"      = "#DC2626",
    "Disables Brake (Locks SASP)" = "#D97706",
    "Immune SASP Shift"           = "#7C3AED",
    "Inhibitory / Rheostat"       = "#2563EB"
  )) +
  scale_y_continuous(
    breaks = 1:5,
    labels = c("STAT1 / STAT3\n(750-770 aa)", "SOCS3\n(225 aa)", "SHP2 (PTPN11)\n(593 aa)",
               "TYK2 / JAK1/2\n(1132-1187 aa)", "gp130 (IL6ST)\n(918 aa)"),
    limits = c(0.65, 5.95)
  ) +
  scale_x_continuous(limits = c(0, 1250), breaks = seq(0, 1200, 200)) +
  theme_bw(base_size = 12) +
  theme(
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = "black", linewidth = 0.9),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey92"),
    panel.grid.minor   = element_blank(),
    legend.position    = "top",
    legend.title       = element_blank(),
    legend.background  = element_rect(fill = "white", color = NA),
    legend.text        = element_text(face = "bold", size = 10, color = "black"),
    axis.text          = element_text(face = "bold", color = "black", size = 10.5),
    axis.title         = element_text(face = "bold", color = "black", size = 11),
    plot.title         = element_text(face = "bold", hjust = 0.5, size = 14, color = "black"),
    plot.subtitle      = element_text(hjust = 0.5, size = 10, color = "grey25")
  ) +
  labs(
    title    = "Extended IL-6 / JAK-STAT Phosphorylation Sites & NF-kB Senescence Impact",
    subtitle = "Activating SASP Switches (Red/Amber) vs. Negative Feedback Rheostats (Blue)",
    x        = "Amino Acid Residue Position (aa)",
    y        = ""
  )

print(p_ext_il6)

out_dir = file.path(Sys.getenv("USERPROFILE"), "Desktop", "Senescence_Results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ext_path = file.path(out_dir, "17_Extended_IL6_JAK_STAT_NFkB_Phospho_Map.png")
ggsave(ext_path, plot = p_ext_il6, width = 13, height = 8.5, dpi = 300, bg = "white")
browseURL(ext_path)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggplot2)
library(dplyr)
library(patchwork)

graphics.off()
set.seed(42)

# 1. 4-Parameter Logistic (4PL) Hill Equation Helper
calc_4pl = function(conc_nM, ic50_nM, hill = 1.3, bottom = 5, top = 100) {
  bottom + (top - bottom) / (1 + (conc_nM / ic50_nM)^hill)
}

# 9-point log10 dose gradient from 1 nM to 10,000 nM (10 uM)
doses_nM = 10^seq(0, 4, length.out = 9)
replicates = 3

# 2. Simulate Viability Dose-Response Data (Proliferating vs. Senescent Cancer Cells)
drug_specs = list(
  list(Drug = "Ruxolitinib (JAK1/2i - Senomorphic)", Prolif_IC50 = 18000, Sen_IC50 = 14000, Prolif_Bot = 82, Sen_Bot = 78),
  list(Drug = "TTI-101 (STAT3i - Senolytic)",        Prolif_IC50 = 3900,  Sen_IC50 = 320,   Prolif_Bot = 15, Sen_Bot = 8),
  list(Drug = "Navitoclax / ABT-263 (BCL-xLi)",      Prolif_IC50 = 2800,  Sen_IC50 = 140,   Prolif_Bot = 12, Sen_Bot = 5)
)

viab_list = list()
idx = 1

for (spec in drug_specs) {
  for (state in c("Proliferating", "Senescent")) {
    ic50_val = ifelse(state == "Proliferating", spec$Prolif_IC50, spec$Sen_IC50)
    bot_val  = ifelse(state == "Proliferating", spec$Prolif_Bot,  spec$Sen_Bot)
    
    for (d in doses_nM) {
      true_v = calc_4pl(d, ic50_nM = ic50_val, hill = 1.35, bottom = bot_val, top = 100)
      obs_v  = pmin(105, pmax(2, rnorm(replicates, mean = true_v, sd = 3.2)))
      
      viab_list[[idx]] = data.frame(
        Drug      = spec$Drug,
        CellState = state,
        Dose_nM   = d,
        Rep       = 1:replicates,
        Viability = obs_v
      )
      idx = idx + 1
    }
  }
}

viab_df = bind_rows(viab_list)
viab_df$Drug = factor(viab_df$Drug, levels = sapply(drug_specs, function(x) x$Drug))

# Compute Mean +/- SD for error bars
viab_sum = viab_df %>%
  group_by(Drug, CellState, Dose_nM) %>%
  summarise(Mean = mean(Viability), SD = sd(Viability), .groups = "drop")

# Generate smooth 4PL curve points for crisp plotting
smooth_doses = 10^seq(0, 4, length.out = 200)
curve_list = list()
c_idx = 1
for (spec in drug_specs) {
  for (state in c("Proliferating", "Senescent")) {
    ic50_val = ifelse(state == "Proliferating", spec$Prolif_IC50, spec$Sen_IC50)
    bot_val  = ifelse(state == "Proliferating", spec$Prolif_Bot,  spec$Sen_Bot)
    curve_list[[c_idx]] = data.frame(
      Drug      = spec$Drug,
      CellState = state,
      Dose_nM   = smooth_doses,
      Viability = calc_4pl(smooth_doses, ic50_nM = ic50_val, hill = 1.35, bottom = bot_val, top = 100)
    )
    c_idx = c_idx + 1
  }
}
curve_df = bind_rows(curve_list)
curve_df$Drug = factor(curve_df$Drug, levels = sapply(drug_specs, function(x) x$Drug))

# IC50 Annotation Labels per Drug Panel
ic50_labels = data.frame(
  Drug  = factor(sapply(drug_specs, function(x) x$Drug), levels = sapply(drug_specs, function(x) x$Drug)),
  Label = c(
    "Senescent IC50: >10,000 nM\nProliferating IC50: >10,000 nM\n(Non-toxic Senomorphic)",
    "Senescent IC50: 320 nM\nProliferating IC50: 3,900 nM\n(12.2x Senolytic Index)",
    "Senescent IC50: 140 nM\nProliferating IC50: 2,800 nM\n(20.0x Senolytic Index)"
  )
)

# 3. Top Panel: Cell Viability IC50 Curves (Proliferating vs. Senescent)
p_viab_ic50 = ggplot() +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey45", linewidth = 0.6) +
  geom_line(data = curve_df, aes(x = Dose_nM, y = Viability, color = CellState), linewidth = 1.2) +
  geom_errorbar(data = viab_sum, aes(x = Dose_nM, ymin = Mean - SD, ymax = Mean + SD, color = CellState), width = 0.12, linewidth = 0.6) +
  geom_point(data = viab_sum, aes(x = Dose_nM, y = Mean, color = CellState), size = 2.5) +
  geom_text(data = ic50_labels, aes(x = 1.2, y = 20, label = Label), hjust = 0, vjust = 0,
            fontface = "bold", size = 3.1, color = "black") +
  facet_wrap(~ Drug, ncol = 3) +
  scale_x_log10(breaks = c(1, 10, 100, 1000, 10000), labels = c("1", "10", "100", "1,000", "10,000")) +
  scale_color_manual(values = c("Proliferating" = "#2563EB", "Senescent" = "#DC2626")) +
  ylim(0, 110) +
  theme_bw(base_size = 11.5) +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = "black", linewidth = 0.85),
    panel.grid.major  = element_line(color = "grey92"),
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill = "#0F172A"),
    strip.text        = element_text(color = "white", face = "bold", size = 10.5),
    legend.position   = "top",
    legend.title      = element_blank(),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text       = element_text(face = "bold", size = 11, color = "black"),
    axis.text         = element_text(face = "bold", color = "black"),
    axis.title        = element_text(face = "bold", color = "black"),
    plot.title        = element_text(face = "bold", hjust = 0.5, size = 14, color = "black")
  ) +
  labs(title = "Cell Viability IC50 Curves: Ruxolitinib (Senomorphic) vs. TTI-101 & Navitoclax (Senolytics)",
       x     = "Drug Concentration (nM, Log10 Scale)",
       y     = "Relative Cell Viability (% of DMSO)")

# 4. Bottom Panel: SASP Suppression EC50 Curves in Senescent Cells (IL-6 / IL-1B Output)
sasp_curve_df = bind_rows(
  data.frame(Drug = "Ruxolitinib (JAK1/2i - Senomorphic)", Dose_nM = smooth_doses, SASP = calc_4pl(smooth_doses, ic50_nM = 18,  hill = 1.4, bottom = 10, top = 100)),
  data.frame(Drug = "TTI-101 (STAT3i - Senolytic)",        Dose_nM = smooth_doses, SASP = calc_4pl(smooth_doses, ic50_nM = 190, hill = 1.3, bottom = 12, top = 100)),
  data.frame(Drug = "Navitoclax / ABT-263 (BCL-xLi)",      Dose_nM = smooth_doses, SASP = calc_4pl(smooth_doses, ic50_nM = 150, hill = 1.4, bottom = 8,  top = 100))
)

p_sasp_ec50 = ggplot(sasp_curve_df, aes(x = Dose_nM, y = SASP, color = Drug)) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey45", linewidth = 0.6) +
  geom_line(linewidth = 1.3) +
  scale_x_log10(breaks = c(1, 10, 100, 1000, 10000), labels = c("1", "10", "100", "1,000", "10,000")) +
  scale_color_manual(values = c(
    "Ruxolitinib (JAK1/2i - Senomorphic)" = "#10B981",
    "TTI-101 (STAT3i - Senolytic)"        = "#DC2626",
    "Navitoclax / ABT-263 (BCL-xLi)"      = "#7C3AED"
  )) +
  annotate("text", x = 22, y = 55, label = "Ruxolitinib SASP EC50 = 18 nM\n(Potent SASP blockade without cell death)",
           color = "#047857", fontface = "bold", size = 3.5, hjust = 0) +
  ylim(0, 108) +
  theme_bw(base_size = 11.5) +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = "black", linewidth = 0.85),
    panel.grid.major  = element_line(color = "grey92"),
    panel.grid.minor  = element_blank(),
    legend.position   = "right",
    legend.title      = element_blank(),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text       = element_text(face = "bold", size = 10, color = "black"),
    axis.text         = element_text(face = "bold", color = "black"),
    axis.title        = element_text(face = "bold", color = "black"),
    plot.title        = element_text(face = "bold", hjust = 0.5, size = 13, color = "black")
  ) +
  labs(title = "Senescent SASP Suppression (IL-6 / IL-1B Secretion EC50 Curves)",
       x     = "Drug Concentration (nM, Log10 Scale)",
       y     = "Secreted SASP Level (% of Control)")

# 5. Combine Both Panels & Save PNG + Raw Dose-Response Data CSV to Desktop
p_ic50_combined = (p_viab_ic50 / p_sasp_ec50) + plot_layout(heights = c(1.25, 1.0))
print(p_ic50_combined)

out_dir = file.path(Sys.getenv("USERPROFILE"), "Desktop", "Senescence_Results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ic50_path = file.path(out_dir, "19_IC50_DoseResponse_Ruxolitinib_TTI101_Navitoclax.png")
ggsave(ic50_path, plot = p_ic50_combined, width = 13.5, height = 9, dpi = 300, bg = "white")

write.csv(viab_df, file.path(out_dir, "19_DoseResponse_Viability_Replicates.csv"), row.names = FALSE)
browseURL(ic50_path)

graphics.off()

# 1. Top Panel: Drug Target Intervention Map across the IL-6 / SHP2 / JAK / STAT3 / SOCS3 Axis
nodes_df = data.frame(
  x     = c(1.0, 2.8, 2.8, 4.8, 4.8, 6.8, 8.8),
  y     = c(2.5, 3.4, 1.6, 3.4, 1.6, 2.5, 2.5),
  label = c("IL-6 / gp130\n(Y759 / Y905)",
            "JAK1 / JAK2\n(p-Y1007/1034)",
            "SHP2 (PTPN11)\n(p-Y542/Y580)",
            "STAT3\n(p-Y705 / p-S727)",
            "IKKb -> NF-kB\n(p-S536 / Ac-K310)",
            "BCL-xL / MCL-1\n& SASP Genes",
            "Cell Fate:\nArrest vs. Death"),
  fill  = c("#F1F5F9", "#DBEAFE", "#DBEAFE", "#FCE7F3", "#FEE2E2", "#FEF3C7", "#F8FAFC"),
  col   = c("#334155", "#1D4ED8", "#1D4ED8", "#BE185D", "#DC2626", "#B45309", "#0F172A")
)

drugs_df = data.frame(
  x     = c(2.8, 2.8, 4.8, 4.8, 1.9),
  y     = c(4.35, 0.65, 4.35, 0.65, 4.35),
  label = c("SENOMORPHIC:\nRuxolitinib / Baricitinib\n(Blocks JAK1/2 -> p-Y705)",
            "SENOMORPHIC:\nSHP099 / TNO155\n(Blocks SHP2 -> IKKb/NF-kB)",
            "SENOLYTIC:\nTTI-101 / Napabucasin\n(Blocks STAT3 + BCL-xL)",
            "SENOMORPHIC:\nRapamycin / Torin-1\n(Blocks mTOR -> p-S727)",
            "SENOLYTIC:\nDasatinib (+ Quercetin)\n(Protects SOCS3 Brake)"),
  fill  = c("#EFF6FF", "#EFF6FF", "#FFF1F2", "#EFF6FF", "#FFF1F2"),
  col   = c("#2563EB", "#2563EB", "#E11D48", "#2563EB", "#E11D48")
)

p_drug_map = ggplot() +
  # Pathway signaling arrows
  geom_segment(aes(x = 1.55, y = 2.75, xend = 2.25, yend = 3.25), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.0, color = "#475569") +
  geom_segment(aes(x = 1.55, y = 2.25, xend = 2.25, yend = 1.75), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.0, color = "#475569") +
  geom_segment(aes(x = 3.35, y = 3.40, xend = 4.20, yend = 3.40), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.0, color = "#475569") +
  geom_segment(aes(x = 3.35, y = 1.60, xend = 4.20, yend = 1.60), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.0, color = "#475569") +
  geom_segment(aes(x = 4.80, y = 2.95, xend = 4.80, yend = 2.05), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.0, color = "#BE185D", linetype = "dashed") +
  geom_segment(aes(x = 5.40, y = 3.25, xend = 6.20, yend = 2.75), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.0, color = "#475569") +
  geom_segment(aes(x = 5.40, y = 1.75, xend = 6.20, yend = 2.25), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.0, color = "#475569") +
  geom_segment(aes(x = 7.40, y = 2.50, xend = 8.20, yend = 2.50), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.0, color = "#475569") +
  # Drug inhibition T-bars (blunt red/blue lines)
  geom_segment(aes(x = 2.8, y = 3.95, xend = 2.8, yend = 3.72), linewidth = 1.3, color = "#2563EB") +
  geom_segment(aes(x = 2.8, y = 1.05, xend = 2.8, yend = 1.28), linewidth = 1.3, color = "#2563EB") +
  geom_segment(aes(x = 4.8, y = 3.95, xend = 4.8, yend = 3.72), linewidth = 1.3, color = "#E11D48") +
  geom_segment(aes(x = 4.8, y = 1.05, xend = 4.8, yend = 1.28), linewidth = 1.3, color = "#2563EB") +
  geom_segment(aes(x = 2.1, y = 3.95, xend = 2.45, yend = 3.65), linewidth = 1.3, color = "#E11D48") +
  # Pathway Nodes
  geom_label(data = nodes_df, aes(x = x, y = y, label = label, fill = fill, color = col),
             fontface = "bold", size = 3.2, label.padding = unit(0.45, "lines"), label.size = 1.1) +
  # Drug Intervention Callouts
  geom_label(data = drugs_df, aes(x = x, y = y, label = label, fill = fill, color = col),
             fontface = "bold", size = 2.9, label.padding = unit(0.38, "lines"), label.size = 1.0) +
  scale_fill_identity() +
  scale_color_identity() +
  xlim(0.3, 9.5) + ylim(0.15, 4.85) +
  theme_void(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = "black", linewidth = 0.85),
    plot.title      = element_text(face = "bold", size = 13.5, hjust = 0.5, color = "black", margin = margin(t = 8, b = 3)),
    plot.subtitle   = element_text(size = 10, hjust = 0.5, color = "grey30", margin = margin(b = 5))
  ) +
  labs(title    = "Senomorphic (Blue) vs. Senolytic (Red) Drug Targets on the IL-6 / SHP2 / JAK / STAT3 Axis",
       subtitle = "Senomorphics Silence SASP While Maintaining Growth Arrest | Senolytics Collapse BCL-xL/STAT3 to Kill Senescent Cells")

# 2. Bottom Panel: Pharmacological Profile Comparison (SASP Output vs. Senescent Cell Survival)
drug_comp = data.frame(
  Treatment = factor(c("Vehicle (SASP Control)", "SHP099 (SHP2i)", "Ruxolitinib (JAK1/2i)",
                       "Rapamycin (mTORC1i)", "TTI-101 (STAT3i)", "Dasatinib + Quercetin"),
                     levels = c("Vehicle (SASP Control)", "SHP099 (SHP2i)", "Ruxolitinib (JAK1/2i)",
                                "Rapamycin (mTORC1i)", "TTI-101 (STAT3i)", "Dasatinib + Quercetin")),
  Strategy  = c("Untreated", "Senomorphic", "Senomorphic", "Senomorphic", "Senolytic", "Senolytic"),
  SASP_Level = c(100, 38, 22, 28, 18, 12),
  Cell_Viability = c(96, 92, 90, 91, 34, 24)
)

drug_long = tidyr::pivot_longer(drug_comp, cols = c("SASP_Level", "Cell_Viability"),
                                names_to = "Metric", values_to = "Percent")
drug_long$Metric = ifelse(drug_long$Metric == "SASP_Level",
                          "SASP Secretion (IL-6 / IL-1B %)",
                          "Senescent Tumor Cell Viability (%)")

p_drug_bar = ggplot(drug_long, aes(x = Treatment, y = Percent, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(0.75), width = 0.65, color = "black", linewidth = 0.4) +
  geom_text(aes(label = paste0(Percent, "%")), position = position_dodge(0.75), vjust = -0.4, fontface = "bold", size = 3.2) +
  scale_fill_manual(values = c("SASP Secretion (IL-6 / IL-1B %)" = "#F59E0B",
                               "Senescent Tumor Cell Viability (%)" = "#2563EB")) +
  ylim(0, 115) +
  theme_bw(base_size = 11.5) +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = "black", linewidth = 0.85),
    panel.grid.major.x = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "top",
    legend.title      = element_blank(),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text       = element_text(face = "bold", size = 10.5, color = "black"),
    axis.text.x       = element_text(angle = 20, hjust = 1, face = "bold", color = "black", size = 10),
    axis.text.y       = element_text(face = "bold", color = "black"),
    axis.title        = element_text(face = "bold", color = "black"),
    plot.title        = element_text(face = "bold", hjust = 0.5, size = 13, color = "black")
  ) +
  labs(title = "Functional Readout: Senomorphics Drop SASP Only vs. Senolytics Eliminate Both SASP & Viability",
       x     = "",
       y     = "Relative Level (% of Vehicle Control)")

# 3. Combine & Save on Pure White Background
p_senotherapeutics = (p_drug_map / p_drug_bar) + plot_layout(heights = c(1.25, 1.0))
print(p_senotherapeutics)

out_dir = file.path(Sys.getenv("USERPROFILE"), "Desktop", "Senescence_Results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

st_path = file.path(out_dir, "18_Senomorphic_vs_Senolytic_Targets_IL6_STAT3.png")
ggsave(st_path, plot = p_senotherapeutics, width = 13, height = 9.5, dpi = 300, bg = "white")
browseURL(st_path)
# ==============================================================================
# Complete File Inventory Check: Desktop -> Senescence_Results
# ==============================================================================

out_dir = file.path(Sys.getenv("USERPROFILE"), "Desktop", "Senescence_Results")

if (!dir.exists(out_dir)) {
  stop("Folder not found at: ", out_dir)
}

# 1. Gather all files and metadata
all_files = list.files(out_dir, full.names = TRUE)
f_info    = file.info(all_files)

inventory_df = data.frame(
  File_No      = seq_along(all_files),
  File_Name    = basename(all_files),
  Type         = toupper(tools::file_ext(all_files)),
  Size_KB      = round(f_info$size / 1024, 1),
  Last_Updated = format(f_info$mtime, "%Y-%m-%d %H:%M:%S"),
  stringsAsFactors = FALSE
)

# 2. Print Summary Counts & Clean Console Table
n_png = sum(inventory_df$Type == "PNG")
n_csv = sum(inventory_df$Type == "CSV")

cat("\n====================================================================\n")
cat("  SENESCENCE & SASP PIPELINE - DESKTOP RESULTS INVENTORY\n")
cat("  Folder:", normalizePath(out_dir), "\n")
cat(sprintf("  Total Files: %d  (%d High-Res PNG Figures | %d CSV Data Tables)\n",
            nrow(inventory_df), n_png, n_csv))
cat("====================================================================\n\n")

print(inventory_df, row.names = FALSE, right = FALSE)

# 3. Save Master Manifest CSV & Open Folder in Windows File Explorer
manifest_path = file.path(out_dir, "00_Project_File_Inventory.csv")
write.csv(inventory_df, manifest_path, row.names = FALSE)

cat("\nSaved master file index to:", manifest_path, "\n")
browseURL(out_dir)
library(r3dmol)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)

graphics.off()

# ==============================================================================
# PART A: INTERACTIVE 3D CRYSTAL STRUCTURES IN RSTUDIO VIEWER (r3dmol)
# ==============================================================================

# 1. Interactive 3D Structure 1: Hexameric IL-6 / IL-6Ra / gp130 Complex (PDB: 1P9M)
#    Highlights Site I (Phe74/Arg179), Site II (Tyr31/Val121), and Site III (Trp157)
view_1p9m = r3dmol() %>%
  m_add_model(data = m_fetch_pdb("1P9M"), format = "pdb") %>%
  m_set_style(style = m_style_cartoon(color = "#CBD5E1", opacity = 0.85)) %>%
  # Color Chains: IL-6 (Red), IL-6Ra (Blue), gp130 (Teal)
  m_set_style(sel = m_sel(chain = "A"), style = m_style_cartoon(color = "#EF4444")) %>%
  m_set_style(sel = m_sel(chain = "B"), style = m_style_cartoon(color = "#3B82F6")) %>%
  m_set_style(sel = m_sel(chain = "C"), style = m_style_cartoon(color = "#10B981")) %>%
  # Highlight Hotspots (Site I: 74,179 | Site II: 31,121 | Site III: 157) as Gold Spheres
  m_set_style(sel = m_sel(chain = "A", resi = c(31, 74, 121, 157, 179)),
              style = m_style_sphere(color = "#F59E0B", scale = 0.9)) %>%
  m_add_label(text = "Site III Hotspot (Trp157)", sel = m_sel(chain = "A", resi = 157),
              style = m_style_label(backgroundColor = "#0F172A", fontColor = "white", fontSize = 11)) %>%
  m_add_label(text = "Site I Hotspot (Arg179)", sel = m_sel(chain = "A", resi = 179),
              style = m_style_label(backgroundColor = "#DC2626", fontColor = "white", fontSize = 11)) %>%
  m_zoom_to() %>%
  m_spin(speed = 0.4)

# Display 1P9M in RStudio Viewer Pane (Rotate/Zoom with mouse)
print(view_1p9m)

# 2. Interactive 3D Structure 2: STAT3 SH2 Domain & Tyr705 Docking Pocket (PDB: 1BG1)
#    Highlights TTI-101 / Stattic SH2 Pocket (Lys591, Arg609, Ser611, Glu612) & p-Tyr705
view_stat3 = r3dmol() %>%
  m_add_model(data = m_fetch_pdb("1BG1"), format = "pdb") %>%
  m_set_style(style = m_style_cartoon(color = "#94A3B8", opacity = 0.85)) %>%
  m_set_style(sel = m_sel(resi = 580:680), style = m_style_cartoon(color = "#8B5CF6")) %>%
  # Highlight SH2 Docking Pocket Residues (Lys591, Arg609, Ser611, Glu612) in Cyan Sticks/Spheres
  m_set_style(sel = m_sel(resi = c(591, 609, 611, 612)),
              style = m_style_stick(color = "#06B6D4", radius = 0.35)) %>%
  # Highlight Phosphorylated Tyr705 in Bright Red Sphere
  m_set_style(sel = m_sel(resi = 705),
              style = m_style_sphere(color = "#DC2626", scale = 1.1)) %>%
  m_add_label(text = "STAT3 SH2 Docking Pocket (Arg609 / Lys591)", sel = m_sel(resi = 609),
              style = m_style_label(backgroundColor = "#0F172A", fontColor = "white", fontSize = 11)) %>%
  m_add_label(text = "p-Tyr705 Dimerization Switch", sel = m_sel(resi = 705),
              style = m_style_label(backgroundColor = "#DC2626", fontColor = "white", fontSize = 11)) %>%
  m_zoom_to(sel = m_sel(resi = 585:710))

# ==============================================================================
# PART B: HIGH-RES 3D-STYLED PPI & ACTIVE DOCKING SITE MAP (SAVED TO DESKTOP)
# ==============================================================================

# Build 3D Structural Architecture & Docking Pocket Diagram across Plasma Membrane
arch_nodes = data.frame(
  x     = c(2.0, 2.0, 2.0, 5.2, 5.2, 8.2, 8.2),
  y     = c(4.3, 3.1, 1.8, 3.1, 1.8, 3.1, 1.8),
  label = c(
    "IL-6 Cytokine (PDB: 1P9M)\nSite I (Phe74/Arg179) -> IL-6Ra\nSite III (Trp157) -> gp130 D1",
    "gp130 + IL-6Ra Hexamer (2:2:2)\nSite II Pocket: gp130 Phe191/Val252\n[Drug: Bazedoxifene / Olamkicept]",
    "gp130 Box1/Box2 Motif\nPro652 / Pro655 / Trp666\nDocks JAK1/2 FERM-SH2 Domain",
    "JAK1/2 JH2 Pseudokinase\nAllosteric Pocket (Val658 / Leu95)\n[Drug: Deucravacitinib / BMS-986165]",
    "JAK1/2 JH1 Active Kinase Cleft\nHinge: Glu957/Leu959 | Gate: Met956\nAct. Loop: p-Y1007/1008 [Ruxolitinib]",
    "STAT3 SH2 Docking Pocket (PDB: 1BG1)\np-Tyr Pocket: Arg609 / Lys591 / Ser611\nSub-pocket: Ile634 [TTI-101 / Stattic]",
    "STAT3 Reciprocal Dimer & TAD\np-Tyr705 (SH2 Lock) + p-Ser727 (p300)\nBinds SASP Promoters (IL6 / BCL2L1)"
  ),
  zone  = c("Extracellular", "Extracellular", "Intracellular", "Intracellular", "Intracellular", "Nuclear / Cytosol", "Nuclear / Cytosol"),
  fill  = c("#FEE2E2", "#DBEAFE", "#E0E7FF", "#FEF3C7", "#FCE7F3", "#CCFBF1", "#F3E8FF"),
  col   = c("#DC2626", "#1D4ED8", "#4338CA", "#B45309", "#BE185D", "#0F766E", "#6D28D9")
)

p_ppi_3d_map = ggplot() +
  # Plasma Membrane Bilayer Band
  annotate("rect", xmin = 0.3, xmax = 9.8, ymin = 2.35, ymax = 2.55, fill = "#94A3B8", alpha = 0.45) +
  annotate("text", x = 9.65, y = 2.45, label = "PLASMA MEMBRANE BILAYER", hjust = 1, fontface = "bold", size = 3.0, color = "#1E293B") +
  # PPI & Catalytic Arrows
  geom_segment(aes(x = 2.0, y = 3.85, xend = 2.0, yend = 3.50), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.1, color = "#334155") +
  geom_segment(aes(x = 2.0, y = 2.70, xend = 2.0, yend = 2.20), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.1, color = "#334155") +
  geom_segment(aes(x = 3.35, y = 1.80, xend = 3.90, yend = 1.80), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.1, color = "#334155") +
  geom_segment(aes(x = 5.20, y = 2.70, xend = 5.20, yend = 2.20), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.1, color = "#B45309", linetype = "dashed") +
  geom_segment(aes(x = 6.55, y = 1.80, xend = 7.05, yend = 3.00), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.1, color = "#BE185D") +
  geom_segment(aes(x = 8.20, y = 2.70, xend = 8.20, yend = 2.20), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.1, color = "#0F766E") +
  # Structural Domain & Docking Cards
  geom_label(data = arch_nodes, aes(x = x, y = y, label = label, fill = fill, color = col),
             fontface = "bold", size = 3.0, label.padding = unit(0.45, "lines"), label.size = 1.1) +
  scale_fill_identity() +
  scale_color_identity() +
  xlim(0.3, 9.8) + ylim(1.1, 4.9) +
  theme_void(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = "black", linewidth = 0.9),
    plot.title      = element_text(face = "bold", size = 14, hjust = 0.5, color = "black", margin = margin(t = 10, b = 4)),
    plot.subtitle   = element_text(size = 10, hjust = 0.5, color = "grey30", margin = margin(b = 6))
  ) +
  labs(title    = "3D Structural PPI Interfaces, Active Catalytic Clefts & Small-Molecule Docking Sites",
       subtitle = "Mapped from Crystal Structures PDB: 1P9M (IL-6 Hexamer), PDB: 7T6F (JAK1 Dimer), and PDB: 1BG1 (STAT3 SH2)")

print(p_ppi_3d_map)

out_dir = file.path(Sys.getenv("USERPROFILE"), "Desktop", "Senescence_Results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ppi_path = file.path(out_dir, "20_IL6_JAK_STAT_3D_Structure_and_Docking_Sites.png")
ggsave(ppi_path, plot = p_ppi_3d_map, width = 13.5, height = 7.5, dpi = 300, bg = "white")
browseURL(ppi_path)
library(r3dmol)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)

graphics.off()

# ==============================================================================
# PART A: INTERACTIVE 3D CRYSTAL STRUCTURES IN RSTUDIO VIEWER (r3dmol)
# ==============================================================================

# 1. Interactive 3D Structure 1: Hexameric IL-6 / IL-6Ra / gp130 Complex (PDB: 1P9M)
#    Highlights Site I (Phe74/Arg179), Site II (Tyr31/Val121), and Site III (Trp157)
view_1p9m = r3dmol() %>%
  m_add_model(data = m_fetch_pdb("1P9M"), format = "pdb") %>%
  m_set_style(style = m_style_cartoon(color = "#CBD5E1", opacity = 0.85)) %>%
  # Color Chains: IL-6 (Red), IL-6Ra (Blue), gp130 (Teal)
  m_set_style(sel = m_sel(chain = "A"), style = m_style_cartoon(color = "#EF4444")) %>%
  m_set_style(sel = m_sel(chain = "B"), style = m_style_cartoon(color = "#3B82F6")) %>%
  m_set_style(sel = m_sel(chain = "C"), style = m_style_cartoon(color = "#10B981")) %>%
  # Highlight Hotspots (Site I: 74,179 | Site II: 31,121 | Site III: 157) as Gold Spheres
  m_set_style(sel = m_sel(chain = "A", resi = c(31, 74, 121, 157, 179)),
              style = m_style_sphere(color = "#F59E0B", scale = 0.9)) %>%
  m_add_label(text = "Site III Hotspot (Trp157)", sel = m_sel(chain = "A", resi = 157),
              style = m_style_label(backgroundColor = "#0F172A", fontColor = "white", fontSize = 11)) %>%
  m_add_label(text = "Site I Hotspot (Arg179)", sel = m_sel(chain = "A", resi = 179),
              style = m_style_label(backgroundColor = "#DC2626", fontColor = "white", fontSize = 11)) %>%
  m_zoom_to() %>%
  m_spin(speed = 0.4)

# Display 1P9M in RStudio Viewer Pane (Rotate/Zoom with mouse)
print(view_1p9m)

# 2. Interactive 3D Structure 2: STAT3 SH2 Domain & Tyr705 Docking Pocket (PDB: 1BG1)
#    Highlights TTI-101 / Stattic SH2 Pocket (Lys591, Arg609, Ser611, Glu612) & p-Tyr705
view_stat3 = r3dmol() %>%
  m_add_model(data = m_fetch_pdb("1BG1"), format = "pdb") %>%
  m_set_style(style = m_style_cartoon(color = "#94A3B8", opacity = 0.85)) %>%
  m_set_style(sel = m_sel(resi = 580:680), style = m_style_cartoon(color = "#8B5CF6")) %>%
  # Highlight SH2 Docking Pocket Residues (Lys591, Arg609, Ser611, Glu612) in Cyan Sticks/Spheres
  m_set_style(sel = m_sel(resi = c(591, 609, 611, 612)),
              style = m_style_stick(color = "#06B6D4", radius = 0.35)) %>%
  # Highlight Phosphorylated Tyr705 in Bright Red Sphere
  m_set_style(sel = m_sel(resi = 705),
              style = m_style_sphere(color = "#DC2626", scale = 1.1)) %>%
  m_add_label(text = "STAT3 SH2 Docking Pocket (Arg609 / Lys591)", sel = m_sel(resi = 609),
              style = m_style_label(backgroundColor = "#0F172A", fontColor = "white", fontSize = 11)) %>%
  m_add_label(text = "p-Tyr705 Dimerization Switch", sel = m_sel(resi = 705),
              style = m_style_label(backgroundColor = "#DC2626", fontColor = "white", fontSize = 11)) %>%
  m_zoom_to(sel = m_sel(resi = 585:710))

# ==============================================================================
# PART B: HIGH-RES 3D-STYLED PPI & ACTIVE DOCKING SITE MAP (SAVED TO DESKTOP)
# ==============================================================================

# Build 3D Structural Architecture & Docking Pocket Diagram across Plasma Membrane
arch_nodes = data.frame(
  x     = c(2.0, 2.0, 2.0, 5.2, 5.2, 8.2, 8.2),
  y     = c(4.3, 3.1, 1.8, 3.1, 1.8, 3.1, 1.8),
  label = c(
    "IL-6 Cytokine (PDB: 1P9M)\nSite I (Phe74/Arg179) -> IL-6Ra\nSite III (Trp157) -> gp130 D1",
    "gp130 + IL-6Ra Hexamer (2:2:2)\nSite II Pocket: gp130 Phe191/Val252\n[Drug: Bazedoxifene / Olamkicept]",
    "gp130 Box1/Box2 Motif\nPro652 / Pro655 / Trp666\nDocks JAK1/2 FERM-SH2 Domain",
    "JAK1/2 JH2 Pseudokinase\nAllosteric Pocket (Val658 / Leu95)\n[Drug: Deucravacitinib / BMS-986165]",
    "JAK1/2 JH1 Active Kinase Cleft\nHinge: Glu957/Leu959 | Gate: Met956\nAct. Loop: p-Y1007/1008 [Ruxolitinib]",
    "STAT3 SH2 Docking Pocket (PDB: 1BG1)\np-Tyr Pocket: Arg609 / Lys591 / Ser611\nSub-pocket: Ile634 [TTI-101 / Stattic]",
    "STAT3 Reciprocal Dimer & TAD\np-Tyr705 (SH2 Lock) + p-Ser727 (p300)\nBinds SASP Promoters (IL6 / BCL2L1)"
  ),
  zone  = c("Extracellular", "Extracellular", "Intracellular", "Intracellular", "Intracellular", "Nuclear / Cytosol", "Nuclear / Cytosol"),
  fill  = c("#FEE2E2", "#DBEAFE", "#E0E7FF", "#FEF3C7", "#FCE7F3", "#CCFBF1", "#F3E8FF"),
  col   = c("#DC2626", "#1D4ED8", "#4338CA", "#B45309", "#BE185D", "#0F766E", "#6D28D9")
)

p_ppi_3d_map = ggplot() +
  # Plasma Membrane Bilayer Band
  annotate("rect", xmin = 0.3, xmax = 9.8, ymin = 2.35, ymax = 2.55, fill = "#94A3B8", alpha = 0.45) +
  annotate("text", x = 9.65, y = 2.45, label = "PLASMA MEMBRANE BILAYER", hjust = 1, fontface = "bold", size = 3.0, color = "#1E293B") +
  # PPI & Catalytic Arrows
  geom_segment(aes(x = 2.0, y = 3.85, xend = 2.0, yend = 3.50), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.1, color = "#334155") +
  geom_segment(aes(x = 2.0, y = 2.70, xend = 2.0, yend = 2.20), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.1, color = "#334155") +
  geom_segment(aes(x = 3.35, y = 1.80, xend = 3.90, yend = 1.80), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.1, color = "#334155") +
  geom_segment(aes(x = 5.20, y = 2.70, xend = 5.20, yend = 2.20), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.1, color = "#B45309", linetype = "dashed") +
  geom_segment(aes(x = 6.55, y = 1.80, xend = 7.05, yend = 3.00), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.1, color = "#BE185D") +
  geom_segment(aes(x = 8.20, y = 2.70, xend = 8.20, yend = 2.20), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.1, color = "#0F766E") +
  # Structural Domain & Docking Cards
  geom_label(data = arch_nodes, aes(x = x, y = y, label = label, fill = fill, color = col),
             fontface = "bold", size = 3.0, label.padding = unit(0.45, "lines"), label.size = 1.1) +
  scale_fill_identity() +
  scale_color_identity() +
  xlim(0.3, 9.8) + ylim(1.1, 4.9) +
  theme_void(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = "black", linewidth = 0.9),
    plot.title      = element_text(face = "bold", size = 14, hjust = 0.5, color = "black", margin = margin(t = 10, b = 4)),
    plot.subtitle   = element_text(size = 10, hjust = 0.5, color = "grey30", margin = margin(b = 6))
  ) +
  labs(title    = "3D Structural PPI Interfaces, Active Catalytic Clefts & Small-Molecule Docking Sites",
       subtitle = "Mapped from Crystal Structures PDB: 1P9M (IL-6 Hexamer), PDB: 7T6F (JAK1 Dimer), and PDB: 1BG1 (STAT3 SH2)")

print(p_ppi_3d_map)

out_dir = file.path(Sys.getenv("USERPROFILE"), "Desktop", "Senescence_Results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ppi_path = file.path(out_dir, "20_IL6_JAK_STAT_3D_Structure_and_Docking_Sites.png")
ggsave(ppi_path, plot = p_ppi_3d_map, width = 13.5, height = 7.5, dpi = 300, bg = "white")
browseURL(ppi_path)
library(ggplot2)
library(ggforce)

graphics.off()

# Helper to draw smooth rounded structural lobes
make_lobe = function(cx, cy, rx, ry, angle_deg, n = 120) {
  t = seq(0, 2 * pi, length.out = n)
  rad = angle_deg * pi / 180
  x = cx + rx * cos(t) * cos(rad) - ry * sin(t) * sin(rad)
  y = cy + rx * cos(t) * sin(rad) + ry * sin(t) * cos(rad)
  data.frame(x = x, y = y)
}

n_lobe = make_lobe(cx = 5.2, cy = 7.3, rx = 3.1, ry = 1.35, angle_deg = 6)
c_lobe = make_lobe(cx = 5.4, cy = 2.6, rx = 3.4, ry = 1.55, angle_deg = -4)

# Structural Elements & Active-Site Residues
motifs_df = data.frame(
  x     = c(1.65, 3.25, 4.95, 6.65, 5.10, 8.65, 8.85),
  y     = c(5.05, 6.20, 6.45, 5.95, 4.15, 4.75, 2.85),
  label = c(
    "HINGE REGION\nJAK1: Glu957 / Leu959\nJAK2: Glu930 / Leu932\n(H-Bond Anchor for Ruxolitinib)",
    "GATEKEEPER\nJAK1: Met956\nJAK2: Met929\n(Hydrophobic Back Pocket)",
    "P-LOOP (Gly-Rich)\nJAK1: Gly882-Gly887\n(Clamps ATP Phosphates)",
    "CATALYTIC LYSINE\nJAK1: Lys908 | JAK2: Lys882\n(Salt Bridge w/ aC-Helix Glu925)",
    "DFG MOTIF (Mg2+ Switch)\nJAK1: Asp994-Phe995-Gly996\nJAK2: Asp967-Phe968-Gly969",
    "ACTIVATION LOOP (Open)\nJAK1: p-Tyr1034 / p-Tyr1035\nJAK2: p-Tyr1007 / p-Tyr1008",
    "STAT3 SUBSTRATE GROOVE\nBinds STAT3 Tyr705 Tail\nfor Phosphorylation"
  ),
  fill  = c("#CFFAFE", "#FEF9C3", "#E0E7FF", "#F3E8FF", "#FCE7F3", "#FEE2E2", "#DCFCE7"),
  col   = c("#0E7490", "#A16207", "#3730A3", "#6B21A8", "#BE185D", "#DC2626", "#15803D")
)

# Central Ligand / Inhibitor Core (Ruxolitinib / ATP Pyrrolopyrimidine Scaffold)
ligand_box = data.frame(
  xmin = 3.55, xmax = 6.15, ymin = 4.75, ymax = 5.45
)

p_jh1_cleft = ggplot() +
  # 1. Draw N-Lobe (Top) and C-Lobe (Bottom)
  geom_polygon(data = n_lobe, aes(x = x, y = y), fill = "#DBEAFE", color = "#1D4ED8", linewidth = 1.4, alpha = 0.85) +
  geom_polygon(data = c_lobe, aes(x = x, y = y), fill = "#D1FAE5", color = "#047857", linewidth = 1.4, alpha = 0.85) +
  
  # Lobe Headers
  annotate("text", x = 5.2, y = 8.15, label = "N-TERMINAL LOBE (Beta-Sheet Rich + aC-Helix)",
           fontface = "bold", size = 4.0, color = "#1E3A8A") +
  annotate("text", x = 5.3, y = 1.65, label = "C-TERMINAL LOBE (Alpha-Helical Catalytic Core)",
           fontface = "bold", size = 4.0, color = "#065F46") +
  
  # 2. Hinge Backbone Bridge connecting N-Lobe & C-Lobe on Left
  geom_curve(aes(x = 2.35, y = 6.65, xend = 2.35, yend = 3.45),
             curvature = 0.35, linewidth = 4.5, color = "#06B6D4") +
  
  # 3. Activation Loop Ribbon swinging open on Right
  geom_curve(aes(x = 7.45, y = 3.65, xend = 8.10, yend = 5.65),
             curvature = -0.45, linewidth = 3.8, color = "#EF4444") +
  
  # 4. Central ATP / Ruxolitinib Inhibitor Docking Pocket
  geom_rect(data = ligand_box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "#FEF3C7", color = "#D97706", linewidth = 1.3, linetype = "solid") +
  annotate("text", x = 4.85, y = 5.10,
           label = "ATP / RUXOLITINIB BINDING CLEFT\n(Pyrrolopyrimidine Core + Cyclopentyl Ring)",
           fontface = "bold", size = 3.2, color = "#92400E") +
  
  # 5. Bidentate Hydrogen Bonds (Hinge <-> Ligand)
  geom_segment(aes(x = 2.15, y = 5.25, xend = 3.55, yend = 5.25),
               linetype = "dashed", linewidth = 1.1, color = "#0E7490") +
  geom_segment(aes(x = 2.15, y = 4.95, xend = 3.55, yend = 4.95),
               linetype = "dashed", linewidth = 1.1, color = "#0E7490") +
  annotate("text", x = 2.85, y = 5.42, label = "2x H-Bonds", fontface = "bold", size = 2.9, color = "#0E7490") +
  
  # 6. Phospho-Tyrosine Spheres on Activation Loop
  geom_point(aes(x = c(7.92, 8.05), y = c(4.95, 5.32)), size = 5.5, color = "#DC2626") +
  annotate("text", x = 7.62, y = 5.15, label = "P", fontface = "bold", size = 3.2, color = "white") +
  
  # 7. Leader Lines from Active Residues to Callout Cards
  geom_segment(aes(x = 3.45, y = 5.80, xend = 3.65, yend = 5.45), color = "#A16207", linewidth = 0.8) +
  geom_segment(aes(x = 4.95, y = 6.00, xend = 4.95, yend = 5.45), color = "#3730A3", linewidth = 0.8) +
  geom_segment(aes(x = 6.35, y = 5.55, xend = 5.95, yend = 5.35), color = "#6B21A8", linewidth = 0.8) +
  geom_segment(aes(x = 5.10, y = 4.55, xend = 5.10, yend = 4.75), color = "#BE185D", linewidth = 0.8) +
  
  # 8. Residue & Motif Callout Boxes
  geom_label(data = motifs_df, aes(x = x, y = y, label = label, fill = fill, color = col),
             fontface = "bold", size = 2.85, label.padding = unit(0.42, "lines"), label.size = 1.0) +
  
  scale_fill_identity() +
  scale_color_identity() +
  xlim(0.2, 10.3) + ylim(0.9, 9.0) +
  theme_void(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = "black", linewidth = 1.0),
    plot.title      = element_text(face = "bold", size = 14.5, hjust = 0.5, color = "black", margin = margin(t = 10, b = 4)),
    plot.subtitle   = element_text(size = 10.5, hjust = 0.5, color = "grey25", margin = margin(b = 8))
  ) +
  labs(
    title    = "2D Structural Architecture of the JAK1 / JAK2 JH1 Catalytic Kinase Cleft",
    subtitle = "Active-Site Residues, Ruxolitinib Hinge Hydrogen-Bonding & Activation-Loop Phosphorylation Switch (PDB: 7T6F / 4E4L)"
  )

print(p_jh1_cleft)

# Save High-Res PNG to Desktop -> Senescence_Results
out_dir = file.path(Sys.getenv("USERPROFILE"), "Desktop", "Senescence_Results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cleft_path = file.path(out_dir, "22_JAK1_JAK2_JH1_Catalytic_Cleft_Diagram.png")
ggsave(cleft_path, plot = p_jh1_cleft, width = 13.5, height = 8.5, dpi = 300, bg = "white")

browseURL(cleft_path)
library(r3dmol)
library(ggplot2)
library(dplyr)

graphics.off()

# ==============================================================================
# PART A: INTERACTIVE 3D CLOSE-UP OF IL-6 BINDING SITES I, II & III (PDB: 1P9M)
# ==============================================================================
view_il6_sites = r3dmol() %>%
  m_add_model(data = m_fetch_pdb("1P9M"), format = "pdb") %>%
  # Hide other chains slightly or keep translucent so IL-6 (Chain A) binding pockets stand out
  m_set_style(sel = m_sel(chain = c("B", "C", "D", "E", "F")),
              style = m_style_cartoon(color = "#E2E8F0", opacity = 0.25)) %>%
  # Color IL-6 Four-Helix Bundle (Chain A) in clean Slate Teal
  m_set_style(sel = m_sel(chain = "A"), style = m_style_cartoon(color = "#0D9488")) %>%
  
  # SITE I (IL-6Ra Binding Pocket: Phe74, Phe78, Leu178, Arg179, Arg182) -> Crimson Spheres
  m_set_style(sel = m_sel(chain = "A", resi = c(74, 78, 178, 179, 182)),
              style = m_style_sphere(color = "#DC2626", scale = 0.85)) %>%
  
  # SITE II (First gp130 D2-D3 Binding Pocket: Tyr31, Gly35, Ser118, Val121, Gln124) -> Royal Blue Spheres
  m_set_style(sel = m_sel(chain = "A", resi = c(31, 35, 118, 121, 124)),
              style = m_style_sphere(color = "#2563EB", scale = 0.85)) %>%
  
  # SITE III (Second gp130 Ig-D1 Hexamer Hotspot: Trp157, Asp160, Thr162) -> Amber Gold Spheres
  m_set_style(sel = m_sel(chain = "A", resi = c(157, 160, 162)),
              style = m_style_sphere(color = "#D97706", scale = 0.95)) %>%
  
  # Add 3D Floating Labels for Each Binding Site
  m_add_label(text = "Site I: IL-6Ra Pocket (Phe74 / Arg179 / Arg182)",
              sel = m_sel(chain = "A", resi = 179),
              style = m_style_label(backgroundColor = "#DC2626", fontColor = "white", fontSize = 11)) %>%
  m_add_label(text = "Site II: gp130 D2-D3 Pocket (Tyr31 / Val121)",
              sel = m_sel(chain = "A", resi = 121),
              style = m_style_label(backgroundColor = "#2563EB", fontColor = "white", fontSize = 11)) %>%
  m_add_label(text = "Site III: gp130 D1 Hexamer Hotspot (Trp157 / Asp160)",
              sel = m_sel(chain = "A", resi = 157),
              style = m_style_label(backgroundColor = "#D97706", fontColor = "white", fontSize = 11)) %>%
  
  m_zoom_to(sel = m_sel(chain = "A")) %>%
  m_spin(speed = 0.4)

# Open Interactive 3D Model in RStudio Viewer
print(view_il6_sites)

# ==============================================================================
# PART B: HIGH-RES 2D STRUCTURAL DIAGRAM OF IL-6 BINDING SITES I, II & III
# ==============================================================================
p_il6_map = ggplot() +
  # Central IL-6 Four-Helix Bundle Core
  annotate("rect", xmin = 3.6, xmax = 6.4, ymin = 2.4, ymax = 5.6,
           fill = "#CCFBF1", color = "#0F766E", linewidth = 1.6) +
  annotate("text", x = 5.0, y = 4.35, label = "HUMAN IL-6\nFour-Helix Bundle\n(Helices A, B, C, D)",
           fontface = "bold", size = 4.3, color = "#115E59") +
  annotate("text", x = 5.0, y = 3.25, label = "184 aa Mature Cytokine\n(PDB: 1P9M Chain A)",
           fontface = "bold", size = 3.3, color = "#0F766E") +
  
  # SITE I (Left): IL-6 <-> IL-6Ra
  annotate("label", x = 1.65, y = 4.65,
           label = "SITE I EPITOPE\n(Helix A/D + AB-Loop)\nHotspots: Phe74, Phe78,\nLeu178, Arg179, Arg182",
           fill = "#FEE2E2", color = "#991B1B", fontface = "bold", size = 3.1, label.size = 1.1) +
  annotate("label", x = 1.65, y = 2.75,
           label = "Partner: IL-6Ra (CD126)\nPocket: Phe229, Phe279, Glu277\n[Blocked by Tocilizumab]",
           fill = "#FFF1F2", color = "#BE123C", fontface = "bold", size = 2.9, label.size = 1.0) +
  geom_segment(aes(x = 3.6, y = 4.5, xend = 2.85, yend = 4.5),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.2, color = "#DC2626") +
  
  # SITE II (Right): IL-6 <-> gp130 (D2-D3 CHR)
  annotate("label", x = 8.35, y = 4.65,
           label = "SITE II EPITOPE\n(Helix A & Helix C)\nHotspots: Tyr31, Gly35,\nSer118, Val121, Gln124",
           fill = "#DBEAFE", color = "#1E3A8A", fontface = "bold", size = 3.1, label.size = 1.1) +
  annotate("label", x = 8.35, y = 2.75,
           label = "Partner: gp130 (D2-D3 CHR)\nAnchor: Phe191, Val252\n[Blocked by Bazedoxifene]",
           fill = "#EFF6FF", color = "#1D4ED8", fontface = "bold", size = 2.9, label.size = 1.0) +
  geom_segment(aes(x = 6.4, y = 4.5, xend = 7.15, yend = 4.5),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.2, color = "#2563EB") +
  
  # SITE III (Top): IL-6 <-> Second gp130 (Ig-like D1)
  annotate("label", x = 5.0, y = 6.85,
           label = "SITE III EPITOPE (Hexamer Assembly Switch)\nAB-Loop & Helix D Tip: Trp157 (Primary Aromatic Hotspot), Asp160, Thr162\nDocks 2nd gp130 Ig-D1 (Tyr116, Phe118, Leu120) -> Forms Active 2:2:2 Hexamer [Olamkicept]",
           fill = "#FEF3C7", color = "#92400E", fontface = "bold", size = 3.1, label.size = 1.1) +
  geom_segment(aes(x = 5.0, y = 5.6, xend = 5.0, yend = 6.15),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.2, color = "#D97706") +
  
  # Bottom Summary Banner
  annotate("label", x = 5.0, y = 1.35,
           label = "Sequential Assembly:   IL-6 + IL-6Ra (Site I)  -->  + gp130 D2/D3 (Site II Trimer)  -->  x2 via gp130 D1 (Site III Active Hexamer)  -->  JAK1/2 + STAT3 Activation",
           fill = "#F1F5F9", color = "#0F172A", fontface = "bold", size = 3.0, label.size = 0.9) +
  
  xlim(0.2, 9.8) + ylim(0.8, 7.6) +
  theme_void(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = "black", linewidth = 1.0),
    plot.title      = element_text(face = "bold", size = 14.5, hjust = 0.5, color = "black", margin = margin(t = 10, b = 4)),
    plot.subtitle   = element_text(size = 10.5, hjust = 0.5, color = "grey25", margin = margin(b = 6))
  ) +
  labs(
    title    = "Structural Epitopes & Receptor Binding Sites of Human Interleukin-6 (IL-6)",
    subtitle = "Hotspot Amino Acids, Receptor Docking Pockets & Therapeutic Inhibitors (Crystal Structure PDB: 1P9M)"
  )

print(p_il6_map)

# Save High-Res PNG to Desktop -> Senescence_Results
out_dir = file.path(Sys.getenv("USERPROFILE"), "Desktop", "Senescence_Results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

il6_site_path = file.path(out_dir, "23_IL6_Binding_Sites_I_II_III_Structure.png")
ggsave(il6_site_path, plot = p_il6_map, width = 13.5, height = 8.0, dpi = 300, bg = "white")

browseURL(il6_site_path)
# 1. Install lightweight screenshot & jpeg packages (run once)
for (pkg in c("webshot2", "jpeg", "png")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, type = "binary")
}

library(r3dmol)
library(htmlwidgets)
library(webshot2)
library(png)
library(jpeg)

# 2. Build 3D IL-6 / IL-6Ra / gp130 Structure with Binding Sites Highlighted
view_il6_3d = r3dmol(width = 1200, height = 900) %>%
  m_add_model(data = m_fetch_pdb("1P9M"), format = "pdb") %>%
  m_set_style(sel = m_sel(chain = "A"), style = m_style_cartoon(color = "#0D9488")) %>%
  m_set_style(sel = m_sel(chain = "B"), style = m_style_cartoon(color = "#3B82F6", opacity = 0.8)) %>%
  m_set_style(sel = m_sel(chain = "C"), style = m_style_cartoon(color = "#10B981", opacity = 0.8)) %>%
  m_set_style(sel = m_sel(chain = "A", resi = c(74, 78, 178, 179, 182)), style = m_style_sphere(color = "#DC2626", scale = 0.9)) %>%
  m_set_style(sel = m_sel(chain = "A", resi = c(31, 35, 118, 121, 124)), style = m_style_sphere(color = "#2563EB", scale = 0.9)) %>%
  m_set_style(sel = m_sel(chain = "A", resi = c(157, 160, 162)),         style = m_style_sphere(color = "#D97706", scale = 1.0)) %>%
  m_zoom_to(sel = m_sel(chain = c("A", "B", "C")))

# 3. Render to Temporary HTML -> PNG -> High-Quality JPG on Desktop
out_dir  = file.path(Sys.getenv("USERPROFILE"), "Desktop", "Senescence_Results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

tmp_html = tempfile(fileext = ".html")
tmp_png  = tempfile(fileext = ".png")
jpg_path = file.path(out_dir, "24_3D_IL6_Receptor_Complex_1P9M.jpg")

saveWidget(view_il6_3d, file = tmp_html, selfcontained = TRUE)
webshot2::webshot(tmp_html, file = tmp_png, vwidth = 1200, vheight = 900, delay = 2)

# Convert PNG to 100%-Quality JPG
img_rgb = readPNG(tmp_png)
writeJPEG(img_rgb, target = jpg_path, quality = 1.0)

message("Saved 3D JPG structure to: ", jpg_path)
shell.exec(normalizePath(jpg_path))
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)

graphics.off()
set.seed(42)

is_gt = function(a, b) { sign(a - b) == 1 }
is_lt = function(a, b) { sign(b - a) == 1 }

# 1. Simulate Realistic LC-MS/MS Quantitative Proteomics (Log2 LFQ Ratio: Senescent vs. Proliferating)
#    Captures both Intracellular Whole-Cell Lysate and Conditioned-Media Secretome (SASP Atlas)
prot_data = data.frame(
  Protein = c(
    # Secreted SASP Factors (Strongly enriched in Conditioned Media / Secretome)
    "IL6", "IL1B", "IL1A", "CXCL8", "CCL2", "CXCL10", "CCL5", "MMP1", "MMP3", "SERPINE1", "IGFBP3", "IGFBP7", "GDF15", "VEGFA", "TGFB1",
    # Intracellular Tumor Suppressors & DDR (Stabilized at Protein Level)
    "TP53", "CDKN1A", "CDKN2A", "CDKN2B", "RB1", "PTEN", "ATM", "CHEK2",
    # Intracellular mTORC1, cGAS-STING, JAK-STAT & NF-kB Signaling
    "MTOR", "RPTOR", "RPS6KB1", "SQSTM1", "GLB1", "CGAS", "STING1", "TBK1", "IRF3", "JAK1", "JAK2", "STAT3", "STAT1", "NFKB1", "RELA", "SHP2",
    # Downregulated Cell Cycle & Nuclear Lamina Proteins
    "MKI67", "TOP2A", "CDK1", "CDK2", "CDK4", "CCNA2", "CCNB1", "E2F1", "PCNA", "MCM2", "MCM6", "AURKB", "PLK1", "LMNB1", "HMGB1",
    # Background Housekeeping / Structural Proteome (~450 proteins)
    paste0("PROT_", 1:450)
  ),
  Compartment = c(
    rep("Secretome (SASP Media)", 15),
    rep("Intracellular (Tumor Suppressors)", 8),
    rep("Intracellular (JAK-STAT/mTOR/cGAS)", 16),
    rep("Intracellular (Cell Cycle / Lamina)", 15),
    rep("Background Proteome", 450)
  ),
  Log2FC_Protein = c(
    # Secretome SASP Log2FC (Massive extracellular accumulation)
    4.8, 4.2, 3.5, 5.1, 4.4, 3.9, 3.6, 5.4, 5.2, 4.6, 3.8, 4.1, 4.5, 3.3, 2.9,
    # Tumor Suppressors (High protein stabilization, especially TP53 & CDKN1A)
    3.4, 4.1, 3.9, 2.8, 0.4, 1.5, 1.8, 1.6,
    # JAK-STAT / mTORC1 / cGAS-STING / NF-kB total protein abundance
    2.1, 1.9, 1.7, 3.2, 3.8, 2.6, 2.8, 1.9, 1.8, 1.5, 1.6, 2.2, 2.4, 2.3, 2.5, 1.4,
    # Cell Cycle & LMNB1 loss (Classic senescence protein loss)
    -4.6, -4.3, -3.9, -3.1, -2.2, -3.8, -4.1, -3.5, -3.2, -2.9, -2.8, -3.7, -3.6, -3.4, -2.7,
    # Background proteome centered at 0
    rnorm(450, mean = 0, sd = 0.38)
  ),
  FDR = c(
    10^(-runif(15, 8, 22)),
    10^(-runif(8, 5, 16)),
    10^(-runif(16, 4, 14)),
    10^(-runif(15, 7, 20)),
    runif(450, 0.08, 0.99)
  )
)

# Add mRNA Log2FC for Top Key Drivers to illustrate Post-Translational Stabilization (e.g., TP53, RB1, IL6)
prot_data$Log2FC_mRNA = prot_data$Log2FC_Protein * runif(nrow(prot_data), 0.75, 1.05)
# TP53 is regulated primarily via Ser15/Ser20 protein stabilization rather than mRNA transcription
prot_data$Log2FC_mRNA[prot_data$Protein == "TP53"]   = 0.35
prot_data$Log2FC_mRNA[prot_data$Protein == "RB1"]    = -0.15
prot_data$Log2FC_mRNA[prot_data$Protein == "LMNB1"]  = -1.40
prot_data$Log2FC_mRNA[prot_data$Protein == "IL6"]    = 3.10
prot_data$Log2FC_mRNA[prot_data$Protein == "MMP3"]   = 3.40

# Label key proteins on Volcano Plot
key_labels = c("IL6", "IL1B", "CXCL8", "MMP3", "SERPINE1", "GDF15", "TP53", "CDKN1A", "CDKN2A",
               "MTOR", "CGAS", "STING1", "STAT3", "RELA", "GLB1", "SQSTM1", "MKI67", "TOP2A", "CDK1", "CCNB1", "LMNB1")
prot_data$Label = ifelse(prot_data$Protein %in% key_labels, prot_data$Protein, "")

# 2. Top Panel: LC-MS/MS Quantitative Proteomics & Secretome Volcano Plot
p_prot_volcano = ggplot(prot_data, aes(x = Log2FC_Protein, y = -log10(FDR + 1e-300), color = Compartment)) +
  geom_point(alpha = 0.82, size = 2.2) +
  geom_vline(xintercept = c(-1.0, 1.0), linetype = "dashed", color = "grey30", linewidth = 0.65) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30", linewidth = 0.65) +
  geom_text_repel(
    aes(label = Label),
    size          = 3.3,
    fontface      = "bold",
    max.overlaps  = 35,
    box.padding   = 0.45,
    point.padding = 0.25,
    segment.color = "grey45",
    segment.size  = 0.35,
    show.legend   = FALSE
  ) +
  scale_color_manual(values = c(
    "Secretome (SASP Media)"              = "#DC2626",
    "Intracellular (Tumor Suppressors)"   = "#10B981",
    "Intracellular (JAK-STAT/mTOR/cGAS)"  = "#D97706",
    "Intracellular (Cell Cycle / Lamina)" = "#2563EB",
    "Background Proteome"                 = "grey80"
  )) +
  theme_bw(base_size = 12) +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = "black", linewidth = 0.9),
    panel.grid.major  = element_line(color = "grey92"),
    panel.grid.minor  = element_blank(),
    legend.position   = "top",
    legend.title      = element_blank(),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text       = element_text(face = "bold", size = 9.5, color = "black"),
    axis.text         = element_text(face = "bold", color = "black"),
    axis.title        = element_text(face = "bold", color = "black"),
    plot.title        = element_text(face = "bold", hjust = 0.5, size = 14, color = "black"),
    plot.subtitle     = element_text(hjust = 0.5, size = 10, color = "grey30")
  ) +
  labs(
    title    = "LC-MS/MS Quantitative Proteomics & SASP Secretome of Senescent Cancer Cells",
    subtitle = "Whole-Cell Lysate vs. Conditioned-Media Secretome (Senescent / Proliferating Log2 LFQ Ratio)",
    x        = "Log2 Protein Fold Change (Senescent vs. Proliferating)",
    y        = "-Log10 Adjusted P-Value (FDR)"
  )

# 3. Bottom Panel: Transcriptome (scRNA-seq) vs. Proteome (LC-MS/MS) Discordance Plot
comp_df = prot_data[prot_data$Compartment != "Background Proteome", ]

p_mrna_prot = ggplot(comp_df, aes(x = Log2FC_mRNA, y = Log2FC_Protein, color = Compartment)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.7) +
  geom_hline(yintercept = 0, color = "grey75", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.4) +
  geom_point(size = 3.0, alpha = 0.88) +
  geom_text_repel(
    aes(label = Label),
    size          = 3.2,
    fontface      = "bold",
    max.overlaps  = 30,
    box.padding   = 0.4,
    segment.color = "grey50",
    show.legend   = FALSE
  ) +
  annotate("text", x = -1.8, y = 4.2, label = "Post-Translational Stabilization & Secretion\n(e.g., TP53 Ser15/20 Stabilization, IL-6/MMP3 Secretome)",
           hjust = 0, fontface = "bold", size = 3.2, color = "#991B1B") +
  annotate("text", x = -4.5, y = -4.6, label = "Coordinated mRNA + Protein Loss\n(MKI67, CDK1, TOP2A, LMNB1)",
           hjust = 0, fontface = "bold", size = 3.2, color = "#1E40AF") +
  scale_color_manual(values = c(
    "Secretome (SASP Media)"              = "#DC2626",
    "Intracellular (Tumor Suppressors)"   = "#10B981",
    "Intracellular (JAK-STAT/mTOR/cGAS)"  = "#D97706",
    "Intracellular (Cell Cycle / Lamina)" = "#2563EB"
  )) +
  theme_bw(base_size = 11.5) +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = "black", linewidth = 0.9),
    panel.grid.major  = element_line(color = "grey92"),
    panel.grid.minor  = element_blank(),
    legend.position   = "none",
    axis.text         = element_text(face = "bold", color = "black"),
    axis.title        = element_text(face = "bold", color = "black"),
    plot.title        = element_text(face = "bold", hjust = 0.5, size = 13, color = "black")
  ) +
  labs(
    title = "Multi-Omics Integration: scRNA-seq Transcriptome (X) vs. LC-MS/MS Proteome (Y)",
    x     = "scRNA-seq Log2 Fold Change (mRNA Transcripts)",
    y     = "LC-MS/MS Log2 Fold Change (Protein LFQ)"
  )

# 4. Combine & Save High-Res PNG + Full Proteomics Table CSV to Desktop
p_proteomics_full = (p_prot_volcano / p_mrna_prot) + plot_layout(heights = c(1.2, 1.0))
print(p_proteomics_full)

out_dir = file.path(Sys.getenv("USERPROFILE"), "Desktop", "Senescence_Results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

prot_png = file.path(out_dir, "25_LCMS_Proteomics_and_MultiOmics_Integration.png")
ggsave(prot_png, plot = p_proteomics_full, width = 12.5, height = 10, dpi = 300, bg = "white")

write.csv(prot_data, file.path(out_dir, "25_LCMS_Proteomics_Secretome_Table.csv"), row.names = FALSE)
browseURL(prot_png)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(patchwork)

graphics.off()

# 1. Top Panel: Epigenetic Marks Across Senescence & SASP Genomic Loci (Z-Score)
epi_mat = matrix(
  c(
    # CDKN2A (p16INK4a) Promoter: Hypermethylated (5mC) in Cancer -> Demethylated & H3K27ac in Senescence
    1.40, -1.10, -1.20,  1.25,
    # CDKN1A (p21) Promoter: Rapidly hyperacetylated (H3K27ac) by p53/p300 & HDACi
    -1.15,  1.25, -1.10,  1.30,
    # Cell Cycle Genes (MKI67 / CCNA2): Locked into H3K9me3 SAHF Heterochromatin
    -1.20, -1.15,  1.45, -1.25,
    # Endogenous Retroviruses (ERV / LINE-1): Loss of 5mC & H3K9me3 -> Triggers cGAS-STING
    1.35, -0.95, -1.15,  1.10,
    # IL6 / IL1B / CXCL8 SASP Super-Enhancers: TET2 Demethylated + H3K27ac + Ac-STAT3/Ac-p65
    1.10, -1.20, -0.95,  1.45
  ),
  nrow = 5, byrow = TRUE,
  dimnames = list(
    c("CDKN2A (p16INK4a) Locus",
      "CDKN1A (p21) Promoter",
      "Cell Cycle (MKI67 / CCNA2)",
      "ERV / LINE-1 Repeats (cGAS)",
      "SASP Super-Enhancers (IL6/IL1B)"),
    c("DNA Methylation\n(5mC / DNMT1)",
      "Active Enhancer\n(H3K27ac / p300)",
      "Repressive SAHF\n(H3K9me3)",
      "TF Acetylation\n(Ac-p65 / Ac-STAT3)")
  )
)

out_dir = file.path(Sys.getenv("USERPROFILE"), "Desktop", "Senescence_Results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

epi_hm_path = file.path(out_dir, "26A_Epigenetic_Chromatin_Landscape_Heatmap.png")

pheatmap(epi_mat,
         cluster_cols    = FALSE,
         cluster_rows    = FALSE,
         display_numbers = TRUE,
         number_format   = "%.2f",
         fontsize_number = 9.5,
         fontsize_row    = 10.5,
         fontsize_col    = 10.5,
         angle_col       = 0,
         color           = colorRampPalette(c("#2563EB", "white", "#DC2626"))(100),
         main            = "Epigenetic Remodeling Across Senescence & SASP Genomic Loci (Senescent State)",
         filename        = epi_hm_path,
         width           = 10.5,
         height          = 5.2)

# 2. Bottom Panel: Epigenetic Drug Interventions on p16/p21 Arrest, IL-6 SASP & Senolytic Death
epi_drugs = data.frame(
  Treatment = factor(c("Untreated Tumor",
                       "Decitabine (DNMTi)",
                       "JQ1 (BRD4 / BETi)",
                       "Panobinostat (HDACi)",
                       "Panobinostat + Ruxolitinib"),
                     levels = c("Untreated Tumor",
                                "Decitabine (DNMTi)",
                                "JQ1 (BRD4 / BETi)",
                                "Panobinostat (HDACi)",
                                "Panobinostat + Ruxolitinib")),
  Senescence_Arrest = c(10, 85, 88, 92, 94),
  IL6_SASP_Output   = c(25, 78, 18, 95, 14),
  Apoptotic_Death   = c(8,  22, 15, 68, 86)
)

epi_long = tidyr::pivot_longer(epi_drugs,
                               cols      = c("Senescence_Arrest", "IL6_SASP_Output", "Apoptotic_Death"),
                               names_to  = "Readout",
                               values_to = "Score")

epi_long$Readout = recode(epi_long$Readout,
                          "Senescence_Arrest" = "p16 / p21 Senescence Arrest (%)",
                          "IL6_SASP_Output"   = "IL-6 / NF-kB SASP Secretion (%)",
                          "Apoptotic_Death"   = "Senolytic Tumor Cell Death (%)")

p_epi_bars = ggplot(epi_long, aes(x = Treatment, y = Score, fill = Readout)) +
  geom_bar(stat = "identity", position = position_dodge(0.78), width = 0.68, color = "black", linewidth = 0.4) +
  geom_text(aes(label = paste0(Score, "%")), position = position_dodge(0.78), vjust = -0.4, fontface = "bold", size = 3.1) +
  scale_fill_manual(values = c(
    "p16 / p21 Senescence Arrest (%)" = "#2563EB",
    "IL-6 / NF-kB SASP Secretion (%)" = "#F59E0B",
    "Senolytic Tumor Cell Death (%)"  = "#DC2626"
  )) +
  ylim(0, 112) +
  theme_bw(base_size = 11.5) +
  theme(
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = "black", linewidth = 0.85),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "top",
    legend.title       = element_blank(),
    legend.background  = element_rect(fill = "white", color = NA),
    legend.text        = element_text(face = "bold", size = 10, color = "black"),
    axis.text.x        = element_text(angle = 18, hjust = 1, face = "bold", color = "black", size = 10),
    axis.text.y        = element_text(face = "bold", color = "black"),
    axis.title         = element_text(face = "bold", color = "black"),
    plot.title         = element_text(face = "bold", hjust = 0.5, size = 13, color = "black"),
    plot.subtitle      = element_text(hjust = 0.5, size = 10, color = "grey30")
  ) +
  labs(
    title    = "Epigenetic Therapy Effects: DNMTi, BETi (JQ1), HDACi & Rational HDACi + JAKi Combination",
    subtitle = "Combining Panobinostat (HDACi) + Ruxolitinib (JAK1/2i) Blocks Ac-STAT3/Ac-p65 SASP Surge While Maximizing Senolysis",
    x        = "",
    y        = "Relative Response (%)"
  )

print(p_epi_bars)

epi_bar_path = file.path(out_dir, "26B_Epigenetic_Drugs_HDACi_DNMTi_SASP_Outcomes.png")
ggsave(epi_bar_path, plot = p_epi_bars, width = 12, height = 6.5, dpi = 300, bg = "white")

browseURL(epi_bar_path)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

graphics.off()
set.seed(42)

# 1. Simulate High-Resolution Pseudo-Bulk scATAC-seq Coverage Tracks (100-bp Genomic Bins)
make_atac_locus = function(locus_name, chr_label, start_bp, end_bp,
                           prolif_peaks, early_peaks, sasp_peaks,
                           tss_bp, enh_bp, motif_label) {
  bins = seq(start_bp, end_bp, length.out = 250)
  
  calc_signal = function(peaks_list, baseline = 4) {
    sig = rep(baseline, length(bins)) + runif(length(bins), 0, 2.5)
    for (pk in peaks_list) {
      # pk = c(center_bp, width_bp, height)
      sig = sig + pk[3] * exp(-0.5 * ((bins - pk[1]) / pk[2])^2)
    }
    return(sig)
  }
  
  df = rbind(
    data.frame(Locus = locus_name, Chr = chr_label, Pos_kb = bins / 1000,
               CellState = "Proliferating",   Accessibility = calc_signal(prolif_peaks)),
    data.frame(Locus = locus_name, Chr = chr_label, Pos_kb = bins / 1000,
               CellState = "Early_Senescent", Accessibility = calc_signal(early_peaks)),
    data.frame(Locus = locus_name, Chr = chr_label, Pos_kb = bins / 1000,
               CellState = "SASP_Senescent",  Accessibility = calc_signal(sasp_peaks))
  )
  df$CellState = factor(df$CellState, levels = c("Proliferating", "Early_Senescent", "SASP_Senescent"))
  
  ann = data.frame(
    Locus = locus_name,
    TSS_kb = tss_bp / 1000,
    Enh_kb = enh_bp / 1000,
    Motif  = motif_label
  )
  list(tracks = df, ann = ann)
}

# Build 3 Key Loci: CDKN2A (p16), IL6 (SASP Super-Enhancer), and MKI67 (SAHF Compaction)
loc_cdkn2a = make_atac_locus(
  locus_name   = "CDKN2A (p16INK4a) Locus [chr9p21.3]",
  chr_label    = "chr9 (kb)",
  start_bp     = 21965000, end_bp = 21998000,
  prolif_peaks = list(c(21975000, 600, 6)),
  early_peaks  = list(c(21975000, 650, 78), c(21988000, 800, 52)),
  sasp_peaks   = list(c(21975000, 650, 84), c(21988000, 800, 60)),
  tss_bp       = 21975000, enh_bp = 21988000,
  motif_label  = "Peaks: Promoter TSS + Enhancer (ETS1 / TP53)"
)

loc_il6 = make_atac_locus(
  locus_name   = "IL6 & SASP Super-Enhancer [chr7p15.3]",
  chr_label    = "chr7 (kb)",
  start_bp     = 22718000, end_bp = 22736000,
  prolif_peaks = list(c(22727200, 500, 5)),
  early_peaks  = list(c(22727200, 550, 18), c(22721500, 600, 15)),
  sasp_peaks   = list(c(22727200, 550, 92), c(22721500, 700, 88), c(22724200, 500, 64)),
  tss_bp       = 22727200, enh_bp = 22721500,
  motif_label  = "Super-Enhancer Footprints: STAT3 + NF-kB (RELA) + BRD4"
)

loc_mki67 = make_atac_locus(
  locus_name   = "MKI67 Proliferation Locus [chr10q26.2]",
  chr_label    = "chr10 (kb)",
  start_bp     = 128095000, end_bp = 128128000,
  prolif_peaks = list(c(128125000, 700, 95), c(128112000, 850, 68)),
  early_peaks  = list(c(128125000, 600, 10)),
  sasp_peaks   = list(c(128125000, 600, 5)),
  tss_bp       = 128125000, enh_bp = 128112000,
  motif_label  = "SAHF Compaction: Loss of E2F1 / MYC Accessibility"
)

# Helper to plot normalized scATAC-seq pileup tracks per locus
plot_atac_panel = function(loc_obj, fill_cols) {
  tr  = loc_obj$tracks
  ann = loc_obj$ann
  
  ggplot(tr, aes(x = Pos_kb, y = Accessibility, fill = CellState, color = CellState)) +
    geom_area(alpha = 0.78, linewidth = 0.6) +
    geom_vline(xintercept = c(ann$TSS_kb, ann$Enh_kb), linetype = "dashed", color = "black", linewidth = 0.55) +
    facet_wrap(~ CellState, ncol = 1, strip.position = "left") +
    scale_fill_manual(values = fill_cols) +
    scale_color_manual(values = fill_cols) +
    ylim(0, 112) +
    theme_bw(base_size = 10.5) +
    theme(
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = "black", linewidth = 0.75),
      panel.grid        = element_blank(),
      strip.background  = element_rect(fill = "#0F172A"),
      strip.text.y.left = element_text(angle = 0, color = "white", face = "bold", size = 8.5),
      legend.position   = "none",
      axis.text.y       = element_blank(),
      axis.ticks.y      = element_blank(),
      axis.text.x       = element_text(face = "bold", color = "black", size = 8.5),
      axis.title        = element_text(face = "bold", color = "black", size = 9.5),
      plot.title        = element_text(face = "bold", hjust = 0.5, size = 11, color = "black"),
      plot.subtitle     = element_text(hjust = 0.5, size = 8.5, color = "grey25")
    ) +
    labs(title = ann$Locus, subtitle = ann$Motif, x = tr$Chr[1], y = "Normalized ATAC Pileup")
}

state_cols = c("Proliferating" = "#2563EB", "Early_Senescent" = "#F59E0B", "SASP_Senescent" = "#DC2626")

p_track_1 = plot_atac_panel(loc_cdkn2a, state_cols)
p_track_2 = plot_atac_panel(loc_il6,    state_cols)
p_track_3 = plot_atac_panel(loc_mki67,  state_cols)

# 2. Bottom Panel: Transcription Factor Motif Accessibility Deviation (chromVAR Z-Score)
motif_df = data.frame(
  Motif = factor(c("E2F1 (Cell Cycle)", "MYC (Proliferation)", "TP53 (DDR Arrest)",
                   "ETS1 / CEBPB", "STAT3 (TTCCNGGAA)", "NF-kB RELA (GGGRNNYYCC)", "IRF3 / ISRE (cGAS)"),
                 levels = rev(c("E2F1 (Cell Cycle)", "MYC (Proliferation)", "TP53 (DDR Arrest)",
                                "ETS1 / CEBPB", "STAT3 (TTCCNGGAA)", "NF-kB RELA (GGGRNNYYCC)", "IRF3 / ISRE (cGAS)"))),
  Proliferating   = c( 2.4,  2.1, -1.6, -1.2, -1.5, -1.7, -1.3),
  Early_Senescent = c(-1.8, -1.5,  2.3,  1.1,  0.4,  0.6,  0.2),
  SASP_Senescent  = c(-2.3, -2.1,  2.0,  2.2,  2.7,  2.9,  2.5)
)

motif_long = tidyr::pivot_longer(motif_df, cols = c("Proliferating", "Early_Senescent", "SASP_Senescent"),
                                 names_to = "CellState", values_to = "chromVAR_Z")
motif_long$CellState = factor(motif_long$CellState, levels = c("Proliferating", "Early_Senescent", "SASP_Senescent"))

p_motif_bar = ggplot(motif_long, aes(x = chromVAR_Z, y = Motif, fill = CellState)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.7) +
  geom_bar(stat = "identity", position = position_dodge(0.75), width = 0.68, color = "black", linewidth = 0.35) +
  scale_fill_manual(values = state_cols) +
  xlim(-2.8, 3.3) +
  theme_bw(base_size = 11) +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = "black", linewidth = 0.85),
    panel.grid.major.y = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "top",
    legend.title      = element_blank(),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text       = element_text(face = "bold", size = 10, color = "black"),
    axis.text         = element_text(face = "bold", color = "black", size = 10),
    axis.title        = element_text(face = "bold", color = "black", size = 10.5),
    plot.title        = element_text(face = "bold", hjust = 0.5, size = 13, color = "black")
  ) +
  labs(title = "scATAC-seq Transcription Factor Motif Accessibility Deviation (chromVAR Z-Score)",
       x     = "Motif Accessibility Z-Score (Closed < 0 | Open > 0)",
       y     = "")

# 3. Combine ATAC Tracks (Top Row) + Motif Accessibility Chart (Bottom Row)
p_scatac_full = (p_track_1 | p_track_2 | p_track_3) / p_motif_bar +
  plot_layout(heights = c(1.35, 1.0)) +
  plot_annotation(
    title    = "Single-Cell Chromatin Accessibility (scATAC-seq) & TF Motif Footprinting in Senescence",
    subtitle = "Dashed Vertical Lines Indicate Promoter TSS & Distal Super-Enhancer Peak Coordinates",
    theme    = theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.title      = element_text(size = 14.5, face = "bold", hjust = 0.5, color = "black"),
      plot.subtitle   = element_text(size = 10.5, hjust = 0.5, color = "grey30")
    )
  )

print(p_scatac_full)

# 4. Save High-Res PNG & chromVAR Motif Table to Desktop -> Senescence_Results
out_dir = file.path(Sys.getenv("USERPROFILE"), "Desktop", "Senescence_Results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

atac_png = file.path(out_dir, "27_scATAC_seq_Chromatin_Accessibility_and_Motifs.png")
ggsave(atac_png, plot = p_scatac_full, width = 14.5, height = 9.5, dpi = 300, bg = "white")

write.csv(motif_df, file.path(out_dir, "27_scATAC_chromVAR_Motif_Scores.csv"), row.names = FALSE)
browseURL(atac_png)
# Install required packages if missing
# install.packages(c("reticulate", "ggplot2", "dplyr", "pheatmap"))
# BiocManager::install(c("GSVA", "msigdbr"))

library(reticulate)
library(GSVA)
library(msigdbr)
library(ggplot2)
library(pheatmap)

# 1. Define Curated Gene Signatures for IL-6 Anterograde Trafficking in Senescence
il6_trafficking_signatures <- list(
  # ER-to-Golgi Anterograde Transport (COPII & ERGIC)
  ER_to_Golgi_Anterograde = c(
    "SEC23A", "SEC23B", "SEC24A", "SEC24B", "SEC24C", "SEC24D",
    "SEC13", "SEC31A", "SAR1A", "SAR1B", "LMAN1", "MCFD2", "RAB1A", "RAB1B", "USO1"
  ),
  # Golgi-to-Plasma Membrane & Secretory Vesicle Exocytosis
  Post_Golgi_Exocytosis = c(
    "ARF1", "ARF4", "GBF1", "RAB6A", "RAB8A", "RAB11A",
    "VAMP3", "VAMP7", "SNAP23", "STX4", "EXOC2", "EXOC4"
  ),
  # Intracrine IL-6 / Senescence / cGAS-STING-NFkB Axis
  Senescence_Intracrine_IL6 = c(
    "IL6", "IL6R", "IL6ST", "CGAS", "STING1", "TBK1",
    "RELA", "NFKB1", "CDKN1A", "CDKN2A", "STAT3", "CEBPB"
  )
)

# 2. Example: Simulate Normalized Log2 Expression Matrix (Replace with your count/TPM matrix)
set.seed(42)
all_genes <- unique(unlist(il6_trafficking_signatures))
samples <- c(paste0("Control_", 1:4), paste0("Senescent_Cancer_", 1:4))
expr_matrix <- matrix(rnorm(length(all_genes) * length(samples), mean = 6, sd = 1.2),
                      nrow = length(all_genes), dimnames = list(all_genes, samples))

# Simulate upregulation of anterograde trafficking & IL-6 intracrine genes in Senescent samples
expr_matrix[all_genes, 5:8] <- expr_matrix[all_genes, 5:8] + runif(length(all_genes), 1.2, 2.8)

# 3. Run ssGSEA Pathway Activity Scoring
gsva_param <- ssgseaParam(exprData = expr_matrix, geneSets = il6_trafficking_signatures)
pathway_scores <- gsva(gsva_param)

# 4. Plot Heatmap of Pathway Activity
pheatmap(pathway_scores,
         main = "IL-6 Anterograde Trafficking & Senescence Pathway Activity",
         cluster_cols = FALSE,
         color = colorRampPalette(c("#2b5c8f", "#f7f7f7", "#b30000"))(50))
install.packages(c("ggplot2", "pheatmap", "BiocManager"))
BiocManager::install("GSVA")
# Run this block inside RStudio to execute the embedded Python script
py_run_string("
import requests
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# 1. Core proteins linking IL-6, anterograde vesicular transport, and senescence
genes = [
    'IL6', 'IL6R', 'IL6ST',          # Cytokine & Receptor complex
    'SEC23A', 'SEC24D', 'SAR1A',     # COPII ER-exit machinery
    'LMAN1', 'RAB1A', 'ARF1',        # ERGIC & Golgi sorting
    'RAB8A', 'VAMP3', 'SNAP23',      # Anterograde secretory vesicles / SNAREs
    'STING1', 'RELA', 'CDKN1A'       # Senescence / SASP regulatory axis
]

# 2. Query STRING-db API for functional & physical protein interactions
string_api_url = 'https://string-db.org'
params = {
    'identifiers': '%0d'.join(genes),
    'species': 9606,          # Homo sapiens
    'required_score': 400     # Medium-high confidence (0.4)
}

response = requests.post(string_api_url, data=params)
interactions = response.json()

# 3. Build Directed/Weighted Network Graph
G = nx.Graph()
for edge in interactions:
    p1 = edge['preferredName_A']
    p2 = edge['preferredName_B']
    score = edge['score']
    G.add_edge(p1, p2, weight=score)

# 4. Calculate Bottleneck / Hub Centrality in Anterograde Trafficking
betweenness = nx.betweenness_centrality(G, weight='weight')
centrality_df = pd.DataFrame(
    list(betweenness.items()), columns=['Protein', 'Betweenness_Centrality']
).sort_values(by='Betweenness_Centrality', ascending=False)

print('Top Hub Regulators in IL-6 Anterograde Trafficking Network:')
print(centrality_df.head(8).to_string(index=False))

# 5. Visualize the Trafficking Network
plt.figure(figsize=(10, 7))
pos = nx.spring_layout(G, seed=42, k=0.6)
nx.draw_networkx_nodes(G, pos, node_size=900, node_color='#4C72B0', alpha=0.9)
nx.draw_networkx_edges(G, pos, width=1.5, alpha=0.5, edge_color='gray')
nx.draw_networkx_labels(G, pos, font_size=9, font_color='white', font_weight='bold')
plt.title('IL-6 Intracellular Anterograde Trafficking & Senescence Network (STRING)', fontsize=12)
plt.axis('off')
plt.tight_layout()
plt.savefig('IL6_anterograde_network.png', dpi=300)
plt.show()
")

# Bring the Python centrality dataframe back into R for downstream analysis
r_centrality_table <- py$centrality_df
print(r_centrality_table)

# ==============================================================================
# RStudio + Python Pipeline:
# Regnase-1 (ZC3H12A) mRNA Decay Brake vs. IL-6 Anterograde Vesicular Trafficking
# ==============================================================================

# 1. Auto-install missing R packages
required_r_pkgs <- c("reticulate", "ggplot2", "pheatmap", "ggpubr", "BiocManager")
missing_r_pkgs <- required_r_pkgs[!(required_r_pkgs %in% installed.packages()[, "Package"])]
if (length(missing_r_pkgs) > 0) install.packages(missing_r_pkgs)

if (!("GSVA" %in% installed.packages()[, "Package"])) {
  BiocManager::install("GSVA", ask = FALSE, update = FALSE)
}

library(reticulate)
library(GSVA)
library(ggplot2)
library(ggpubr)
library(pheatmap)

# 2. Connect to your existing Python installation
# use_python("C:/Path/To/Your/python.exe", required = TRUE) # Uncomment if needed
required_py_pkgs <- c("requests", "pandas", "networkx", "matplotlib")
for (pkg in required_py_pkgs) {
  if (!py_module_available(pkg)) py_install(pkg)
}

# ==============================================================================
# PART A (R): Pathway Scoring & Negative Correlation Modeling
# ==============================================================================

pathway_signatures <- list(
  # Post-transcriptional mRNA Decay Brake (Active in proliferating, lost in senescence)
  Regnase1_mRNA_Decay_Brake = c(
    "ZC3H12A", "UPF1", "RC3H1", "RC3H2", "SMG1", "SMG5", "SMG7", "ZFP36"
  ),
  # ER-to-Golgi Anterograde Transport (COPII & ERGIC)
  ER_to_Golgi_COPII = c(
    "SEC23A", "SEC23B", "SEC24A", "SEC24B", "SEC24C", "SEC24D",
    "SEC13", "SEC31A", "SAR1A", "SAR1B", "LMAN1", "MCFD2", "RAB1A", "RAB1B", "USO1"
  ),
  # Post-Golgi Anterograde Secretory & Intracrine Vesicles
  Golgi_to_PM_Anterograde = c(
    "ARF1", "ARF4", "GBF1", "RAB6A", "RAB8A", "RAB11A",
    "VAMP3", "VAMP7", "SNAP23", "STX4", "EXOC2", "EXOC4"
  ),
  # Senescence & Intracrine IL-6 / IKK-NFkB Axis
  Senescence_Intracrine_IL6 = c(
    "IL6", "IL6R", "IL6ST", "IKBKB", "IKBKG", "CGAS", "STING1",
    "RELA", "NFKB1", "CDKN1A", "CDKN2A", "STAT3"
  )
)

# Simulated log2 expression matrix (Replace with your actual expression matrix)
set.seed(42)
all_genes <- unique(unlist(pathway_signatures))
samples <- c(paste0("Proliferating_", 1:6), paste0("Senescent_Cancer_", 1:6))

expr_matrix <- matrix(
  rnorm(length(all_genes) * length(samples), mean = 6, sd = 0.8),
  nrow = length(all_genes),
  dimnames = list(all_genes, samples)
)

# Model biological switch in Senescent_Cancer samples:
# 1. Downregulate Regnase-1 / UPF1 / Roquin-1 mRNA decay complex
decay_genes <- pathway_signatures$Regnase1_mRNA_Decay_Brake
expr_matrix[decay_genes, 7:12] <- expr_matrix[decay_genes, 7:12] - runif(length(decay_genes), 1.8, 3.0)

# 2. Upregulate IL-6, IKK, and Anterograde Vesicular Trafficking machinery
trafficking_genes <- setdiff(all_genes, decay_genes)
expr_matrix[trafficking_genes, 7:12] <- expr_matrix[trafficking_genes, 7:12] + runif(length(trafficking_genes), 1.6, 2.8)

# Run ssGSEA Pathway Activity Scoring
gsva_param <- ssgseaParam(exprData = expr_matrix, geneSets = pathway_signatures)
pathway_scores <- gsva(gsva_param)

# Plot 1: Heatmap showing inverse activity between Regnase-1 Brake and Anterograde Trafficking
pheatmap(
  pathway_scores,
  main = "Regnase-1 mRNA Decay vs. IL-6 Anterograde Trafficking in Senescence",
  cluster_cols = FALSE,
  fontsize = 10,
  color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(50)
)

# Plot 2: Regression Scatter Plot of Negative Correlation
cor_df <- data.frame(
  Sample = colnames(pathway_scores),
  Condition = rep(c("Proliferating", "Senescent Cancer"), each = 6),
  Regnase1_Brake_Score = pathway_scores["Regnase1_mRNA_Decay_Brake", ],
  ER_Golgi_Trafficking_Score = pathway_scores["ER_to_Golgi_COPII", ],
  IL6_Expression = expr_matrix["IL6", ]
)

p_corr <- ggplot(cor_df, aes(x = Regnase1_Brake_Score, y = ER_Golgi_Trafficking_Score)) +
  geom_point(aes(color = Condition, size = IL6_Expression), alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, color = "#333333", linetype = "dashed") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 4.5) +
  scale_color_manual(values = c("Proliferating" = "#1f78b4", "Senescent Cancer" = "#e31a1c")) +
  labs(
    title = "Inverse Coupling: Regnase-1 Decay vs. IL-6 Anterograde Transport",
    subtitle = "Loss of ZC3H12A/UPF1/RC3H1 stabilizes IL6 mRNA and drives ER-to-Golgi flux",
    x = "Regnase-1 / UPF1 / Roquin-1 Activity (ssGSEA Score)",
    y = "ER-to-Golgi COPII Anterograde Trafficking (ssGSEA Score)",
    size = "IL6 mRNA (log2)"
  ) +
  theme_minimal(base_size = 12)

print(p_corr)

# ==============================================================================
# PART B (Python via Reticulate): Integrated Regulatory & Trafficking Network
# ==============================================================================

py_run_string("
import requests
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

genes = [
    # Post-transcriptional decay module
    'ZC3H12A', 'UPF1', 'RC3H1', 'RC3H2',
    # IKK / Senescence switch
    'IKBKB', 'RELA', 'STING1', 'CDKN1A',
    # Cytokine & Receptor
    'IL6', 'IL6R', 'IL6ST',
    # Anterograde COPII / Golgi / Secretory vesicles
    'SEC23A', 'SEC24D', 'SAR1A', 'LMAN1', 'RAB1A', 'ARF1', 'RAB8A', 'VAMP3'
]

url = 'https://string-db.org'
params = {
    'identifiers': '%0d'.join(genes),
    'species': 9606,
    'required_score': 400
}

interactions = requests.post(url, data=params).json()

G = nx.Graph()
for edge in interactions:
    G.add_edge(edge['preferredName_A'], edge['preferredName_B'], weight=edge['score'])

# Color-code nodes by functional module
decay_module = {'ZC3H12A', 'UPF1', 'RC3H1', 'RC3H2'}
senescence_module = {'IKBKB', 'RELA', 'STING1', 'CDKN1A', 'IL6', 'IL6R', 'IL6ST'}

node_colors = []
for node in G.nodes():
    if node in decay_module:
        node_colors.append('#2b83ba')   # Blue: Regnase-1 mRNA Decay Brake
    elif node in senescence_module:
        node_colors.append('#d7191c')   # Red: IKK-Senescence / IL-6 Axis
    else:
        node_colors.append('#fdae61')   # Orange: Anterograde Vesicular Trafficking

centrality = nx.betweenness_centrality(G, weight='weight')
hub_df = pd.DataFrame(
    list(centrality.items()), columns=['Gene', 'Betweenness_Centrality']
).sort_values(by='Betweenness_Centrality', ascending=False)

print('\\n--- Network Bridge Regulators (Decay <-> Senescence <-> Anterograde Transport) ---')
print(hub_df.head(10).to_string(index=False))

plt.figure(figsize=(11, 8))
pos = nx.spring_layout(G, seed=42, k=0.7)
nx.draw_networkx_nodes(G, pos, node_size=1000, node_color=node_colors, alpha=0.92)
nx.draw_networkx_edges(G, pos, width=1.6, alpha=0.45, edge_color='#555555')
nx.draw_networkx_labels(G, pos, font_size=8.5, font_color='black', font_weight='bold')
plt.title('Regnase-1 (Blue) vs. IKK/IL-6 (Red) & Anterograde Trafficking (Orange)', fontsize=12)
plt.axis('off')
plt.tight_layout()
plt.savefig('Regnase1_IL6_Anterograde_Network.png', dpi=300)
plt.close()
")

bridge_regulators <- py$hub_df
View(bridge_regulators)

# ==============================================================================
# Pure-R RStudio Pipeline (No Python Required):
# Regnase-1 (ZC3H12A) mRNA Decay Brake vs. IL-6 Anterograde Trafficking
# ==============================================================================

required_pkgs <- c("httr", "jsonlite", "igraph", "ggplot2", "ggpubr", "pheatmap", "BiocManager")
missing_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

if (!("GSVA" %in% installed.packages()[, "Package"])) {
  BiocManager::install("GSVA", ask = FALSE, update = FALSE)
}

library(httr)
library(jsonlite)
library(igraph)
library(GSVA)
library(ggplot2)
library(ggpubr)
library(pheatmap)

# 1. Pathway Signatures
pathway_signatures <- list(
  Regnase1_mRNA_Decay_Brake = c("ZC3H12A", "UPF1", "RC3H1", "RC3H2", "SMG1", "SMG5", "SMG7", "ZFP36"),
  ER_to_Golgi_COPII = c("SEC23A", "SEC23B", "SEC24A", "SEC24B", "SEC24C", "SEC24D", "SEC13", "SEC31A", "SAR1A", "SAR1B", "LMAN1", "MCFD2", "RAB1A", "RAB1B", "USO1"),
  Golgi_to_PM_Anterograde = c("ARF1", "ARF4", "GBF1", "RAB6A", "RAB8A", "RAB11A", "VAMP3", "VAMP7", "SNAP23", "STX4", "EXOC2", "EXOC4"),
  Senescence_Intracrine_IL6 = c("IL6", "IL6R", "IL6ST", "IKBKB", "IKBKG", "CGAS", "STING1", "RELA", "NFKB1", "CDKN1A", "CDKN2A", "STAT3")
)

# 2. Simulated Expression Matrix (Replace with your own matrix)
set.seed(42)
all_genes <- unique(unlist(pathway_signatures))
samples <- c(paste0("Proliferating_", 1:6), paste0("Senescent_Cancer_", 1:6))

expr_matrix <- matrix(rnorm(length(all_genes) * length(samples), mean = 6, sd = 0.8),
                      nrow = length(all_genes), dimnames = list(all_genes, samples))

decay_genes <- pathway_signatures$Regnase1_mRNA_Decay_Brake
expr_matrix[decay_genes, 7:12] <- expr_matrix[decay_genes, 7:12] - runif(length(decay_genes), 1.8, 3.0)
trafficking_genes <- setdiff(all_genes, decay_genes)
expr_matrix[trafficking_genes, 7:12] <- expr_matrix[trafficking_genes, 7:12] + runif(length(trafficking_genes), 1.6, 2.8)

# 3. ssGSEA Pathway Scoring & Correlation Plot
gsva_param <- ssgseaParam(exprData = expr_matrix, geneSets = pathway_signatures)
pathway_scores <- gsva(gsva_param)

cor_df <- data.frame(
  Sample = colnames(pathway_scores),
  Condition = rep(c("Proliferating", "Senescent Cancer"), each = 6),
  Regnase1_Brake_Score = pathway_scores["Regnase1_mRNA_Decay_Brake", ],
  ER_Golgi_Trafficking_Score = pathway_scores["ER_to_Golgi_COPII", ],
  IL6_Expression = expr_matrix["IL6", ]
)

p_corr <- ggplot(cor_df, aes(x = Regnase1_Brake_Score, y = ER_Golgi_Trafficking_Score)) +
  geom_point(aes(color = Condition, size = IL6_Expression), alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, color = "#333333", linetype = "dashed") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 4.5) +
  scale_color_manual(values = c("Proliferating" = "#1f78b4", "Senescent Cancer" = "#e31a1c")) +
  labs(
    title = "Inverse Coupling: Regnase-1 Decay vs. IL-6 Anterograde Transport",
    x = "Regnase-1 / UPF1 / Roquin-1 Activity (ssGSEA Score)",
    y = "ER-to-Golgi COPII Anterograde Trafficking (ssGSEA Score)",
    size = "IL6 mRNA (log2)"
  ) +
  theme_minimal(base_size = 12)

print(p_corr)

# 4. Query STRING-db API & Plot Network Natively in R
genes <- c("ZC3H12A", "UPF1", "RC3H1", "RC3H2", "IKBKB", "RELA", "STING1", "CDKN1A",
           "IL6", "IL6R", "IL6ST", "SEC23A", "SEC24D", "SAR1A", "LMAN1", "RAB1A", "ARF1", "RAB8A", "VAMP3")

res <- POST(
  url = "https://string-db.org",
  body = list(identifiers = paste(genes, collapse = "%0d"), species = 9606, required_score = 400),
  encode = "form"
)

interactions <- fromJSON(content(res, "text", encoding = "UTF-8"))
edges <- interactions[, c("preferredName_A", "preferredName_B", "score")]
g <- graph_from_data_frame(edges, directed = FALSE)
g <- simplify(g)

# Color nodes by functional module
V(g)$color <- ifelse(V(g)$name %in% c("ZC3H12A", "UPF1", "RC3H1", "RC3H2"), "#2b83ba",
                     ifelse(V(g)$name %in% c("IKBKB", "RELA", "STING1", "CDKN1A", "IL6", "IL6R", "IL6ST"), "#d7191c", "#fdae61"))

set.seed(42)
plot(g,
     vertex.size = 24,
     vertex.label.cex = 0.75,
     vertex.label.color = "black",
     edge.width = E(g)$score * 3,
     main = "Regnase-1 (Blue) vs. IKK/IL-6 (Red) & Anterograde Trafficking (Orange)")

# Table of top bridge regulators
bridge_regulators <- data.frame(
  Gene = V(g)$name,
  Betweenness_Centrality = betweenness(g, normalized = TRUE)
)
bridge_regulators <- bridge_regulators[order(-bridge_regulators$Betweenness_Centrality), ]
View(bridge_regulators)
# ==============================================================================
# 100% OFFLINE PURE-R SCRIPT (No Python, No JSON, No External APIs)
# Regnase-1 (ZC3H12A) mRNA Decay vs. IL-6 Anterograde Trafficking in Senescence
# ==============================================================================

# 1. Install only basic CRAN packages if missing
pkgs <- c("igraph", "ggplot2", "pheatmap")
new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

library(igraph)
library(ggplot2)
library(pheatmap)

# 2. Define Functional Gene Signatures
signatures <- list(
  Regnase1_Decay_Brake = c("ZC3H12A", "UPF1", "RC3H1", "RC3H2", "SMG1", "SMG5", "SMG7"),
  ER_Golgi_COPII       = c("SEC23A", "SEC23B", "SEC24A", "SEC24D", "SEC13", "SEC31A", "SAR1A", "LMAN1", "RAB1A"),
  PostGolgi_Exocytosis = c("ARF1", "ARF4", "RAB6A", "RAB8A", "RAB11A", "VAMP3", "SNAP23", "STX4"),
  Senescence_IL6_Axis  = c("IL6", "IL6R", "IL6ST", "IKBKB", "CGAS", "STING1", "RELA", "CDKN1A", "CDKN2A")
)

# 3. Simulate Expression Matrix (Proliferating vs. Senescent Cancer)
set.seed(42)
all_genes <- unique(unlist(signatures))
samples <- c(paste0("Proliferating_", 1:6), paste0("Senescent_", 1:6))

expr <- matrix(
  rnorm(length(all_genes) * length(samples), mean = 6, sd = 0.7),
  nrow = length(all_genes),
  dimnames = list(all_genes, samples)
)

# In Senescence: Regnase-1 decay complex drops; IL-6 & Anterograde trafficking rise
decay_genes <- signatures$Regnase1_Decay_Brake
other_genes <- setdiff(all_genes, decay_genes)

expr[decay_genes, 7:12] <- expr[decay_genes, 7:12] - runif(length(decay_genes), 2.0, 3.2)
expr[other_genes, 7:12] <- expr[other_genes, 7:12] + runif(length(other_genes), 1.8, 3.0)

# 4. Calculate Pathway Activity Scores (Mean Z-score per pathway - no Bioconductor needed)
z_expr <- t(scale(t(expr)))
pathway_scores <- sapply(signatures, function(g) colMeans(z_expr[g, , drop = FALSE]))
pathway_scores <- t(pathway_scores)

# PLOT 1: Heatmap
pheatmap(
  pathway_scores,
  main = "Regnase-1 Decay Brake vs. IL-6 Anterograde Trafficking",
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(50)
)

# 5. Negative Correlation Plot (Regnase-1 Decay vs. ER-to-Golgi COPII Transport)
df <- data.frame(
  Condition = rep(c("Proliferating", "Senescent"), each = 6),
  Regnase1_Score = pathway_scores["Regnase1_Decay_Brake", ],
  COPII_Trafficking_Score = pathway_scores["ER_Golgi_COPII", ],
  IL6_mRNA = expr["IL6", ]
)

cor_test <- cor.test(df$Regnase1_Score, df$COPII_Trafficking_Score)
r_val <- round(cor_test$estimate, 2)
p_val <- signif(cor_test$p.value, 2)

p_scatter <- ggplot(df, aes(x = Regnase1_Score, y = COPII_Trafficking_Score)) +
  geom_point(aes(color = Condition, size = IL6_mRNA), alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  annotate("text", x = min(df$Regnase1_Score), y = max(df$COPII_Trafficking_Score),
           label = paste0("Pearson R = ", r_val, ", p = ", p_val),
           hjust = 0, vjust = 1, size = 4.5, fontface = "bold") +
  scale_color_manual(values = c("Proliferating" = "#1f78b4", "Senescent" = "#e31a1c")) +
  labs(
    title = "Negative Correlation: Regnase-1 Decay vs. IL-6 Anterograde Transport",
    x = "Regnase-1 / UPF1 / Roquin-1 Activity (Z-score)",
    y = "ER-to-Golgi COPII Trafficking Activity (Z-score)",
    size = "IL6 mRNA (log2)"
  ) +
  theme_minimal(base_size = 12)

print(p_scatter)

# 6. Curated Mechanistic Interaction Network (Offline - Zero API calls)
edges <- data.frame(
  from = c("ZC3H12A", "ZC3H12A", "ZC3H12A", "IKBKB",   "STING1", "CGAS",   "RELA",
           "IL6",     "IL6",     "IL6R",    "SEC24D",  "SEC23A", "LMAN1",  "RAB1A",
           "ARF1",    "RAB8A",   "VAMP3",   "IL6",     "RELA"),
  to   = c("UPF1",    "RC3H1",   "IL6",     "ZC3H12A", "IKBKB",  "STING1", "IL6",
           "IL6R",    "SEC24D",  "IL6ST",   "SEC23A",  "SAR1A",  "SEC24D", "LMAN1",
           "RAB1A",   "ARF1",    "RAB8A",   "STAT3",   "CDKN1A")
)

g <- graph_from_data_frame(edges, directed = FALSE)

V(g)$color <- ifelse(V(g)$name %in% c("ZC3H12A", "UPF1", "RC3H1"), "#2b83ba",
                     ifelse(V(g)$name %in% c("IKBKB", "CGAS", "STING1", "RELA", "IL6", "IL6R", "IL6ST", "STAT3", "CDKN1A"), "#d7191c", "#fdae61"))

# PLOT 3: Network Diagram
set.seed(42)
plot(
  g,
  vertex.size = 24,
  vertex.label.cex = 0.75,
  vertex.label.color = "black",
  edge.width = 2,
  edge.color = "#666666",
  main = "Regnase-1 (Blue) -> IKK/IL-6 (Red) -> Anterograde Trafficking (Orange)"
)

# View Hub Ranking Table
hubs <- data.frame(
  Gene = V(g)$name,
  Betweenness_Centrality = round(betweenness(g, normalized = TRUE), 3)
)
hubs <- hubs[order(-hubs$Betweenness_Centrality), ]
print(hubs)
# ==============================================================================
# Clean Pastel Plots & Auto-Export (Heatmap, Scatter Plot, Unmixed Network)
# ==============================================================================

library(igraph)
library(ggplot2)
library(pheatmap)

# 1. Define Functional Gene Signatures
signatures <- list(
  Regnase1_Decay_Brake = c("ZC3H12A", "UPF1", "RC3H1", "RC3H2", "SMG1", "SMG5", "SMG7"),
  ER_Golgi_COPII       = c("SEC23A", "SEC23B", "SEC24A", "SEC24D", "SEC13", "SEC31A", "SAR1A", "LMAN1", "RAB1A"),
  PostGolgi_Exocytosis = c("ARF1", "ARF4", "RAB6A", "RAB8A", "RAB11A", "VAMP3", "SNAP23", "STX4"),
  Senescence_IL6_Axis  = c("IL6", "IL6R", "IL6ST", "IKBKB", "CGAS", "STING1", "RELA", "CDKN1A", "CDKN2A")
)

# 2. Simulate Expression Matrix (Proliferating vs. Senescent Cancer)
set.seed(42)
all_genes <- unique(unlist(signatures))
samples <- c(paste0("Proliferating_", 1:6), paste0("Senescent_", 1:6))

expr <- matrix(
  rnorm(length(all_genes) * length(samples), mean = 6, sd = 0.7),
  nrow = length(all_genes),
  dimnames = list(all_genes, samples)
)

decay_genes <- signatures$Regnase1_Decay_Brake
other_genes <- setdiff(all_genes, decay_genes)
expr[decay_genes, 7:12] <- expr[decay_genes, 7:12] - runif(length(decay_genes), 2.0, 3.2)
expr[other_genes, 7:12] <- expr[other_genes, 7:12] + runif(length(other_genes), 1.8, 3.0)

z_expr <- t(scale(t(expr)))
pathway_scores <- t(sapply(signatures, function(g) colMeans(z_expr[g, , drop = FALSE])))

# ==============================================================================
# PLOT 1: Soft Pastel Heatmap (Saves as Plot1_Heatmap.png & shows in RStudio)
# ==============================================================================
pastel_heatmap <- colorRampPalette(c("#89CFF0", "#FFFFFF", "#F49AC2"))(50)

pheatmap(
  pathway_scores,
  main = "Regnase-1 Decay vs. IL-6 Anterograde Trafficking",
  cluster_cols = FALSE,
  color = pastel_heatmap,
  border_color = "#E0E0E0"
)

# Save high-res copy
pheatmap(
  pathway_scores,
  main = "Regnase-1 Decay vs. IL-6 Anterograde Trafficking",
  cluster_cols = FALSE,
  color = pastel_heatmap,
  border_color = "#E0E0E0",
  filename = "Plot1_Heatmap.png",
  width = 8, height = 5
)

# ==============================================================================
# PLOT 2: Negative Correlation Scatter Plot (Saves as Plot2_Correlation.png)
# ==============================================================================
df <- data.frame(
  Condition = rep(c("Proliferating", "Senescent"), each = 6),
  Regnase1_Score = pathway_scores["Regnase1_Decay_Brake", ],
  COPII_Trafficking_Score = pathway_scores["ER_Golgi_COPII", ],
  IL6_mRNA = expr["IL6", ]
)

cor_test <- cor.test(df$Regnase1_Score, df$COPII_Trafficking_Score)
r_val <- round(cor_test$estimate, 2)
p_val <- signif(cor_test$p.value, 2)

p_scatter <- ggplot(df, aes(x = Regnase1_Score, y = COPII_Trafficking_Score)) +
  geom_point(aes(color = Condition, size = IL6_mRNA), alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, fill = "#EAEAEA", color = "#555555", linetype = "dashed") +
  annotate("text", x = min(df$Regnase1_Score), y = max(df$COPII_Trafficking_Score),
           label = paste0("Pearson R = ", r_val, ", p = ", p_val),
           hjust = 0, vjust = 1, size = 4.5, fontface = "bold", color = "#333333") +
  scale_color_manual(values = c("Proliferating" = "#6BAED6", "Senescent" = "#FB6A4A")) +
  labs(
    title = "Negative Correlation: Regnase-1 Decay vs. IL-6 Anterograde Transport",
    subtitle = "Loss of ZC3H12A/UPF1/RC3H1 stabilizes IL6 mRNA and drives ER-to-Golgi flux",
    x = "Regnase-1 / UPF1 / Roquin-1 Activity (Z-score)",
    y = "ER-to-Golgi COPII Trafficking Activity (Z-score)",
    size = "IL6 mRNA (log2)"
  ) +
  theme_classic(base_size = 12)

print(p_scatter)
ggsave("Plot2_Correlation.png", plot = p_scatter, width = 8, height = 5.5, dpi = 300)

# ==============================================================================
# PLOT 3: Clean, Unmixed Pastel Network Diagram (Saves as Plot3_Network.png)
# ==============================================================================
edges <- data.frame(
  from = c("ZC3H12A", "ZC3H12A", "ZC3H12A", "IKBKB",   "STING1", "CGAS",   "RELA",
           "IL6",     "IL6",     "IL6R",    "SEC24D",  "SEC23A", "LMAN1",  "RAB1A",
           "ARF1",    "RAB8A",   "VAMP3",   "IL6",     "RELA"),
  to   = c("UPF1",    "RC3H1",   "IL6",     "ZC3H12A", "IKBKB",  "STING1", "IL6",
           "IL6R",    "SEC24D",  "IL6ST",   "SEC23A",  "SAR1A",  "SEC24D", "LMAN1",
           "RAB1A",   "ARF1",    "RAB8A",   "STAT3",   "CDKN1A")
)

g <- graph_from_data_frame(edges, directed = FALSE)

# Light pastel colors for clear readability
V(g)$color <- ifelse(V(g)$name %in% c("ZC3H12A", "UPF1", "RC3H1"), "#AED6F1",         # Soft Sky Blue
                     ifelse(V(g)$name %in% c("IKBKB", "CGAS", "STING1", "RELA", "IL6",
                                             "IL6R", "IL6ST", "STAT3", "CDKN1A"), "#F5B7B1",  # Soft Rose Pink
                            "#FAD7A0"))                                                       # Soft Peach

# Use Kamada-Kawai force-directed layout so nodes never overlap
set.seed(123)
unmixed_layout <- layout_with_kk(g)

# Helper function to draw clean network
draw_clean_network <- function() {
  par(mar = c(1, 1, 3, 1))
  plot(
    g,
    layout = unmixed_layout,
    vertex.size = 22,
    vertex.frame.color = "#7F8C8D",
    vertex.label.cex = 0.72,
    vertex.label.family = "sans",
    vertex.label.font = 2,
    vertex.label.color = "#1C2833",
    edge.width = 1.8,
    edge.color = "#BDC3C7",
    main = "Regnase-1 (Blue) | IKK/IL-6 Senescence (Pink) | Anterograde Transport (Peach)"
  )
  legend(
    "bottomleft",
    legend = c("Regnase-1 mRNA Decay", "IKK / Senescence / IL-6", "Anterograde Trafficking"),
    pch = 21,
    pt.bg = c("#AED6F1", "#F5B7B1", "#FAD7A0"),
    col = "#7F8C8D",
    pt.cex = 2,
    cex = 0.85,
    bty = "n"
  )
}

# Show in RStudio Plots pane
draw_clean_network()

# Save high-res PNG to your folder
png("Plot3_Network.png", width = 2600, height = 2000, res = 300)
draw_clean_network()
dev.off()

cat("\nDone! Saved all 3 plots to:", getwd(), "\n")
# ==============================================================================
# Focused RStudio Pipeline:
# IL-6 Intracellular Vesicular Anterograde Trafficking Pathways in Senescence
# ==============================================================================

library(igraph)
library(ggplot2)
library(pheatmap)

# 1. Core Anterograde Vesicular Trafficking Signatures (4 Sequential Stages)
il6_trafficking_stages <- list(
  Stage1_ER_Exit_COPII = c(
    "IL6", "IL6R", "IL6ST", "SAR1A", "SAR1B",
    "SEC23A", "SEC23B", "SEC24A", "SEC24B", "SEC24C", "SEC24D", "SEC13", "SEC31A"
  ),
  Stage2_ERGIC_to_CisGolgi = c(
    "LMAN1", "MCFD2", "RAB1A", "RAB1B", "USO1", "GOLGA2", "BET1", "SEC22B"
  ),
  Stage3_IntraGolgi_TGN_Sorting = c(
    "ARF1", "ARF4", "GBF1", "RAB6A", "GOLGA1", "GOLGA4", "SORT1", "AP1G1"
  ),
  Stage4_PostGolgi_Secretory_Intracrine = c(
    "RAB8A", "RAB11A", "VAMP3", "VAMP7", "SNAP23", "STX4", "EXOC2", "EXOC4"
  )
)

# 2. Simulate Expression Matrix (Proliferating vs. Senescent Cells)
set.seed(101)
all_genes <- unique(unlist(il6_trafficking_stages))
samples <- c(paste0("NonSenescent_", 1:6), paste0("Senescent_", 1:6))

expr <- matrix(
  rnorm(length(all_genes) * length(samples), mean = 6, sd = 0.65),
  nrow = length(all_genes),
  dimnames = list(all_genes, samples)
)

# Upregulate IL-6 cargo and anterograde vesicular machinery in Senescence
expr[, 7:12] <- expr[, 7:12] + runif(length(all_genes), 1.8, 3.0)

z_expr <- t(scale(t(expr)))
stage_scores <- t(sapply(il6_trafficking_stages, function(g) colMeans(z_expr[g, , drop = FALSE])))

# ==============================================================================
# PLOT 1: Pastel Heatmap of the 4 Anterograde Trafficking Stages
# ==============================================================================
pastel_palette <- colorRampPalette(c("#AED6F1", "#FFFFFF", "#F5B7B1"))(50)

pheatmap(
  stage_scores,
  main = "IL-6 Anterograde Vesicular Trafficking Stages in Senescence",
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  color = pastel_palette,
  border_color = "#EAECEE"
)

pheatmap(
  stage_scores,
  main = "IL-6 Anterograde Vesicular Trafficking Stages in Senescence",
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  color = pastel_palette,
  border_color = "#EAECEE",
  filename = "Plot1_IL6_Trafficking_Heatmap.png",
  width = 8.5, height = 4.8
)

# ==============================================================================
# PLOT 2: Coupling of ER-to-Golgi COPII Transport & Post-Golgi Vesicle Release
# ==============================================================================
df <- data.frame(
  State = rep(c("Non-Senescent", "Senescent"), each = 6),
  COPII_ER_Exit = stage_scores["Stage1_ER_Exit_COPII", ],
  PostGolgi_Exocytosis = stage_scores["Stage4_PostGolgi_Secretory_Intracrine", ],
  IL6_Level = expr["IL6", ]
)

cor_res <- cor.test(df$COPII_ER_Exit, df$PostGolgi_Exocytosis)
r_val <- round(cor_res$estimate, 2)
p_val <- signif(cor_res$p.value, 2)

p_scatter <- ggplot(df, aes(x = COPII_ER_Exit, y = PostGolgi_Exocytosis)) +
  geom_point(aes(color = State, size = IL6_Level), alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, fill = "#F2F4F4", color = "#566573", linetype = "dashed") +
  annotate("text", x = min(df$COPII_ER_Exit), y = max(df$PostGolgi_Exocytosis),
           label = paste0("Pearson R = ", r_val, ", p = ", p_val),
           hjust = 0, vjust = 1, size = 4.5, fontface = "bold", color = "#2C3E50") +
  scale_color_manual(values = c("Non-Senescent" = "#7FB3D5", "Senescent" = "#F1948A")) +
  labs(
    title = "Anterograde Flux: ER-to-Golgi COPII vs. Post-Golgi IL-6 Vesicles",
    subtitle = "Coordinated induction of COPII budding and SNARE/Rab exocytosis in senescence",
    x = "Stage 1: ER Exit & COPII Budding Activity (Z-score)",
    y = "Stage 4: Post-Golgi Vesicular Transport & Exocytosis (Z-score)",
    size = "IL6 Cargo (log2)"
  ) +
  theme_classic(base_size = 12)

print(p_scatter)
ggsave("Plot2_ER_to_Golgi_Coupling.png", plot = p_scatter, width = 8, height = 5.5, dpi = 300)

# ==============================================================================
# PLOT 3: Clean, Unmixed Pastel Network of IL-6 Anterograde Trafficking
# ==============================================================================
edges <- data.frame(
  from = c(
    # ER Lumen Cargo -> COPII Coat Assembly
    "IL6",    "IL6R",   "IL6",    "SAR1A",  "SEC23A", "SEC24D", "SEC13",
    # ER-Golgi Intermediate Compartment (ERGIC) -> Cis-Golgi
    "SEC24D", "LMAN1",  "LMAN1",  "RAB1A",  "USO1",
    # Intra-Golgi & Trans-Golgi Network (TGN) Sorting
    "GOLGA2", "ARF1",   "GBF1",   "RAB6A",
    # Post-Golgi Secretory / Intracrine Vesicles & SNARE Fusion
    "RAB6A",  "RAB8A",  "RAB11A", "VAMP3",  "SNAP23", "RAB8A"
  ),
  to   = c(
    "IL6R",   "IL6ST",  "SEC24D", "SEC23A", "SEC24D", "SEC13",  "SEC31A",
    "LMAN1",  "MCFD2",  "RAB1A",  "USO1",   "GOLGA2",
    "ARF1",   "GBF1",   "RAB6A",  "SORT1",
    "RAB8A",  "VAMP3",  "RAB8A",  "SNAP23", "STX4",   "EXOC4"
  )
)

g <- graph_from_data_frame(edges, directed = TRUE)

# Assign soft pastel colors by subcellular compartment
V(g)$color <- ifelse(V(g)$name %in% c("IL6", "IL6R", "IL6ST"), "#F5B7B1",                 # Soft Pink: IL-6 Cargo Complex
                     ifelse(V(g)$name %in% c("SAR1A", "SEC23A", "SEC24D", "SEC13", "SEC31A"), "#FAD7A0", # Soft Peach: ER COPII Budding
                            ifelse(V(g)$name %in% c("LMAN1", "MCFD2", "RAB1A", "USO1", "GOLGA2"), "#F9E79F",    # Soft Cream Yellow: ERGIC / Cis-Golgi
                                   ifelse(V(g)$name %in% c("ARF1", "GBF1", "RAB6A", "SORT1"), "#D7BDE2",               # Soft Lavender: TGN Sorting
                                          "#AED6F1"))))                                                                 # Soft Sky Blue: Post-Golgi / SNAREs

set.seed(42)
unmixed_layout <- layout_with_kk(g)

draw_il6_trafficking_network <- function() {
  par(mar = c(1, 1, 3, 1))
  plot(
    g,
    layout = unmixed_layout,
    vertex.size = 21,
    vertex.frame.color = "#85929E",
    vertex.label.cex = 0.72,
    vertex.label.family = "sans",
    vertex.label.font = 2,
    vertex.label.color = "#1C2833",
    edge.width = 1.7,
    edge.arrow.size = 0.45,
    edge.color = "#AEB6BF",
    main = "IL-6 Intracellular Vesicular Anterograde Trafficking Pathway"
  )
  legend(
    "bottomleft",
    legend = c(
      "IL-6 / IL-6R Cargo Complex",
      "ER Exit (COPII Coat)",
      "ERGIC & Cis-Golgi Tether",
      "Trans-Golgi Network (TGN)",
      "Post-Golgi Vesicles & SNAREs"
    ),
    pch = 21,
    pt.bg = c("#F5B7B1", "#FAD7A0", "#F9E79F", "#D7BDE2", "#AED6F1"),
    col = "#85929E",
    pt.cex = 1.9,
    cex = 0.8,
    bty = "n"
  )
}

draw_il6_trafficking_network()

png("Plot3_IL6_Anterograde_Network.png", width = 2800, height = 2100, res = 300)
draw_il6_trafficking_network()
dev.off()

cat("\nDone! All 3 IL-6 anterograde trafficking plots saved in:", getwd(), "\n")
# ==============================================================================
# Complete RStudio Pipeline:
# IL-6 Intracellular Vesicular Anterograde Trafficking Pathways in Senescence
# (Includes Heatmap, Scatter Coupling, Unmixed Network, & Directional Bar Plot)
# ==============================================================================

library(igraph)
library(ggplot2)
library(pheatmap)

# 1. Core Anterograde Vesicular Trafficking Signatures (4 Sequential Stages)
il6_trafficking_stages <- list(
  Stage1_ER_Exit_COPII = c(
    "IL6", "IL6R", "IL6ST", "SAR1A", "SAR1B",
    "SEC23A", "SEC23B", "SEC24A", "SEC24B", "SEC24C", "SEC24D", "SEC13", "SEC31A"
  ),
  Stage2_ERGIC_to_CisGolgi = c(
    "LMAN1", "MCFD2", "RAB1A", "RAB1B", "USO1", "GOLGA2", "BET1", "SEC22B"
  ),
  Stage3_IntraGolgi_TGN_Sorting = c(
    "ARF1", "ARF4", "GBF1", "RAB6A", "GOLGA1", "GOLGA4", "SORT1", "AP1G1"
  ),
  Stage4_PostGolgi_Secretory_Intracrine = c(
    "RAB8A", "RAB11A", "VAMP3", "VAMP7", "SNAP23", "STX4", "EXOC2", "EXOC4"
  )
)

# 2. Simulate Expression Matrix (Proliferating vs. Senescent Cells)
set.seed(101)
all_genes <- unique(unlist(il6_trafficking_stages))
samples <- c(paste0("NonSenescent_", 1:6), paste0("Senescent_", 1:6))

expr <- matrix(
  rnorm(length(all_genes) * length(samples), mean = 6, sd = 0.65),
  nrow = length(all_genes),
  dimnames = list(all_genes, samples)
)

# Upregulate IL-6 cargo and anterograde vesicular machinery in Senescence
expr[, 7:12] <- expr[, 7:12] + runif(length(all_genes), 1.8, 3.0)

z_expr <- t(scale(t(expr)))
stage_scores <- t(sapply(il6_trafficking_stages, function(g) colMeans(z_expr[g, , drop = FALSE])))

# ==============================================================================
# PLOT 1: Pastel Heatmap of the 4 Anterograde Trafficking Stages
# ==============================================================================
pastel_palette <- colorRampPalette(c("#AED6F1", "#FFFFFF", "#F5B7B1"))(50)

pheatmap(
  stage_scores,
  main = "IL-6 Anterograde Vesicular Trafficking Stages in Senescence",
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  color = pastel_palette,
  border_color = "#EAECEE"
)

pheatmap(
  stage_scores,
  main = "IL-6 Anterograde Vesicular Trafficking Stages in Senescence",
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  color = pastel_palette,
  border_color = "#EAECEE",
  filename = "Plot1_IL6_Trafficking_Heatmap.png",
  width = 8.5, height = 4.8
)

# ==============================================================================
# PLOT 2: Coupling of ER-to-Golgi COPII Transport & Post-Golgi Vesicle Release
# ==============================================================================
df <- data.frame(
  State = rep(c("Non-Senescent", "Senescent"), each = 6),
  COPII_ER_Exit = stage_scores["Stage1_ER_Exit_COPII", ],
  PostGolgi_Exocytosis = stage_scores["Stage4_PostGolgi_Secretory_Intracrine", ],
  IL6_Level = expr["IL6", ]
)

cor_res <- cor.test(df$COPII_ER_Exit, df$PostGolgi_Exocytosis)
r_val <- round(cor_res$estimate, 2)
p_val <- signif(cor_res$p.value, 2)

p_scatter <- ggplot(df, aes(x = COPII_ER_Exit, y = PostGolgi_Exocytosis)) +
  geom_point(aes(color = State, size = IL6_Level), alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, fill = "#F2F4F4", color = "#566573", linetype = "dashed") +
  annotate("text", x = min(df$COPII_ER_Exit), y = max(df$PostGolgi_Exocytosis),
           label = paste0("Pearson R = ", r_val, ", p = ", p_val),
           hjust = 0, vjust = 1, size = 4.5, fontface = "bold", color = "#2C3E50") +
  scale_color_manual(values = c("Non-Senescent" = "#7FB3D5", "Senescent" = "#F1948A")) +
  labs(
    title = "Anterograde Flux: ER-to-Golgi COPII vs. Post-Golgi IL-6 Vesicles",
    subtitle = "Coordinated induction of COPII budding and SNARE/Rab exocytosis in senescence",
    x = "Stage 1: ER Exit & COPII Budding Activity (Z-score)",
    y = "Stage 4: Post-Golgi Vesicular Transport & Exocytosis (Z-score)",
    size = "IL6 Cargo (log2)"
  ) +
  theme_classic(base_size = 12)

print(p_scatter)
ggsave("Plot2_ER_to_Golgi_Coupling.png", plot = p_scatter, width = 8, height = 5.5, dpi = 300)

# ==============================================================================
# PLOT 3: Clean, Unmixed Pastel Network of IL-6 Anterograde Trafficking
# ==============================================================================
edges <- data.frame(
  from = c(
    "IL6",    "IL6R",   "IL6",    "SAR1A",  "SEC23A", "SEC24D", "SEC13",
    "SEC24D", "LMAN1",  "LMAN1",  "RAB1A",  "USO1",
    "GOLGA2", "ARF1",   "GBF1",   "RAB6A",
    "RAB6A",  "RAB8A",  "RAB11A", "VAMP3",  "SNAP23", "RAB8A"
  ),
  to   = c(
    "IL6R",   "IL6ST",  "SEC24D", "SEC23A", "SEC24D", "SEC13",  "SEC31A",
    "LMAN1",  "MCFD2",  "RAB1A",  "USO1",   "GOLGA2",
    "ARF1",   "GBF1",   "RAB6A",  "SORT1",
    "RAB8A",  "VAMP3",  "RAB8A",  "SNAP23", "STX4",   "EXOC4"
  )
)

g <- graph_from_data_frame(edges, directed = TRUE)

V(g)$color <- ifelse(V(g)$name %in% c("IL6", "IL6R", "IL6ST"), "#F5B7B1",
                     ifelse(V(g)$name %in% c("SAR1A", "SEC23A", "SEC24D", "SEC13", "SEC31A"), "#FAD7A0",
                            ifelse(V(g)$name %in% c("LMAN1", "MCFD2", "RAB1A", "USO1", "GOLGA2"), "#F9E79F",
                                   ifelse(V(g)$name %in% c("ARF1", "GBF1", "RAB6A", "SORT1"), "#D7BDE2",
                                          "#AED6F1"))))

set.seed(42)
unmixed_layout <- layout_with_kk(g)

draw_il6_trafficking_network <- function() {
  par(mar = c(1, 1, 3, 1))
  plot(
    g,
    layout = unmixed_layout,
    vertex.size = 21,
    vertex.frame.color = "#85929E",
    vertex.label.cex = 0.72,
    vertex.label.family = "sans",
    vertex.label.font = 2,
    vertex.label.color = "#1C2833",
    edge.width = 1.7,
    edge.arrow.size = 0.45,
    edge.color = "#AEB6BF",
    main = "IL-6 Intracellular Vesicular Anterograde Trafficking Pathway"
  )
  legend(
    "bottomleft",
    legend = c(
      "IL-6 / IL-6R Cargo Complex",
      "ER Exit (COPII Coat)",
      "ERGIC & Cis-Golgi Tether",
      "Trans-Golgi Network (TGN)",
      "Post-Golgi Vesicles & SNAREs"
    ),
    pch = 21,
    pt.bg = c("#F5B7B1", "#FAD7A0", "#F9E79F", "#D7BDE2", "#AED6F1"),
    col = "#85929E",
    pt.cex = 1.9,
    cex = 0.8,
    bty = "n"
  )
}

draw_il6_trafficking_network()

png("Plot3_IL6_Anterograde_Network.png", width = 2800, height = 2100, res = 300)
draw_il6_trafficking_network()
dev.off()

# ==============================================================================
# PLOT 4: Directional Pathway Enrichment Bar Plot (Senescent vs. Non-Senescent)
# ==============================================================================
enrichment_list <- lapply(rownames(stage_scores), function(stage) {
  non_sen <- stage_scores[stage, 1:6]
  sen     <- stage_scores[stage, 7:12]
  t_res   <- t.test(sen, non_sen)
  delta_z <- mean(sen) - mean(non_sen)
  p_val   <- t_res$p.value
  # Directional score: sign(Delta Z) * -log10(p-value)
  dir_score <- sign(delta_z) * (-log10(p_val))
  
  data.frame(
    Stage = stage,
    Delta_Z = delta_z,
    P_Value = p_val,
    Directional_Score = dir_score
  )
})

enrich_df <- do.call(rbind, enrichment_list)

# Clean readable stage labels ordered sequentially from ER -> Plasma Membrane
enrich_df$Clean_Label <- factor(
  c(
    "Stage 1: ER Exit & COPII Budding",
    "Stage 2: ERGIC to Cis-Golgi Tethering",
    "Stage 3: Intra-Golgi & TGN Sorting",
    "Stage 4: Post-Golgi Vesicles & Exocytosis"
  ),
  levels = rev(c(
    "Stage 1: ER Exit & COPII Budding",
    "Stage 2: ERGIC to Cis-Golgi Tethering",
    "Stage 3: Intra-Golgi & TGN Sorting",
    "Stage 4: Post-Golgi Vesicles & Exocytosis"
  ))
)

enrich_df$Direction <- ifelse(enrich_df$Directional_Score >= 0, "Activated in Senescence", "Suppressed in Senescence")
enrich_df$Sig_Stars <- ifelse(enrich_df$P_Value < 0.001, "***",
                              ifelse(enrich_df$P_Value < 0.01, "**",
                                     ifelse(enrich_df$P_Value < 0.05, "*", "ns")))

p_bar <- ggplot(enrich_df, aes(x = Directional_Score, y = Clean_Label, fill = Direction)) +
  geom_col(width = 0.62, color = "#7F8C8D", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "#2C3E50", linewidth = 0.6) +
  geom_vline(xintercept = c(-1.301, 1.301), linetype = "dashed", color = "#95A5A6", linewidth = 0.5) +
  geom_text(
    aes(label = paste0("ΔZ = +", round(Delta_Z, 2), " (", Sig_Stars, ")")),
    hjust = ifelse(enrich_df$Directional_Score >= 0, -0.08, 1.08),
    size = 3.8, fontface = "bold", color = "#2C3E50"
  ) +
  scale_fill_manual(values = c(
    "Activated in Senescence" = "#F5B7B1",
    "Suppressed in Senescence" = "#AED6F1"
  )) +
  xlim(
    min(0, min(enrich_df$Directional_Score) * 1.35),
    max(enrich_df$Directional_Score) * 1.35
  ) +
  labs(
    title = "Directional Pathway Enrichment Across IL-6 Anterograde Trafficking",
    subtitle = "Signed -log10(p-value) comparing Senescent vs. Non-Senescent states (dashed line: p = 0.05)",
    x = "Directional Enrichment Score  [ sign(ΔZ) × -log10(p) ]",
    y = NULL,
    fill = "Pathway State"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(face = "bold", color = "#1C2833", size = 10.5)
  )

print(p_bar)
ggsave("Plot4_Directional_Enrichment_Barplot.png", plot = p_bar, width = 9.5, height = 5, dpi = 300)

cat("\nDone! All 4 plots saved in:", getwd(), "\n")
# ==============================================================================
# Complete RStudio Pipeline:
# IL-6 Intracellular Vesicular Anterograde Trafficking Pathways in Senescence
# (Includes Heatmap, Scatter Coupling, Network, Bar Plot, & Gene Lollipop Plot)
# ==============================================================================

library(igraph)
library(ggplot2)
library(pheatmap)

# 1. Core Anterograde Vesicular Trafficking Signatures (4 Sequential Stages)
il6_trafficking_stages <- list(
  Stage1_ER_Exit_COPII = c(
    "IL6", "IL6R", "IL6ST", "SAR1A", "SAR1B",
    "SEC23A", "SEC23B", "SEC24A", "SEC24B", "SEC24C", "SEC24D", "SEC13", "SEC31A"
  ),
  Stage2_ERGIC_to_CisGolgi = c(
    "LMAN1", "MCFD2", "RAB1A", "RAB1B", "USO1", "GOLGA2", "BET1", "SEC22B"
  ),
  Stage3_IntraGolgi_TGN_Sorting = c(
    "ARF1", "ARF4", "GBF1", "RAB6A", "GOLGA1", "GOLGA4", "SORT1", "AP1G1"
  ),
  Stage4_PostGolgi_Secretory_Intracrine = c(
    "RAB8A", "RAB11A", "VAMP3", "VAMP7", "SNAP23", "STX4", "EXOC2", "EXOC4"
  )
)

# 2. Simulate Expression Matrix (Non-Senescent vs. Senescent Cells)
set.seed(101)
all_genes <- unique(unlist(il6_trafficking_stages))
samples <- c(paste0("NonSenescent_", 1:6), paste0("Senescent_", 1:6))

expr <- matrix(
  rnorm(length(all_genes) * length(samples), mean = 6, sd = 0.65),
  nrow = length(all_genes),
  dimnames = list(all_genes, samples)
)

# Upregulate IL-6 cargo and anterograde vesicular machinery in Senescence
expr[, 7:12] <- expr[, 7:12] + runif(length(all_genes), 1.8, 3.0)

z_expr <- t(scale(t(expr)))
stage_scores <- t(sapply(il6_trafficking_stages, function(g) colMeans(z_expr[g, , drop = FALSE])))

# ==============================================================================
# PLOT 1: Heatmap of the 4 Anterograde Trafficking Stages
# ==============================================================================
pastel_palette <- colorRampPalette(c("#AED6F1", "#FFFFFF", "#F5B7B1"))(50)

pheatmap(
  stage_scores,
  main = "IL-6 Anterograde Vesicular Trafficking Stages in Senescence",
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  color = pastel_palette,
  border_color = "#EAECEE"
)

pheatmap(
  stage_scores,
  main = "IL-6 Anterograde Vesicular Trafficking Stages in Senescence",
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  color = pastel_palette,
  border_color = "#EAECEE",
  filename = "Plot1_IL6_Trafficking_Heatmap.png",
  width = 8.5, height = 4.8
)

# ==============================================================================
# PLOT 2: Coupling of ER-to-Golgi COPII Transport & Post-Golgi Vesicle Release
# ==============================================================================
df <- data.frame(
  State = rep(c("Non-Senescent", "Senescent"), each = 6),
  COPII_ER_Exit = stage_scores["Stage1_ER_Exit_COPII", ],
  PostGolgi_Exocytosis = stage_scores["Stage4_PostGolgi_Secretory_Intracrine", ],
  IL6_Level = expr["IL6", ]
)

cor_res <- cor.test(df$COPII_ER_Exit, df$PostGolgi_Exocytosis)
r_val <- round(cor_res$estimate, 2)
p_val <- signif(cor_res$p.value, 2)

p_scatter <- ggplot(df, aes(x = COPII_ER_Exit, y = PostGolgi_Exocytosis)) +
  geom_point(aes(color = State, size = IL6_Level), alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, fill = "#F2F4F4", color = "#566573", linetype = "dashed") +
  annotate("text", x = min(df$COPII_ER_Exit), y = max(df$PostGolgi_Exocytosis),
           label = paste0("Pearson R = ", r_val, ", p = ", p_val),
           hjust = 0, vjust = 1, size = 4.5, fontface = "bold", color = "#2C3E50") +
  scale_color_manual(values = c("Non-Senescent" = "#377EB8", "Senescent" = "#E41A1C")) +
  labs(
    title = "Anterograde Flux: ER-to-Golgi COPII vs. Post-Golgi IL-6 Vesicles",
    subtitle = "Coordinated induction of COPII budding and SNARE/Rab exocytosis in senescence",
    x = "Stage 1: ER Exit & COPII Budding Activity (Z-score)",
    y = "Stage 4: Post-Golgi Vesicular Transport & Exocytosis (Z-score)",
    size = "IL6 Cargo (log2)"
  ) +
  theme_classic(base_size = 12)

print(p_scatter)
ggsave("Plot2_ER_to_Golgi_Coupling.png", plot = p_scatter, width = 8, height = 5.5, dpi = 300)

# ==============================================================================
# PLOT 3: Clean, Unmixed Network of IL-6 Anterograde Trafficking
# ==============================================================================
edges <- data.frame(
  from = c(
    "IL6",    "IL6R",   "IL6",    "SAR1A",  "SEC23A", "SEC24D", "SEC13",
    "SEC24D", "LMAN1",  "LMAN1",  "RAB1A",  "USO1",
    "GOLGA2", "ARF1",   "GBF1",   "RAB6A",
    "RAB6A",  "RAB8A",  "RAB11A", "VAMP3",  "SNAP23", "RAB8A"
  ),
  to   = c(
    "IL6R",   "IL6ST",  "SEC24D", "SEC23A", "SEC24D", "SEC13",  "SEC31A",
    "LMAN1",  "MCFD2",  "RAB1A",  "USO1",   "GOLGA2",
    "ARF1",   "GBF1",   "RAB6A",  "SORT1",
    "RAB8A",  "VAMP3",  "RAB8A",  "SNAP23", "STX4",   "EXOC4"
  )
)

g <- graph_from_data_frame(edges, directed = TRUE)

V(g)$color <- ifelse(V(g)$name %in% c("IL6", "IL6R", "IL6ST"), "#F5B7B1",
                     ifelse(V(g)$name %in% c("SAR1A", "SEC23A", "SEC24D", "SEC13", "SEC31A"), "#FAD7A0",
                            ifelse(V(g)$name %in% c("LMAN1", "MCFD2", "RAB1A", "USO1", "GOLGA2"), "#F9E79F",
                                   ifelse(V(g)$name %in% c("ARF1", "GBF1", "RAB6A", "SORT1"), "#D7BDE2",
                                          "#AED6F1"))))

set.seed(42)
unmixed_layout <- layout_with_kk(g)

draw_il6_trafficking_network <- function() {
  par(mar = c(1, 1, 3, 1))
  plot(
    g,
    layout = unmixed_layout,
    vertex.size = 21,
    vertex.frame.color = "#85929E",
    vertex.label.cex = 0.72,
    vertex.label.family = "sans",
    vertex.label.font = 2,
    vertex.label.color = "#1C2833",
    edge.width = 1.7,
    edge.arrow.size = 0.45,
    edge.color = "#AEB6BF",
    main = "IL-6 Intracellular Vesicular Anterograde Trafficking Pathway"
  )
  legend(
    "bottomleft",
    legend = c(
      "IL-6 / IL-6R Cargo Complex",
      "ER Exit (COPII Coat)",
      "ERGIC & Cis-Golgi Tether",
      "Trans-Golgi Network (TGN)",
      "Post-Golgi Vesicles & SNAREs"
    ),
    pch = 21,
    pt.bg = c("#F5B7B1", "#FAD7A0", "#F9E79F", "#D7BDE2", "#AED6F1"),
    col = "#85929E",
    pt.cex = 1.9,
    cex = 0.8,
    bty = "n"
  )
}

draw_il6_trafficking_network()

png("Plot3_IL6_Anterograde_Network.png", width = 2800, height = 2100, res = 300)
draw_il6_trafficking_network()
dev.off()

# ==============================================================================
# PLOT 4: Directional Pathway Enrichment Bar Plot
# ==============================================================================
enrichment_list <- lapply(rownames(stage_scores), function(stage) {
  non_sen <- stage_scores[stage, 1:6]
  sen     <- stage_scores[stage, 7:12]
  t_res   <- t.test(sen, non_sen)
  delta_z <- mean(sen) - mean(non_sen)
  p_val   <- t_res$p.value
  dir_score <- sign(delta_z) * (-log10(p_val))
  
  data.frame(
    Stage = stage,
    Delta_Z = delta_z,
    P_Value = p_val,
    Directional_Score = dir_score
  )
})

enrich_df <- do.call(rbind, enrichment_list)

enrich_df$Clean_Label <- factor(
  c(
    "Stage 1: ER Exit & COPII Budding",
    "Stage 2: ERGIC to Cis-Golgi Tethering",
    "Stage 3: Intra-Golgi & TGN Sorting",
    "Stage 4: Post-Golgi Vesicles & Exocytosis"
  ),
  levels = rev(c(
    "Stage 1: ER Exit & COPII Budding",
    "Stage 2: ERGIC to Cis-Golgi Tethering",
    "Stage 3: Intra-Golgi & TGN Sorting",
    "Stage 4: Post-Golgi Vesicles & Exocytosis"
  ))
)

enrich_df$Direction <- ifelse(enrich_df$Directional_Score >= 0, "Activated in Senescence", "Suppressed in Senescence")
enrich_df$Sig_Stars <- ifelse(enrich_df$P_Value < 0.001, "***",
                              ifelse(enrich_df$P_Value < 0.01, "**",
                                     ifelse(enrich_df$P_Value < 0.05, "*", "ns")))

p_bar <- ggplot(enrich_df, aes(x = Directional_Score, y = Clean_Label, fill = Direction)) +
  geom_col(width = 0.62, color = "#333333", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "#1C2833", linewidth = 0.6) +
  geom_vline(xintercept = c(-1.301, 1.301), linetype = "dashed", color = "#7F8C8D", linewidth = 0.5) +
  geom_text(
    aes(label = paste0("ΔZ = +", round(Delta_Z, 2), " (", Sig_Stars, ")")),
    hjust = ifelse(enrich_df$Directional_Score >= 0, -0.08, 1.08),
    size = 3.8, fontface = "bold", color = "#1C2833"
  ) +
  scale_fill_manual(values = c(
    "Activated in Senescence" = "#E41A1C",
    "Suppressed in Senescence" = "#377EB8"
  )) +
  xlim(
    min(0, min(enrich_df$Directional_Score) * 1.35),
    max(enrich_df$Directional_Score) * 1.35
  ) +
  labs(
    title = "Directional Pathway Enrichment Across IL-6 Anterograde Trafficking",
    subtitle = "Signed -log10(p-value) comparing Senescent vs. Non-Senescent states (dashed line: p = 0.05)",
    x = "Directional Enrichment Score  [ sign(ΔZ) × -log10(p) ]",
    y = NULL,
    fill = "Pathway State"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(face = "bold", color = "#1C2833", size = 10.5)
  )

print(p_bar)
ggsave("Plot4_Directional_Enrichment_Barplot.png", plot = p_bar, width = 9.5, height = 5, dpi = 300)

# ==============================================================================
# PLOT 5: Gene-Level Lollipop Plot (COPII / ERGIC / Golgi / Rab-SNARE Drivers)
# Bright Crimson Red (#E41A1C) & Royal Blue (#377EB8) - No Pink
# ==============================================================================
stage_clean_names <- c(
  Stage1_ER_Exit_COPII                  = "Stage 1: ER Exit & COPII",
  Stage2_ERGIC_to_CisGolgi              = "Stage 2: ERGIC to Cis-Golgi",
  Stage3_IntraGolgi_TGN_Sorting         = "Stage 3: Intra-Golgi & TGN",
  Stage4_PostGolgi_Secretory_Intracrine = "Stage 4: Post-Golgi Rab & SNARE"
)

gene_stats_list <- lapply(names(il6_trafficking_stages), function(st) {
  genes_in_stage <- il6_trafficking_stages[[st]]
  do.call(rbind, lapply(genes_in_stage, function(g) {
    non_sen_vals <- expr[g, 1:6]
    sen_vals     <- expr[g, 7:12]
    # ==============================================================================
    # RStudio Pipeline: 4-Layer SASP Regulatory Crosstalk Mapping
    # Layer 1: cGAS-STING & DDR Sensors
    # Layer 2: NF-kB / C/EBPb / JAK-STAT Transcription
    # Layer 3: mTORC1 & RNA-Binding Protein (RBP) Stability Switch
    # Layer 4: Anterograde Vesicular & Exosome Trafficking
    # ==============================================================================
    
    # Auto-install required CRAN packages if missing
    pkgs <- c("igraph", "ggplot2", "pheatmap", "reshape2")
    new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
    if (length(new_pkgs) > 0) install.packages(new_pkgs)
    
    library(igraph)
    library(ggplot2)
    library(pheatmap)
    library(reshape2)
    
    # 1. Define the 4 SASP Regulatory Layer Signatures
    sasp_4layers <- list(
      L1_Sensors_cGAS_STING_DDR = c(
        "CGAS", "STING1", "TBK1", "IRF3", "ATM", "ATR", "CHEK2", "MAPK14"
      ),
      L2_Transcription_NFkB_STAT3 = c(
        "RELA", "NFKB1", "CEBPB", "STAT3", "JAK1", "JAK2", "BRD4", "IL6", "CXCL8", "IL1A"
      ),
      L3_PostTranscriptional_mTOR_RBP = c(
        "MTOR", "RPTOR", "RPS6KB1", "MAPKAPK2", "IKBKB", "ZC3H12A", "UPF1", "ZFP36", "RC3H1"
      ),
      L4_Trafficking_Anterograde_Exosome = c(
        "SEC24D", "SAR1A", "LMAN1", "RAB1A", "ARF1", "RAB6A", "RAB8A", "VAMP3", "SNAP23", "RAB27A", "PDCD6IP"
      )
    )
    
    # 2. Simulate Multi-Stage Expression Matrix (Proliferating -> Early SASP -> Mature SASP)
    set.seed(2026)
    all_genes <- unique(unlist(sasp_4layers))
    samples <- c(
      paste0("Proliferating_", 1:4),
      paste0("Early_Senescent_", 1:4),
      paste0("Mature_SASP_", 1:4)
    )
    
    expr <- matrix(
      rnorm(length(all_genes) * length(samples), mean = 6.0, sd = 0.55),
      nrow = length(all_genes),
      dimnames = list(all_genes, samples)
    )
    
    # Progressive activation across senescence stages
    rbp_brakes <- c("ZC3H12A", "UPF1", "ZFP36", "RC3H1") # Degraded/inhibited during mature SASP
    active_sasp <- setdiff(all_genes, rbp_brakes)
    
    # Early senescence: Sensors (L1) & initial transcription rise; RBP brakes begin falling
    expr[active_sasp, 5:8]  <- expr[active_sasp, 5:8]  + runif(length(active_sasp), 0.9, 1.6)
    expr[rbp_brakes, 5:8]   <- expr[rbp_brakes, 5:8]   - runif(length(rbp_brakes), 0.8, 1.4)
    
    # Mature SASP: Full cascade across all 4 layers; RBP brakes strongly suppressed
    expr[active_sasp, 9:12] <- expr[active_sasp, 9:12] + runif(length(active_sasp), 2.2, 3.4)
    expr[rbp_brakes, 9:12]  <- expr[rbp_brakes, 9:12]  - runif(length(rbp_brakes), 2.0, 3.0)
    
    # Compute Z-score activity per layer (inverting RBP brakes inside Layer 3 so score reflects net SASP licensing)
    z_expr <- t(scale(t(expr)))
    z_expr_oriented <- z_expr
    z_expr_oriented[rbp_brakes, ] <- -z_expr_oriented[rbp_brakes, ]
    
    layer_scores <- t(sapply(sasp_4layers, function(g) colMeans(z_expr_oriented[g, , drop = FALSE])))
    rownames(layer_scores) <- c(
      "Layer 1: cGAS-STING / DDR Sensors",
      "Layer 2: NF-kB / C/EBPb / STAT3",
      "Layer 3: mTORC1 & RBP Derepression",
      "Layer 4: Anterograde & Exosome Release"
    )
    
    # ==============================================================================
    # PLOT 1: 4-Layer SASP Cascade Heatmap
    # ==============================================================================
    heat_pal <- colorRampPalette(c("#377EB8", "#FFFFFF", "#E41A1C"))(60)
    
    pheatmap(
      layer_scores,
      main = "4-Layer SASP Regulatory Cascade Across Senescence Progression",
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      color = heat_pal,
      border_color = "#E5E7E9",
      fontsize = 10.5
    )
    
    pheatmap(
      layer_scores,
      main = "4-Layer SASP Regulatory Cascade Across Senescence Progression",
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      color = heat_pal,
      border_color = "#E5E7E9",
      fontsize = 10.5,
      filename = "SASP_Plot1_4Layer_Heatmap.png",
      width = 9.5, height = 4.8
    )
    
    # ==============================================================================
    # PLOT 2: Inter-Layer Crosstalk Correlation Matrix
    # ==============================================================================
    cor_mat <- cor(t(layer_scores), method = "pearson")
    cor_melted <- melt(cor_mat)
    
    p_cormat <- ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
      geom_tile(color = "white", linewidth = 1.2) +
      geom_text(aes(label = sprintf("R = %.2f", value)), color = "black", fontface = "bold", size = 4.2) +
      scale_fill_gradient2(
        low = "#377EB8", mid = "#F8F9F9", high = "#E41A1C",
        midpoint = 0, limit = c(-1, 1), name = "Pearson R"
      ) +
      labs(
        title = "Inter-Layer SASP Crosstalk Coupling Matrix",
        subtitle = "Coordinated feed-forward amplification across sensing, transcription, mRNA stability & secretion",
        x = NULL, y = NULL
      ) +
      theme_minimal(base_size = 11.5) +
      theme(
        axis.text.x = element_text(angle = 25, hjust = 1, face = "bold", color = "#1C2833"),
        axis.text.y = element_text(face = "bold", color = "#1C2833"),
        panel.grid = element_blank()
      )
    
    print(p_cormat)
    ggsave("SASP_Plot2_Crosstalk_Correlation_Matrix.png", plot = p_cormat, width = 8.5, height = 6.5, dpi = 300)
    
    # ==============================================================================
    # PLOT 3: Unmixed 4-Layer SASP Crosstalk Signaling Network
    # ==============================================================================
    crosstalk_edges <- data.frame(
      from = c(
        # Layer 1 Internal & L1 -> L2 Crosstalk
        "CGAS",    "STING1",  "STING1",  "ATM",     "MAPK14",
        # Layer 1 -> L3 Crosstalk (IKK/MK2 inhibit RBP brakes)
        "STING1",  "MAPK14",  "IKBKB",   "MAPKAPK2",
        # Layer 2 Internal & Autocrine Loop
        "RELA",    "CEBPB",   "RELA",    "IL1A",    "IL6",     "JAK1",    "STAT3",   "BRD4",
        # Layer 3 mTOR / RBP -> IL6/IL1A Stabilization & L4 Translation
        "MTOR",    "MTOR",    "ZC3H12A", "ZFP36",   "MTOR",
        # Layer 2/3 -> Layer 4 Anterograde & Exosome Trafficking
        "IL6",     "SEC24D",  "LMAN1",   "RAB1A",   "ARF1",    "RAB6A",   "RAB8A",   "VAMP3",
        "IL6",     "RAB27A"
      ),
      to = c(
        "STING1",  "TBK1",    "RELA",    "RELA",    "MAPKAPK2",
        "IKBKB",   "MTOR",    "ZC3H12A", "ZFP36",
        "IL6",     "IL6",     "CXCL8",   "RELA",    "JAK1",    "STAT3",   "IL6",     "RELA",
        "IL1A",    "MAPKAPK2","IL6",     "CXCL8",   "SEC24D",
        "SEC24D",  "LMAN1",   "RAB1A",   "ARF1",    "RAB6A",   "RAB8A",   "VAMP3",   "SNAP23",
        "RAB27A",  "PDCD6IP"
      )
    )
    
    g <- graph_from_data_frame(crosstalk_edges, directed = TRUE)
    
    # Bright distinct colors per SASP regulatory layer
    V(g)$color <- ifelse(V(g)$name %in% sasp_4layers$L1_Sensors_cGAS_STING_DDR, "#377EB8",        # Bright Royal Blue (Layer 1)
                         ifelse(V(g)$name %in% sasp_4layers$L2_Transcription_NFkB_STAT3, "#E41A1C",      # Bright Crimson Red (Layer 2)
                                ifelse(V(g)$name %in% sasp_4layers$L3_PostTranscriptional_mTOR_RBP, "#FF7F00",  # Bright Amber Orange (Layer 3)
                                       "#4DAF4A")))                                                              # Bright Emerald Green (Layer 4)
    
    set.seed(77)
    kk_layout <- layout_with_kk(g)
    
    draw_4layer_network <- function() {
      par(mar = c(1, 1, 3, 1))
      plot(
        g,
        layout = kk_layout,
        vertex.size = 19,
        vertex.frame.color = "#1C2833",
        vertex.label.cex = 0.68,
        vertex.label.family = "sans",
        vertex.label.font = 2,
        vertex.label.color = "white",
        edge.width = 1.6,
        edge.arrow.size = 0.42,
        edge.color = "#7F8C8D",
        main = "4-Layer SASP Regulatory Crosstalk Network"
      )
      legend(
        "bottomleft",
        legend = c(
          "Layer 1: cGAS-STING & DDR Sensors",
          "Layer 2: NF-kB / C/EBPb / JAK-STAT3",
          "Layer 3: mTORC1 & RBP Stability Switch",
          "Layer 4: Anterograde & Exosome Trafficking"
        ),
        pch = 21,
        pt.bg = c("#377EB8", "#E41A1C", "#FF7F00", "#4DAF4A"),
        col = "#1C2833",
        pt.cex = 1.9,
        cex = 0.8,
        bty = "n"
      )
    }
    
    draw_4layer_network()
    
    png("SASP_Plot3_4Layer_Crosstalk_Network.png", width = 3000, height = 2200, res = 300)
    draw_4layer_network()
    dev.off()
    
    # ==============================================================================
    # PLOT 4: Crosstalk Bridge Centrality Lollipop Plot (Top Hub Regulators)
    # ==============================================================================
    bet_vals <- betweenness(g, directed = FALSE, normalized = TRUE)
    hub_df <- data.frame(
      Gene = names(bet_vals),
      Centrality = as.numeric(bet_vals)
    )
    
    hub_df$Layer <- ifelse(hub_df$Gene %in% sasp_4layers$L1_Sensors_cGAS_STING_DDR, "Layer 1: Sensors",
                           ifelse(hub_df$Gene %in% sasp_4layers$L2_Transcription_NFkB_STAT3, "Layer 2: Transcription",
                                  ifelse(hub_df$Gene %in% sasp_4layers$L3_PostTranscriptional_mTOR_RBP, "Layer 3: mTOR / RBP",
                                         "Layer 4: Trafficking")))
    
    hub_df <- hub_df[order(hub_df$Centrality), ]
    top_hubs <- tail(hub_df, 14)
    top_hubs$Gene <- factor(top_hubs$Gene, levels = top_hubs$Gene)
    
    p_hubs <- ggplot(top_hubs, aes(x = Centrality, y = Gene, color = Layer)) +
      geom_segment(aes(x = 0, xend = Centrality, y = Gene, yend = Gene), linewidth = 1.2) +
      geom_point(size = 5) +
      scale_color_manual(values = c(
        "Layer 1: Sensors"       = "#377EB8",
        "Layer 2: Transcription" = "#E41A1C",
        "Layer 3: mTOR / RBP"    = "#FF7F00",
        "Layer 4: Trafficking"   = "#4DAF4A"
      )) +
      labs(
        title = "Key Inter-Layer Bridge Hubs Connecting the 4 SASP Modules",
        subtitle = "Normalized network betweenness centrality across sensing, transcription, stability & secretion",
        x = "Betweenness Centrality (Crosstalk Bridge Score)",
        y = NULL,
        color = "SASP Layer"
      ) +
      theme_classic(base_size = 12) +
      theme(
        axis.text.y = element_text(face = "bold", color = "#1C2833", size = 10.5),
        legend.position = "bottom"
      )
    
    print(p_hubs)
    ggsave("SASP_Plot4_Hub_Centrality_Lollipop.png", plot = p_hubs, width = 9, height = 6, dpi = 300)
    
    cat("\nDone! All 4 SASP crosstalk plots saved to:", getwd(), "\n")
    # ==============================================================================
    # RStudio Pipeline: 4-Layer SASP Regulatory Crosstalk Mapping
    # Layer 1: cGAS-STING & DDR Sensors
    # Layer 2: NF-kB / C/EBPb / JAK-STAT Transcription
    # Layer 3: mTORC1 & RNA-Binding Protein (RBP) Stability Switch
    # Layer 4: Anterograde Vesicular & Exosome Trafficking
    # ==============================================================================
    
    # Auto-install required CRAN packages if missing
    pkgs <- c("igraph", "ggplot2", "pheatmap", "reshape2")
    new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
    if (length(new_pkgs) > 0) install.packages(new_pkgs)
    
    library(igraph)
    library(ggplot2)
    library(pheatmap)
    library(reshape2)
    
    # 1. Define the 4 SASP Regulatory Layer Signatures
    sasp_4layers <- list(
      L1_Sensors_cGAS_STING_DDR = c(
        "CGAS", "STING1", "TBK1", "IRF3", "ATM", "ATR", "CHEK2", "MAPK14"
      ),
      L2_Transcription_NFkB_STAT3 = c(
        "RELA", "NFKB1", "CEBPB", "STAT3", "JAK1", "JAK2", "BRD4", "IL6", "CXCL8", "IL1A"
      ),
      L3_PostTranscriptional_mTOR_RBP = c(
        "MTOR", "RPTOR", "RPS6KB1", "MAPKAPK2", "IKBKB", "ZC3H12A", "UPF1", "ZFP36", "RC3H1"
      ),
      L4_Trafficking_Anterograde_Exosome = c(
        "SEC24D", "SAR1A", "LMAN1", "RAB1A", "ARF1", "RAB6A", "RAB8A", "VAMP3", "SNAP23", "RAB27A", "PDCD6IP"
      )
    )
    
    # 2. Simulate Multi-Stage Expression Matrix (Proliferating -> Early SASP -> Mature SASP)
    set.seed(2026)
    all_genes <- unique(unlist(sasp_4layers))
    samples <- c(
      paste0("Proliferating_", 1:4),
      paste0("Early_Senescent_", 1:4),
      paste0("Mature_SASP_", 1:4)
    )
    
    expr <- matrix(
      rnorm(length(all_genes) * length(samples), mean = 6.0, sd = 0.55),
      nrow = length(all_genes),
      dimnames = list(all_genes, samples)
    )
    
    # Progressive activation across senescence stages
    rbp_brakes <- c("ZC3H12A", "UPF1", "ZFP36", "RC3H1") # Degraded/inhibited during mature SASP
    active_sasp <- setdiff(all_genes, rbp_brakes)
    
    # Early senescence: Sensors (L1) & initial transcription rise; RBP brakes begin falling
    expr[active_sasp, 5:8]  <- expr[active_sasp, 5:8]  + runif(length(active_sasp), 0.9, 1.6)
    expr[rbp_brakes, 5:8]   <- expr[rbp_brakes, 5:8]   - runif(length(rbp_brakes), 0.8, 1.4)
    
    # Mature SASP: Full cascade across all 4 layers; RBP brakes strongly suppressed
    expr[active_sasp, 9:12] <- expr[active_sasp, 9:12] + runif(length(active_sasp), 2.2, 3.4)
    expr[rbp_brakes, 9:12]  <- expr[rbp_brakes, 9:12]  - runif(length(rbp_brakes), 2.0, 3.0)
    
    # Compute Z-score activity per layer (inverting RBP brakes inside Layer 3 so score reflects net SASP licensing)
    z_expr <- t(scale(t(expr)))
    z_expr_oriented <- z_expr
    z_expr_oriented[rbp_brakes, ] <- -z_expr_oriented[rbp_brakes, ]
    
    layer_scores <- t(sapply(sasp_4layers, function(g) colMeans(z_expr_oriented[g, , drop = FALSE])))
    rownames(layer_scores) <- c(
      "Layer 1: cGAS-STING / DDR Sensors",
      "Layer 2: NF-kB / C/EBPb / STAT3",
      "Layer 3: mTORC1 & RBP Derepression",
      "Layer 4: Anterograde & Exosome Release"
    )
    
    # ==============================================================================
    # PLOT 1: 4-Layer SASP Cascade Heatmap
    # ==============================================================================
    heat_pal <- colorRampPalette(c("#377EB8", "#FFFFFF", "#E41A1C"))(60)
    
    pheatmap(
      layer_scores,
      main = "4-Layer SASP Regulatory Cascade Across Senescence Progression",
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      color = heat_pal,
      border_color = "#E5E7E9",
      fontsize = 10.5
    )
    
    pheatmap(
      layer_scores,
      main = "4-Layer SASP Regulatory Cascade Across Senescence Progression",
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      color = heat_pal,
      border_color = "#E5E7E9",
      fontsize = 10.5,
      filename = "SASP_Plot1_4Layer_Heatmap.png",
      width = 9.5, height = 4.8
    )
    
    # ==============================================================================
    # PLOT 2: Inter-Layer Crosstalk Correlation Matrix
    # ==============================================================================
    cor_mat <- cor(t(layer_scores), method = "pearson")
    cor_melted <- melt(cor_mat)
    
    p_cormat <- ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
      geom_tile(color = "white", linewidth = 1.2) +
      geom_text(aes(label = sprintf("R = %.2f", value)), color = "black", fontface = "bold", size = 4.2) +
      scale_fill_gradient2(
        low = "#377EB8", mid = "#F8F9F9", high = "#E41A1C",
        midpoint = 0, limit = c(-1, 1), name = "Pearson R"
      ) +
      labs(
        title = "Inter-Layer SASP Crosstalk Coupling Matrix",
        subtitle = "Coordinated feed-forward amplification across sensing, transcription, mRNA stability & secretion",
        x = NULL, y = NULL
      ) +
      theme_minimal(base_size = 11.5) +
      theme(
        axis.text.x = element_text(angle = 25, hjust = 1, face = "bold", color = "#1C2833"),
        axis.text.y = element_text(face = "bold", color = "#1C2833"),
        panel.grid = element_blank()
      )
    
    print(p_cormat)
    ggsave("SASP_Plot2_Crosstalk_Correlation_Matrix.png", plot = p_cormat, width = 8.5, height = 6.5, dpi = 300)
    
    # ==============================================================================
    # PLOT 3: Unmixed 4-Layer SASP Crosstalk Signaling Network
    # ==============================================================================
    crosstalk_edges <- data.frame(
      from = c(
        # Layer 1 Internal & L1 -> L2 Crosstalk
        "CGAS",    "STING1",  "STING1",  "ATM",     "MAPK14",
        # Layer 1 -> L3 Crosstalk (IKK/MK2 inhibit RBP brakes)
        "STING1",  "MAPK14",  "IKBKB",   "MAPKAPK2",
        # Layer 2 Internal & Autocrine Loop
        "RELA",    "CEBPB",   "RELA",    "IL1A",    "IL6",     "JAK1",    "STAT3",   "BRD4",
        # Layer 3 mTOR / RBP -> IL6/IL1A Stabilization & L4 Translation
        "MTOR",    "MTOR",    "ZC3H12A", "ZFP36",   "MTOR",
        # Layer 2/3 -> Layer 4 Anterograde & Exosome Trafficking
        "IL6",     "SEC24D",  "LMAN1",   "RAB1A",   "ARF1",    "RAB6A",   "RAB8A",   "VAMP3",
        "IL6",     "RAB27A"
      ),
      to = c(
        "STING1",  "TBK1",    "RELA",    "RELA",    "MAPKAPK2",
        "IKBKB",   "MTOR",    "ZC3H12A", "ZFP36",
        "IL6",     "IL6",     "CXCL8",   "RELA",    "JAK1",    "STAT3",   "IL6",     "RELA",
        "IL1A",    "MAPKAPK2","IL6",     "CXCL8",   "SEC24D",
        "SEC24D",  "LMAN1",   "RAB1A",   "ARF1",    "RAB6A",   "RAB8A",   "VAMP3",   "SNAP23",
        "RAB27A",  "PDCD6IP"
      )
    )
    
    g <- graph_from_data_frame(crosstalk_edges, directed = TRUE)
    
    # Bright distinct colors per SASP regulatory layer
    V(g)$color <- ifelse(V(g)$name %in% sasp_4layers$L1_Sensors_cGAS_STING_DDR, "#377EB8",        # Bright Royal Blue (Layer 1)
                         ifelse(V(g)$name %in% sasp_4layers$L2_Transcription_NFkB_STAT3, "#E41A1C",      # Bright Crimson Red (Layer 2)
                                ifelse(V(g)$name %in% sasp_4layers$L3_PostTranscriptional_mTOR_RBP, "#FF7F00",  # Bright Amber Orange (Layer 3)
                                       "#4DAF4A")))                                                              # Bright Emerald Green (Layer 4)
    
    set.seed(77)
    kk_layout <- layout_with_kk(g)
    
    draw_4layer_network <- function() {
      par(mar = c(1, 1, 3, 1))
      plot(
        g,
        layout = kk_layout,
        vertex.size = 19,
        vertex.frame.color = "#1C2833",
        vertex.label.cex = 0.68,
        vertex.label.family = "sans",
        vertex.label.font = 2,
        vertex.label.color = "white",
        edge.width = 1.6,
        edge.arrow.size = 0.42,
        edge.color = "#7F8C8D",
        main = "4-Layer SASP Regulatory Crosstalk Network"
      )
      legend(
        "bottomleft",
        legend = c(
          "Layer 1: cGAS-STING & DDR Sensors",
          "Layer 2: NF-kB / C/EBPb / JAK-STAT3",
          "Layer 3: mTORC1 & RBP Stability Switch",
          "Layer 4: Anterograde & Exosome Trafficking"
        ),
        pch = 21,
        pt.bg = c("#377EB8", "#E41A1C", "#FF7F00", "#4DAF4A"),
        col = "#1C2833",
        pt.cex = 1.9,
        cex = 0.8,
        bty = "n"
      )
    }
    
    draw_4layer_network()
    
    png("SASP_Plot3_4Layer_Crosstalk_Network.png", width = 3000, height = 2200, res = 300)
    draw_4layer_network()
    dev.off()
    
    # ==============================================================================
    # PLOT 4: Crosstalk Bridge Centrality Lollipop Plot (Top Hub Regulators)
    # ==============================================================================
    bet_vals <- betweenness(g, directed = FALSE, normalized = TRUE)
    hub_df <- data.frame(
      Gene = names(bet_vals),
      Centrality = as.numeric(bet_vals)
    )
    
    hub_df$Layer <- ifelse(hub_df$Gene %in% sasp_4layers$L1_Sensors_cGAS_STING_DDR, "Layer 1: Sensors",
                           ifelse(hub_df$Gene %in% sasp_4layers$L2_Transcription_NFkB_STAT3, "Layer 2: Transcription",
                                  ifelse(hub_df$Gene %in% sasp_4layers$L3_PostTranscriptional_mTOR_RBP, "Layer 3: mTOR / RBP",
                                         "Layer 4: Trafficking")))
    
    hub_df <- hub_df[order(hub_df$Centrality), ]
    top_hubs <- tail(hub_df, 14)
    top_hubs$Gene <- factor(top_hubs$Gene, levels = top_hubs$Gene)
    
    p_hubs <- ggplot(top_hubs, aes(x = Centrality, y = Gene, color = Layer)) +
      geom_segment(aes(x = 0, xend = Centrality, y = Gene, yend = Gene), linewidth = 1.2) +
      geom_point(size = 5) +
      scale_color_manual(values = c(
        "Layer 1: Sensors"       = "#377EB8",
        "Layer 2: Transcription" = "#E41A1C",
        "Layer 3: mTOR / RBP"    = "#FF7F00",
        "Layer 4: Trafficking"   = "#4DAF4A"
      )) +
      labs(
        title = "Key Inter-Layer Bridge Hubs Connecting the 4 SASP Modules",
        subtitle = "Normalized network betweenness centrality across sensing, transcription, stability & secretion",
        x = "Betweenness Centrality (Crosstalk Bridge Score)",
        y = NULL,
        color = "SASP Layer"
      ) +
      theme_classic(base_size = 12) +
      theme(
        axis.text.y = element_text(face = "bold", color = "#1C2833", size = 10.5),
        legend.position = "bottom"
      )
    
    print(p_hubs)
    ggsave("SASP_Plot4_Hub_Centrality_Lollipop.png", plot = p_hubs, width = 9, height = 6, dpi = 300)
    
    cat("\nDone! All 4 SASP crosstalk plots saved to:", getwd(), "\n")
    # ==============================================================================
    # RStudio Pipeline: Pro-Inflammatory IL-6 / IL-8 SASP Sub-Type Analysis
    # (Feed-Forward Loop, Post-Transcriptional Licensing & Anterograde Secretion)
    # ==============================================================================
    
    library(igraph)
    library(ggplot2)
    library(pheatmap)
    
    # 1. Define Sub-Type Specific SASP Gene Modules
    il6_il8_sasp_modules <- list(
      # Contrasting Early / Immunosuppressive SASP (Suppressed during mature inflammatory switch)
      Early_Notch1_TGFb_SASP = c(
        "NOTCH1", "JAG1", "HES1", "TGFB1", "TGFB2", "SERPINE1", "COL1A1", "ACTA2"
      ),
      # Upstream Amplifier: IL-1a / cGAS-STING / NF-kB / C/EBPb / BRD4
      IL1a_NFkB_CEBPb_Driver = c(
        "IL1A", "IL1R1", "IRAK1", "CGAS", "STING1", "RELA", "NFKB1", "CEBPB", "BRD4"
      ),
      # Core Pro-Inflammatory IL-6 / IL-8 Cytokine & Chemokine Output
      Core_IL6_IL8_Effectors = c(
        "IL6", "IL6R", "IL6ST", "CXCL8", "CXCR1", "CXCR2", "CXCL1", "CXCL2", "CCL2", "STAT3"
      ),
      # Anterograde Vesicular Transport & Exocytosis of IL-6 / IL-8
      IL6_IL8_Anterograde_Export = c(
        "SEC24D", "SAR1A", "LMAN1", "RAB1A", "ARF1", "RAB6A", "RAB8A", "VAMP3", "SNAP23", "STX4"
      )
    )
    
    # 2. Simulate Expression Across Proliferating, Early TGF-b SASP, and Mature IL-6/IL-8 SASP
    set.seed(303)
    all_genes <- unique(unlist(il6_il8_sasp_modules))
    samples <- c(
      paste0("Proliferating_", 1:4),
      paste0("Early_TGFb_SASP_", 1:4),
      paste0("Mature_IL6_IL8_SASP_", 1:4)
    )
    
    expr <- matrix(
      rnorm(length(all_genes) * length(samples), mean = 6.0, sd = 0.55),
      nrow = length(all_genes),
      dimnames = list(all_genes, samples)
    )
    
    notch_genes <- il6_il8_sasp_modules$Early_Notch1_TGFb_SASP
    inflam_genes <- setdiff(all_genes, notch_genes)
    
    # Early SASP: Notch1/TGF-b peaks; IL-6/IL-8 axis is still low
    expr[notch_genes, 5:8]  <- expr[notch_genes, 5:8]  + runif(length(notch_genes), 2.0, 3.0)
    expr[inflam_genes, 5:8] <- expr[inflam_genes, 5:8] + runif(length(inflam_genes), 0.3, 0.8)
    
    # Mature Pro-Inflammatory SASP: Notch1 drops; IL-6/IL-8 & Anterograde Secretion surge
    expr[notch_genes, 9:12]  <- expr[notch_genes, 9:12]  - runif(length(notch_genes), 1.0, 1.8)
    expr[inflam_genes, 9:12] <- expr[inflam_genes, 9:12] + runif(length(inflam_genes), 2.4, 3.6)
    
    z_expr <- t(scale(t(expr)))
    module_scores <- t(sapply(il6_il8_sasp_modules, function(g) colMeans(z_expr[g, , drop = FALSE])))
    rownames(module_scores) <- c(
      "Early Notch1 / TGF-b SASP (Suppressed)",
      "IL-1a / NF-kB / C/EBPb Master Driver",
      "Core IL-6 / IL-8 / STAT3 Effectors",
      "COPII / Golgi / SNARE Vesicle Export"
    )
    
    # ==============================================================================
    # PLOT 1: Heatmap of the SASP Sub-Type Switch (Bright Royal Blue to Crimson Red)
    # ==============================================================================
    bright_pal <- colorRampPalette(c("#377EB8", "#FFFFFF", "#E41A1C"))(60)
    
    pheatmap(
      module_scores,
      main = "Temporal Switch to Pro-Inflammatory IL-6 / IL-8 SASP",
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      color = bright_pal,
      border_color = "#E5E7E9",
      fontsize = 10.5
    )
    
    pheatmap(
      module_scores,
      main = "Temporal Switch to Pro-Inflammatory IL-6 / IL-8 SASP",
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      color = bright_pal,
      border_color = "#E5E7E9",
      fontsize = 10.5,
      filename = "IL6_IL8_SASP_Plot1_SubType_Switch.png",
      width = 9.5, height = 4.8
    )
    
    # ==============================================================================
    # PLOT 2: Unmixed Pro-Inflammatory IL-6/IL-8 SASP Feed-Forward & Trafficking Graph
    # ==============================================================================
    edges <- data.frame(
      from = c(
        # Upstream IL-1a / STING -> NF-kB & C/EBPb
        "IL1A",   "IL1R1",  "IRAK1",  "STING1", "BRD4",
        # Transcription -> IL-6 & IL-8
        "RELA",   "CEBPB",  "RELA",   "CEBPB",  "RELA",
        # Autocrine / Paracrine Receptor & STAT3 Loop
        "IL6",    "IL6R",   "IL6ST",  "STAT3",  "CXCL8",  "CXCL8",
        # Anterograde Vesicular Secretion of IL-6 & IL-8
        "IL6",    "CXCL8",  "SEC24D", "LMAN1",  "RAB1A",  "ARF1",   "RAB8A",  "VAMP3"
      ),
      to = c(
        "IL1R1",  "IRAK1",  "RELA",   "RELA",   "RELA",
        "IL6",    "IL6",    "CXCL8",  "CXCL8",  "IL1A",
        "IL6R",   "IL6ST",  "STAT3",  "IL6",    "CXCR1",  "CXCR2",
        "SEC24D", "SEC24D", "LMAN1",  "RAB1A",  "ARF1",   "RAB8A",  "VAMP3",  "SNAP23"
      )
    )
    
    g <- graph_from_data_frame(edges, directed = TRUE)
    
    V(g)$color <- ifelse(V(g)$name %in% il6_il8_sasp_modules$IL1a_NFkB_CEBPb_Driver, "#377EB8",     # Royal Blue: Upstream Drivers
                         ifelse(V(g)$name %in% il6_il8_sasp_modules$Core_IL6_IL8_Effectors, "#E41A1C",     # Crimson Red: Core IL-6/IL-8 Axis
                                "#FF7F00"))                                                                  # Bright Amber: Anterograde Export
    
    set.seed(55)
    kk_layout <- layout_with_kk(g)
    
    draw_il6_il8_network <- function() {
      par(mar = c(1, 1, 3, 1))
      plot(
        g,
        layout = kk_layout,
        vertex.size = 20,
        vertex.frame.color = "#1C2833",
        vertex.label.cex = 0.70,
        vertex.label.family = "sans",
        vertex.label.font = 2,
        vertex.label.color = "white",
        edge.width = 1.7,
        edge.arrow.size = 0.45,
        edge.color = "#7F8C8D",
        main = "Pro-Inflammatory IL-6 / IL-8 SASP Feed-Forward & Trafficking Circuit"
      )
      legend(
        "bottomleft",
        legend = c(
          "IL-1a / STING / NF-kB / C/EBPb Drivers",
          "Core IL-6 / IL-8 / STAT3 Effectors",
          "COPII / Golgi / SNARE Vesicular Export"
        ),
        pch = 21,
        pt.bg = c("#377EB8", "#E41A1C", "#FF7F00"),
        col = "#1C2833",
        pt.cex = 1.9,
        cex = 0.82,
        bty = "n"
      )
    }
    
    draw_il6_il8_network()
    
    png("IL6_IL8_SASP_Plot2_FeedForward_Network.png", width = 2900, height = 2150, res = 300)
    draw_il6_il8_network()
    dev.off()
    
    # ==============================================================================
    # PLOT 3: Gene-Level Lollipop Plot (Mature IL-6/IL-8 SASP vs. Proliferating)
    # ==============================================================================
    gene_df <- do.call(rbind, lapply(all_genes, function(g) {
      prolif <- expr[g, 1:4]
      mature <- expr[g, 9:12]
      lfc    <- mean(mature) - mean(prolif)
      pval   <- t.test(mature, prolif)$p.value
      
      mod <- names(il6_il8_sasp_modules)[sapply(il6_il8_sasp_modules, function(x) g %in% x)][1]
      data.frame(
        Gene = g,
        Module = mod,
        Log2FC = lfc,
        NegLog10P = -log10(pval),
        Direction = ifelse(lfc >= 0, "Activated in IL-6/IL-8 SASP", "Suppressed (Notch1/TGF-b Switch)")
      )
    }))
    
    gene_df <- gene_df[order(gene_df$Log2FC), ]
    gene_df$Gene <- factor(gene_df$Gene, levels = gene_df$Gene)
    
    p_lolli <- ggplot(gene_df, aes(x = Log2FC, y = Gene, color = Direction)) +
      geom_vline(xintercept = 0, color = "#1C2833", linewidth = 0.6) +
      geom_segment(aes(x = 0, xend = Log2FC, y = Gene, yend = Gene), linewidth = 1.0) +
      geom_point(aes(size = NegLog10P)) +
      scale_color_manual(values = c(
        "Activated in IL-6/IL-8 SASP"    = "#E41A1C",  # Bright Crimson Red
        "Suppressed (Notch1/TGF-b Switch)" = "#377EB8"   # Bright Royal Blue
      )) +
      scale_size_continuous(range = c(2.8, 6.2), name = "-log10(p-value)") +
      labs(
        title = "Gene Signature of the Pro-Inflammatory IL-6 / IL-8 SASP",
        subtitle = "Log2 Fold Change (Mature IL-6/IL-8 SASP vs. Proliferating) showing Notch1 suppression & IL-6/IL-8 surge",
        x = "Log2 Fold Change",
        y = NULL,
        color = "SASP Sub-Type Shift"
      ) +
      theme_classic(base_size = 11.5) +
      theme(
        axis.text.y = element_text(face = "bold", color = "#1C2833", size = 9),
        legend.position = "bottom"
      )
    
    print(p_lolli)
    ggsave("IL6_IL8_SASP_Plot3_Effector_Lollipop.png", plot = p_lolli, width = 9.5, height = 8.5, dpi = 300)
    
    cat("\nDone! Saved all 3 IL-6/IL-8 SASP plots to:", getwd(), "\n")
    