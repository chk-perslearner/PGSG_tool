#### v1.1支持条件 ####
# 支持每组单样本或多样本的常规qpcr，单样本的三次重复强制单独计算以保留误差棒，多样本在PGSG界面可选将每个样本以三次重复或均值进行计算

library(shiny)
library(plotly)
library(dplyr)
library(DT)
library(RColorBrewer)

# ==== 构建 384 孔板坐标 ====
rows <- LETTERS[1:16]
cols <- 1:24
plate_layout <- expand.grid(Row = rows, Col = cols)
plate_layout <- plate_layout %>%
  mutate(Pos = paste0(Row, Col), Group = NA, Sample = NA, Gene = NA)

# ==== UI ====
ui <- fluidPage(
  titlePanel("PSG · 384布板工具 Pos-Sample-Gene"),
  sidebarLayout(
    sidebarPanel(
      HTML("<b>🟦 配色说明：</b><br>
           <span style='color:#333;'>孔的填充颜色：代表样本（Sample）<br>
           孔的边框颜色：代表基因（Gene）<br>
           孔中文字颜色：代表组别（Group）</span><br><br>"),
      textInput("group", "组别名", value = ""),
      textInput("sample", "样本名", value = ""),
      selectizeInput("gene", "基因名", choices = NULL, selected = "", options = list(
        create = TRUE,
        placeholder = '输入或选择基因名'
      )),
      checkboxInput("use_replicates", "保留多样本组的重复Ct值", FALSE),
      textInput("ref_group", "对照组组别名（Group）", value = ""),
      textInput("ref_gene", "内参基因名", value = ""),
      actionButton("apply", "应用到选中孔位"),
      actionButton("clear", "清除所有标注"),
      br(), br(),
      downloadButton("downloadData", "下载标注为CSV"),
      downloadButton("downloadParams", "下载分析参数"),
      uiOutput("sample_legend"),
      uiOutput("gene_legend"),
      uiOutput("group_legend")
    ),
    mainPanel(
      plotlyOutput("platePlot", height = "800px"),
      DTOutput("preview")
    )
  )
)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ==== SERVER ====
server <- function(input, output, session) {
  plate_data <- reactiveVal(plate_layout)
  selected_points <- reactiveVal(data.frame())
  
  observe({
    gene_choices <- plate_data() %>%
      filter(!is.na(Gene)) %>%
      pull(Gene) %>%
      unique() %>%
      sort()
    updateSelectizeInput(session, "gene", choices = c("", gene_choices), server = TRUE)
  })
  
  observe({
    select_evt <- event_data("plotly_selected")
    click_evt  <- event_data("plotly_click")
    
    if (!is.null(select_evt) && !is.null(select_evt$key)) {
      selected_points(select_evt)
      return()
    }
    
    if (!is.null(click_evt) && !is.null(click_evt$key)) {
      prev <- selected_points()
      k <- click_evt$key
      if (k %in% prev$key) {
        updated <- prev %>% filter(key != k)
      } else {
        updated <- bind_rows(prev, click_evt) %>% distinct(key, .keep_all = TRUE)
      }
      selected_points(updated)
      return()
    }
    
    selected_points(data.frame())
  })
  
  observeEvent(input$apply, {
    brush <- selected_points()
    if (nrow(brush) > 0) {
      selected_pos <- unique(brush$key)
      updated <- plate_data()
      idx <- which(updated$Pos %in% selected_pos)
      
      if (nzchar(input$group))  updated$Group[idx]  <- input$group
      if (nzchar(input$sample)) updated$Sample[idx] <- input$sample
      if (nzchar(input$gene))   updated$Gene[idx]   <- input$gene
      
      plate_data(updated)
      updateTextInput(session, "group", value = "")
      updateTextInput(session, "sample", value = "")
      updateSelectizeInput(session, "gene", selected = "")
    }
  })
  
  observeEvent(input$clear, {
    plate_data(plate_layout)
    selected_points(data.frame())
  })
  
  output$platePlot <- renderPlotly({
    pdata <- plate_data()
    selected <- selected_points()
    pdata$Selected <- pdata$Pos %in% selected$key
    
    fill_palette   <- brewer.pal(8, "Pastel1")
    border_palette <- brewer.pal(8, "Dark2")
    text_palette   <- brewer.pal(8, "Set1")
    
    # Sample填充
    unique_samples <- unique(na.omit(pdata$Sample))
    sample_colors <- setNames(fill_palette[seq_along(unique_samples)], unique_samples)
    pdata$FillColor <- ifelse(is.na(pdata$Sample), "white", sample_colors[pdata$Sample])
    
    # Gene边框
    unique_genes <- unique(na.omit(pdata$Gene))
    gene_colors <- setNames(border_palette[seq_along(unique_genes)], unique_genes)
    pdata$BorderColor <- ifelse(is.na(pdata$Gene), "black", gene_colors[pdata$Gene])
    
    # Group文字
    unique_groups <- unique(na.omit(pdata$Group))
    group_colors <- setNames(text_palette[seq_along(unique_groups)], unique_groups)
    pdata$FontColor <- ifelse(is.na(pdata$Group), "black", group_colors[pdata$Group])
    
    # 警示标志
    pdata$IsIncomplete <- (is.na(pdata$Sample) != is.na(pdata$Gene))
    pdata$DisplayText <- ifelse(pdata$IsIncomplete, paste0(pdata$Pos, " ⚠️"), pdata$Pos)
    pdata$HoverText <- paste("Pos:", pdata$Pos,
                             "\nGroup:", pdata$Group,
                             "\nSample:", pdata$Sample,
                             "\nGene:", pdata$Gene)
    
    plot_ly(
      data = pdata,
      x = ~Col,
      y = ~Row,
      type = "scatter",
      mode = "markers+text",
      text = ~DisplayText,
      textposition = "middle center",
      textfont = list(color = ~FontColor, size = 12),
      hovertext = ~HoverText,
      hoverinfo = "text",
      marker = list(
        size = 30,
        symbol = "circle",
        color = ~FillColor,
        line = list(width = ifelse(pdata$Selected, 4, 2), color = ~BorderColor)
      ),
      key = ~Pos
    ) %>%
      layout(
        dragmode = "select",
        xaxis = list(side = "top", range = c(0.5, 24.5), tickvals = 1:24, fixedrange = TRUE),
        yaxis = list(categoryorder = "array", categoryarray = rev(LETTERS[1:16]), fixedrange = TRUE),
        paper_bgcolor = "white", plot_bgcolor = "white", showlegend = FALSE
      )
  })
  
  output$preview <- renderDT({
    datatable(plate_data() %>% filter(!is.na(Sample) | !is.na(Gene)))
  })
  
  output$sample_legend <- renderUI({
    pdata <- plate_data()
    unique_samples <- unique(na.omit(pdata$Sample))
    if (length(unique_samples) == 0) return(NULL)
    fill_palette <- brewer.pal(8, "Pastel1")
    sample_colors <- setNames(fill_palette[seq_along(unique_samples)], unique_samples)
    
    HTML(paste0(
      "<b>样本颜色：</b><br>",
      paste0(
        "<div style='display:inline-block;margin-right:12px;'>",
        "<span style='display:inline-block;width:12px;height:12px;background:", sample_colors, ";border:1px solid #000;margin-right:4px;'></span>",
        names(sample_colors),
        "</div>", collapse = ""
      )
    ))
  })
  
  output$gene_legend <- renderUI({
    pdata <- plate_data()
    unique_genes <- unique(na.omit(pdata$Gene))
    if (length(unique_genes) == 0) return(NULL)
    border_palette <- brewer.pal(8, "Dark2")
    gene_colors <- setNames(border_palette[seq_along(unique_genes)], unique_genes)
    
    HTML(paste0(
      "<br><b>基因边框颜色：</b><br>",
      paste0(
        "<div style='display:inline-block;margin-right:12px;'>",
        "<span style='display:inline-block;width:12px;height:12px;background:white;border:2px solid ", gene_colors, ";margin-right:4px;'></span>",
        names(gene_colors),
        "</div>", collapse = ""
      )
    ))
  })
  
  output$group_legend <- renderUI({
    pdata <- plate_data()
    unique_groups <- unique(na.omit(pdata$Group))
    if (length(unique_groups) == 0) return(NULL)
    text_palette <- brewer.pal(8, "Set1")
    group_colors <- setNames(text_palette[seq_along(unique_groups)], unique_groups)
    
    HTML(paste0(
      "<br><b>组别文字颜色：</b><br>",
      paste0(
        "<div style='display:inline-block;margin-right:12px;'>",
        "<span style='display:inline-block;width:12px;height:12px;background:", group_colors, ";color:", group_colors, ";border:1px solid #000;margin-right:4px;'></span>",
        names(group_colors),
        "</div>", collapse = ""
      )
    ))
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("PGSG_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(plate_data(), file, row.names = FALSE)
    }
  )
  
  output$downloadParams <- downloadHandler(
    filename = function() {
      paste0("PGSG_parameters_", Sys.Date(), ".txt")
    },
    content = function(file) {
      lines <- c(
        paste0("ref_gene=", input$ref_gene),
        paste0("ref_group=", input$ref_group),
        paste0("use_replicates_for_multi_sample_group=", input$use_replicates)
      )
      writeLines(lines, file)
    }
  )
}

# ==== 启动 App ====
shinyApp(ui, server)

















# ==== 设置路径 ====
param_path    <- "C:/Users/c/Desktop/PGSG_parameters_2025-04-17.txt"  # 参数文件
app_csv_path  <- "C:/Users/c/Desktop/PGSG_2025-04-17.csv"             # PSG标注文件
qpcr_txt_path <- "C:/Users/c/Desktop/250402 heplb_nc+day1-3.txt"      # 原始数据文件
output_path   <- "C:/Users/c/Desktop/analysised.xlsx"                 # 输出路径

# ==== 加载包 ====
library(readr)
library(dplyr)
library(writexl)
library(stringr)

# ==== 读取参数 ====
params <- readLines(param_path)
get_param <- function(key) {
  val <- grep(paste0("^", key, "="), params, value = TRUE)
  if (length(val) == 0) return("")
  sub(paste0("^", key, "="), "", val)
}
ref_gene <- get_param("ref_gene")
ref_group <- get_param("ref_group")
use_replicates_raw <- get_param("use_replicates_for_multi_sample_group")
use_replicates <- tolower(use_replicates_raw) == "true"

cat("🔧 参数读取完成：\n")
cat("ref_gene =", ref_gene, "\n")
cat("ref_group =", ref_group, "\n")
cat("use_replicates_for_multi_sample_group =", use_replicates, "\n\n")

# ==== 读取数据 ====
cat("📥 正在读取数据...\n")
raw_data <- read.delim(qpcr_txt_path, skip = 1)
raw_data <- raw_data[, c("Pos", "Cp")]
annotation <- read.csv(app_csv_path)

# ==== 合并并预处理 ====
cat("🔗 合并 PSG 标注...\n")
merged_data <- raw_data %>%
  left_join(annotation[, c("Pos", "Group", "Sample", "Gene")], by = "Pos") %>%
  filter(!is.na(Group), !is.na(Sample), !is.na(Gene), !is.na(Cp))

# ==== 判断组别样本数 ====
sample_counts <- merged_data %>%
  group_by(Group) %>%
  summarise(n_sample = n_distinct(Sample), .groups = "drop")

merged_data <- merged_data %>%
  left_join(sample_counts, by = "Group")

# ==== ΔCt计算逻辑 ====
cat("📊 正在计算 ΔCt...\n")
# 计算每个样本的内参Ct（平均）
ref_gene_mean <- merged_data %>%
  filter(Gene == ref_gene) %>%
  group_by(Sample, Group) %>%
  summarise(ref_Ct = mean(Cp), .groups = "drop")

# 合并内参Ct
data <- merged_data %>%
  left_join(ref_gene_mean, by = c("Sample", "Group")) %>%
  mutate(delta_Ct = Cp - ref_Ct)

# ==== 计算对照组 ΔCt 均值 ====
ref_dCt <- data %>%
  filter(Group == ref_group) %>%
  group_by(Gene) %>%
  summarise(ref_dCt = mean(delta_Ct), .groups = "drop")

# ==== ΔΔCt 和 2^-ΔΔCt ====
cat("📈 正在计算 ΔΔCt 和 2^-ΔΔCt...\n")
data <- data %>%
  left_join(ref_dCt, by = "Gene") %>%
  mutate(delta_delta_Ct = delta_Ct - ref_dCt,
         `2^-ΔΔCt` = 2^-delta_delta_Ct)

# ==== 如果为多样本组且设置为不保留重复，则对样本Ct取平均 ====
if (!use_replicates) {
  data <- data %>%
    mutate(replicate = row_number()) %>%
    group_by(Group, Sample, Gene) %>%
    summarise(Ct_1 = Cp[1], Ct_2 = Cp[2], Ct_3 = Cp[3],
              `2^-ΔΔCt_1` = `2^-ΔΔCt`[1],
              `2^-ΔΔCt_2` = `2^-ΔΔCt`[2],
              `2^-ΔΔCt_3` = `2^-ΔΔCt`[3],
              .groups = "drop")
} else {
  data <- data %>%
    mutate(replicate = rep(1:3, length.out = nrow(data))) %>%
    select(Group, Sample, Gene, replicate, Cp, `2^-ΔΔCt`) %>%
    arrange(Group, Sample, Gene, replicate)
}

# ==== 输出 Excel ====
cat("💾 写入 Excel 文件：", output_path, "\n")
write_xlsx(data, path = output_path)
cat("✅ 分析完成。\n")
