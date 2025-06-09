#### v1.0支持条件 ####
#仅支持每组单样本的常规qpcr，样本的三次重复单独计算以保留误差棒

# ==== 加载包 ====
library(shiny)
library(plotly)
library(dplyr)
library(DT)
library(RColorBrewer)

# ==== 构建 384 孔板坐标 ====
rows <- LETTERS[1:16]  # A-P
cols <- 1:24
plate_layout <- expand.grid(Row = rows, Col = cols)
plate_layout <- plate_layout %>%
  mutate(Pos = paste0(Row, Col), Sample = NA, Gene = NA)

# ==== UI 界面 ====
ui <- fluidPage(
  titlePanel("384布板工具 Pos-Sample-Gene"),
  sidebarLayout(
    sidebarPanel(
      textInput("sample", "样本名", value = ""),
      selectizeInput("gene", "基因名", choices = NULL, selected = "", options = list(
        create = TRUE,
        placeholder = '输入或选择基因名'
      )),
      actionButton("apply", "应用到选中孔位"),
      actionButton("clear", "清除所有标注"),
      br(), br(),
      downloadButton("downloadData", "下载结果为CSV"),
      uiOutput("sample_legend"),
      uiOutput("gene_legend")
    ),
    mainPanel(
      plotlyOutput("platePlot", height = "800px"),
      DTOutput("preview")
    )
  )
)

# ==== Server ====
`%||%` <- function(a, b) if (!is.null(a)) a else b

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
      
      if (nzchar(input$sample)) {
        updated$Sample[idx] <- input$sample
      }
      if (nzchar(input$gene)) {
        updated$Gene[idx] <- input$gene
      }
      
      plate_data(updated)
      
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
    
    fill_palette <- brewer.pal(8, "Pastel1")
    border_palette <- brewer.pal(8, "Dark2")
    
    unique_samples <- unique(na.omit(pdata$Sample))
    sample_colors <- setNames(rep("white", length(pdata$Sample)), pdata$Sample)
    if (length(unique_samples) > 0) {
      sample_colors <- setNames(fill_palette[seq_along(unique_samples)], unique_samples)
    }
    pdata$FillColor <- ifelse(is.na(pdata$Sample), "white", sample_colors[pdata$Sample])
    
    unique_genes <- unique(na.omit(pdata$Gene))
    gene_colors <- setNames(rep("black", length(pdata$Gene)), pdata$Gene)
    if (length(unique_genes) > 0) {
      gene_colors <- setNames(border_palette[seq_along(unique_genes)], unique_genes)
    }
    pdata$BorderColor <- ifelse(is.na(pdata$Gene), "black", gene_colors[pdata$Gene])
    
    # ⚠️ 信息不完整标识
    pdata$IsIncomplete <- (is.na(pdata$Sample) != is.na(pdata$Gene))
    pdata$DisplayText <- ifelse(
      pdata$IsIncomplete,
      paste0(pdata$Pos, " ⚠️"),
      pdata$Pos
    )
    pdata$HoverText <- ifelse(
      pdata$IsIncomplete,
      paste("⚠️ 信息不完整\nPos:", pdata$Pos,
            "\nSample:", pdata$Sample,
            "\nGene:", pdata$Gene),
      paste("Pos:", pdata$Pos,
            "\nSample:", pdata$Sample,
            "\nGene:", pdata$Gene)
    )
    
    plot_ly(
      data = pdata,
      x = ~Col,
      y = ~Row,
      type = "scatter",
      mode = "markers+text",
      text = ~DisplayText,
      textposition = "middle center",
      textfont = list(size = 12),
      hovertext = ~HoverText,
      hoverinfo = "text",
      marker = list(
        size = 30,
        symbol = "circle",
        color = pdata$FillColor,
        line = list(width = ifelse(pdata$Selected, 4, 2), color = pdata$BorderColor)
      ),
      key = ~Pos
    ) %>%
      layout(
        paper_bgcolor = "white",
        plot_bgcolor = "white",
        newshape = list(line = list(color = "black")),
        uirevision = "fixed",
        dragmode = "select",
        xaxis = list(title = "", side = "top", range = c(0.5, 24.5), tickvals = 1:24, fixedrange = TRUE),
        yaxis = list(title = "", categoryorder = "array", categoryarray = rev(LETTERS[1:16]), fixedrange = TRUE),
        title = "框选孔位 → 填写样本/基因 → 应用",
        showlegend = FALSE
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
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("plate_annotation_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(plate_data(), file, row.names = FALSE)
    }
  )
}

# ==== 启动 App ====
shinyApp(ui, server)


# ==== 后续整合分析脚本（从 App 输出 CSV + Cp 原始数据进行 2^-ddCt 分析） ====

# ==== 参数设置 ====
ref_gene <- "actin"           # 内参基因
ref_sample <- "HEP0"          # 对照样本名
app_csv_path <- "C:/Users/c/Desktop/plate_annotation_2025-04-15.csv"  # App导出CSV路径
qpcr_txt_path <- "C:/Users/c/Desktop/250402 heplb_nc+day1-3.txt"      # qPCR原始TXT路径

# ==== 加载必要包 ====
library(readr)
library(dplyr)
library(readxl)
library(writexl)
library(stringr)

# ==== 读取数据 ====
raw_data <- read.delim(qpcr_txt_path, skip = 1, header = TRUE, stringsAsFactors = FALSE)
raw_data <- raw_data[, c("Pos", "Cp")]

annotation <- read.csv(app_csv_path)

# ==== 合并标注信息 ====
data <- raw_data %>%
  left_join(annotation[, c("Pos", "Sample", "Gene")], by = "Pos") %>%
  filter(!is.na(Sample) & !is.na(Gene))

# ==== 赋予 replicate 编号 ====
data <- data %>%
  group_by(Sample, Gene) %>%
  mutate(Replicate = row_number()) %>%
  ungroup()

# ==== ΔΔCt 分析准备 ====
samples <- unique(data$Sample)
genes <- unique(data$Gene)

results_matrix <- matrix(NA, nrow = length(samples) * length(genes), ncol = 5)
dimnames(results_matrix) <- list(NULL, c("Sample", "Gene", "2^-ddCt_1", "2^-ddCt_2", "2^-ddCt_3"))
counter <- 0

# 每个样品的内参平均 Ct
ref_gene_mean <- data %>%
  filter(Gene == ref_gene) %>%
  group_by(Sample) %>%
  summarise(ref_gene_mean = mean(Cp, na.rm = TRUE)) %>%
  ungroup()

# 对照样本中每个基因的平均 Ct
ref_sample_mean <- data %>%
  filter(Sample == ref_sample) %>%
  group_by(Gene) %>%
  summarise(ref_sample_mean = mean(Cp, na.rm = TRUE)) %>%
  ungroup()

# ==== ΔΔCt 计算 ====
for (sample in samples) {
  for (gene in genes) {
    target_ct <- data %>%
      filter(Sample == sample & Gene == gene) %>%
      pull(Cp)
    
    delta_ct <- target_ct - ref_gene_mean$ref_gene_mean[ref_gene_mean$Sample == sample]
    ref_ct <- ref_sample_mean %>% filter(Gene == gene) %>% pull(ref_sample_mean)
    ref_ref_ct <- ref_sample_mean %>% filter(Gene == ref_gene) %>% pull(ref_sample_mean)
    
    delta_delta_ct <- delta_ct - (ref_ct - ref_ref_ct)
    two_to_minus_delta_delta_ct <- 2^-delta_delta_ct
    
    counter <- counter + 1
    results_matrix[counter, "Sample"] <- sample
    results_matrix[counter, "Gene"] <- gene
    results_matrix[counter, "2^-ddCt_1"] <- two_to_minus_delta_delta_ct[1]
    results_matrix[counter, "2^-ddCt_2"] <- two_to_minus_delta_delta_ct[2]
    results_matrix[counter, "2^-ddCt_3"] <- two_to_minus_delta_delta_ct[3]
  }
}

# ==== 输出分析结果 ====
write_xlsx(as.data.frame(results_matrix), path = "C:/Users/c/Desktop/analysised.xlsx")
